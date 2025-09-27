nextflow.enable.dsl=2

// Import modules
include { CELLBENDER; CELLBENDER_GPU; CELLBENDER_H5_CONVERT } from './modules/cellbender/cellbender'
include { GENERATE_REPORTS; COMBINE_REPORTS; GENERATE_COMBINED_REPORT } from './modules/reports/reports'
include { CREATE_SEURAT } from './modules/seurat/seurat'
include { DROPLETQC } from './modules/dropletqc/dropletqc'
include { SCDBL } from './modules/scdbl/scdbl'

// Default parameters
params.mapping_dirs = "${projectDir}/personal/mapping_dirs.csv"
params.outputDir = "results"
params.gpu = false
params.report = true
params.book = false
params.cellbender = false
params.help = false
params.max_mito = 10.0
params.min_nuclear = 0.4
params.metadata = null // Optional metadata CSV file

// Help message
def helpMessage() {
    log.info """
    ===================================
    scQC-flow
    ===================================
    
    Usage:
    nextflow run main.nf --mapping_dirs <mapping_dirs.csv> --outputDir <output_directory>
    
    Required arguments:
      --mapping_dirs        Path to mapping directories CSV file
      --outputDir           Output directory for results
      
    Optional arguments:
      --cellbender           Run CellBender ambient RNA removal
      --gpu                  Use GPU acceleration for cellbender (requires --cellbender)
      --report               Generate Quarto QC reports for each sample
      --book                 Combine all reports into a single Quarto book
      --max_mito             Maximum mitochondrial percentage threshold (default: 20.0)
      --min_nuclear          Minimum nuclear fraction threshold (default: 0.0)
      --help                 Show this help message
      
    Examples:
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results --report
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results --report --book
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results --cellbender
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results --cellbender --gpu
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results --cellbender --report --book
    """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Main workflow
workflow {
    // Log the pipeline parameters
    log.info """
    ===================================
    scQC-flow pipeline
    ===================================
    Mapping directories: ${params.mapping_dirs}
    Metadata file: ${params.metadata}
    Output directory: ${params.outputDir}
    CellBender: ${params.cellbender}
    GPU acceleration: ${params.gpu}
    Generate reports: ${params.report}
    Generate book: ${params.book}
    QC Thresholds:
      Max mitochondrial %: ${params.max_mito}
      Min nuclear fraction: ${params.min_nuclear}
    ===================================
    """
    
    // Validate required parameters
    if (!params.mapping_dirs) {
        error "Mapping directories file not provided. Please specify --mapping_dirs"
    }
    
    // Read mapping directory information from the CSV file
    Channel
        .fromPath(params.mapping_dirs)
        .splitCsv(header: true)
        .map { row -> tuple(row.samplename, file(row.path)) }
        .set { sampleChannelBase }

    // Prepare script/template files as value channels
    dropletqc_script_path = file("${projectDir}/modules/dropletqc/run_dropletqc.R")
    scdbl_script_path = file("${projectDir}/modules/scdbl/run_scdbl.R")
    report_template_path = file("${projectDir}/modules/reports/seurat_template.qmd")
    combined_template_path = file("${projectDir}/modules/reports/combined_template.qmd")
    book_template_path = file("${projectDir}/modules/reports/book_template/")
    seurat_script_path = file("${projectDir}/modules/seurat/make_seurat.R")

    // Expand sampleChannel to include script inputs using map
    dropletqc_input_ch = sampleChannelBase.map { sampleName, mappingDir -> tuple(sampleName, mappingDir, dropletqc_script_path) }
    scdbl_input_ch = sampleChannelBase.map { sampleName, mappingDir -> tuple(sampleName, mappingDir, scdbl_script_path) }

    // Always run DropletQC nuclear fraction analysis
    dropletqc_results = DROPLETQC(dropletqc_input_ch)

    // Always run scDblFinder doublet detection
    scdbl_results = SCDBL(scdbl_input_ch)

    // Link CellBender and Seurat processes so Seurat waits for CellBender/H5 conversion if requested
    if (params.cellbender) {
        log.info "Running CellBender for all samples"
        if (params.gpu) {
            log.info "GPU acceleration enabled for CellBender"
            cellbender_results = CELLBENDER_GPU(sampleChannelBase)
        } else {
            cellbender_results = CELLBENDER(sampleChannelBase)
        }
        // Run H5 conversion after CellBender (CPU or GPU)
        cellbender_h5_results = CELLBENDER_H5_CONVERT(cellbender_results.cellbender_output)

        // After DropletQC and scDbl, create Seurat objects using CellBender H5
        seurat_input_ch = sampleChannelBase
            .join(dropletqc_results.metrics)
            .join(scdbl_results.metrics)
            .join(cellbender_h5_results.seurat_h5)
            // structure: [sampleName, mappingDir, dropletqc_metrics, scdbl_metrics, cellbender_h5]
            .map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) }

        seurat_input_with_script = seurat_input_ch.map { sampleName, mappingDir, dropletqc, scdbl, h5_path -> tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, params.max_mito, params.min_nuclear, params.metadata, h5_path) }
        seurat_results = CREATE_SEURAT(seurat_input_with_script)
    } else {
        // Use default 10X H5 if CellBender is not run
        h5_path_ch = sampleChannelBase.map { sampleName, mappingDir ->
            def default_h5 = file("${mappingDir}/outs/filtered_feature_bc_matrix.h5")
            tuple(sampleName, default_h5)
        }
        seurat_input_ch = sampleChannelBase
            .join(dropletqc_results.metrics)
            .join(scdbl_results.metrics)
            .join(h5_path_ch)
            // structure: [sampleName, mappingDir, dropletqc_metrics, scdbl_metrics, h5_path]
            .map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) }

        seurat_input_with_script = seurat_input_ch.map { sampleName, mappingDir, dropletqc, scdbl, h5_path -> tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, params.max_mito, params.min_nuclear, params.metadata, h5_path) }
        seurat_results = CREATE_SEURAT(seurat_input_with_script)
    }

    // Prepare GENERATE_REPORTS input channel to use Seurat objects
    // seurat_results emits: (sampleName, pre_rds, post_rds)
    // Need to rejoin with mappingDir from sampleChannelBase for reports
    report_input_ch = sampleChannelBase
        .join(seurat_results)
        .map { sampleName, mappingDir, pre_rds, post_rds -> tuple(sampleName, mappingDir, pre_rds, post_rds, report_template_path, params.max_mito, params.min_nuclear) }

    // Debug: Print contents of report_input_ch

    // Conditionally run Quarto report generation
    if (params.report) {
        

        reports_output = GENERATE_REPORTS(report_input_ch)

        // // Generate combined report for side-by-side comparison
        // log.info "Generating combined QC report for all samples"
        // combined_report_output = GENERATE_COMBINED_REPORT(

        // Conditionally run CellBender if requested
        if (params.cellbender) {
            log.info "Running CellBender for all samples"
            if (params.gpu) {
                log.info "GPU acceleration enabled for CellBender"
                cellbender_results = CELLBENDER_GPU(sampleChannelBase)
            } else {
                cellbender_results = CELLBENDER(sampleChannelBase)
            }
            // Run H5 conversion after CellBender (CPU or GPU)
            cellbender_h5_results = CELLBENDER_H5_CONVERT(cellbender_results.cellbender_output)
        }
        //     report_input_ch.map { it -> it[0] }.collect(),
        //     report_input_ch.map { it -> it[1] }.collect(),
        //     report_input_ch.map { it -> it[2] }.collect(),
        //     report_input_ch.map { it -> it[3] }.collect(),
        //     combined_template_path
        // )

        // Optionally combine all reports into a single book
        if (params.book) {
            log.info "Combining reports into a Quarto book"

            // Collect all HTML reports and QMD sources (extract just the file paths)
            all_html_reports = reports_output.html_report.map { sampleName, htmlFile -> htmlFile }.collect()
            all_qmd_sources = reports_output.qmd_source.map { sampleName, qmdFile -> qmdFile }.collect()

            // Generate combined book
            combined_book = COMBINE_REPORTS(all_html_reports, all_qmd_sources, book_template_path)
        }

        // For backward compatibility, also run the original report if available
        // quarto_sc_report(sampleChannelBase)
    }
}
