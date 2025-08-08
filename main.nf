nextflow.enable.dsl=2

// Import modules
include { CELLBENDER; CELLBENDER_GPU; CELLBENDER_H5_CONVERT } from './modules/cellbender/cellbender'
include { GENERATE_REPORTS; COMBINE_REPORTS; GENERATE_COMBINED_REPORT } from './modules/reports/reports'
include { DROPLETQC } from './modules/dropletqc/dropletqc'
include { SCDBL } from './modules/scdbl/scdbl'

// Default parameters
params.mapping_dirs = "${projectDir}/personal/mapping_dirs.csv"
params.outputDir = "results"
params.gpu = false
params.report = false
params.book = false
params.cellbender = false
params.help = false

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
    Output directory: ${params.outputDir}
    CellBender: ${params.cellbender}
    GPU acceleration: ${params.gpu}
    Generate reports: ${params.report}
    Generate book: ${params.book}
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
        .set { sampleChannel }

    // Conditionally run CellBender on each sample
    if (params.cellbender) {
        if (params.gpu) {
            log.info "Running CellBender with GPU acceleration"
            cellbender_results = CELLBENDER_GPU(sampleChannel)
        } else {
            log.info "Running CellBender with CPU"
            cellbender_results = CELLBENDER(sampleChannel)
        }
        
        // Convert H5 files to Seurat-compatible format
        seurat_h5_results = CELLBENDER_H5_CONVERT(cellbender_results.cellbender_output)
    } else {
        log.info "Skipping CellBender - processing raw Cell Ranger outputs"
    }

    // Always run DropletQC nuclear fraction analysis
    log.info "Running DropletQC nuclear fraction analysis"
    dropletqc_results = DROPLETQC(sampleChannel)

    // Always run scDblFinder doublet detection
    log.info "Running scDblFinder doublet detection"
    scdbl_results = SCDBL(sampleChannel)

    // Conditionally run Quarto report generation
    if (params.report) {
        log.info "Generating Quarto QC reports (including DropletQC and scDblFinder data)"
        
        // Wait for both DropletQC and scDblFinder to complete before generating reports
        combined_qc_data = sampleChannel.join(dropletqc_results.metrics).join(scdbl_results.metrics)
        reports_output = GENERATE_REPORTS(combined_qc_data)
        
        // Generate combined report for side-by-side comparison
        log.info "Generating combined QC report for all samples"
        combined_report_output = GENERATE_COMBINED_REPORT(
            combined_qc_data.map { it -> it[0] }.collect(),
            combined_qc_data.map { it -> it[1] }.collect(),
            combined_qc_data.map { it -> it[2] }.collect(),
            combined_qc_data.map { it -> it[3] }.collect()
        )
        
        // Optionally combine all reports into a single book
        if (params.book) {
            log.info "Combining reports into a Quarto book"
            
            // Collect all HTML reports and QMD sources (extract just the file paths)
            all_html_reports = reports_output.html_report.map { sampleName, htmlFile -> htmlFile }.collect()
            all_qmd_sources = reports_output.qmd_source.map { sampleName, qmdFile -> qmdFile }.collect()
            
            // Generate combined book
            combined_book = COMBINE_REPORTS(all_html_reports, all_qmd_sources)
        }
        
        // For backward compatibility, also run the original report if available
        // quarto_sc_report(sampleChannel)
    }
}
