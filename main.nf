nextflow.enable.dsl=2

// Import workflows
include { STANDARD_WORKFLOW } from './workflows/workflows'
include { MULTIOME_WORKFLOW } from './workflows/workflows'
include { REPORTING } from './workflows/workflows'

// =============================================================================
// PARAMETERS
// =============================================================================
params.mapping_dirs = "${projectDir}/personal/mapping_dirs.csv"
params.outputDir = "results"
params.gpu = false
params.report = true
params.book = false
params.cellbender = false
params.multiome = false
params.help = false
params.max_mito = 20.0
params.min_nuclear = 0.4
params.metadata = null

// =============================================================================
// HELP MESSAGE
// =============================================================================
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
      --multiome             Run multiome workflow for 10x Multiome data
      --report               Generate Quarto QC reports for each sample
      --book                 Combine all reports into a single Quarto book
      --max_mito             Maximum mitochondrial percentage threshold (default: 20.0)
      --min_nuclear          Minimum nuclear fraction threshold (default: 0.4)
      --metadata             Optional metadata CSV file
      --help                 Show this help message
      
    Examples:
      nextflow run main.nf --mapping_dirs mapping_dirs.csv --outputDir results
      nextflow run main.nf --mapping_dirs mapping_dirs.csv --outputDir results --cellbender
      nextflow run main.nf --mapping_dirs mapping_dirs.csv --outputDir results --cellbender --gpu
      nextflow run main.nf --mapping_dirs mapping_dirs.csv --outputDir results --multiome
      nextflow run main.nf --mapping_dirs mapping_dirs.csv --outputDir results --report --book
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

// =============================================================================
// MAIN WORKFLOW
// =============================================================================
workflow {
    // Log pipeline parameters
    log.info """
    ===================================
    scQC-flow pipeline
    ===================================
    Mapping directories : ${params.mapping_dirs}
    Output directory    : ${params.outputDir}
    Multiome mode       : ${params.multiome}
    CellBender          : ${params.cellbender}
    GPU acceleration    : ${params.gpu}
    Generate reports    : ${params.report}
    Generate book       : ${params.book}
    Metadata file       : ${params.metadata}
    QC Thresholds:
      Max mitochondrial % : ${params.max_mito}
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
    report_template_path = file("${projectDir}/modules/reports/seurat_template.qmd")
    combined_template_path = file("${projectDir}/modules/reports/combined_template.qmd")
    book_template_path = file("${projectDir}/modules/reports/book_template/")

    // Standard module scripts
    scdbl_script_path = file("${projectDir}/modules/scdbl/run_scdbl.R")
    seurat_script_path = file("${projectDir}/modules/seurat/make_seurat.R")

    // Multiome-specific scripts (extract Gene Expression modality from H5)
    scdbl_multiome_script_path = file("${projectDir}/modules/multiome/run_scdbl_multiome.R")
    seurat_multiome_script_path = file("${projectDir}/modules/multiome/make_seurat_multiome.R")
    extract_gex_script_path = file("${projectDir}/modules/multiome/extract_gex_h5.R")

    // =========================================================================
    // RUN APPROPRIATE WORKFLOW
    // =========================================================================
    if (params.multiome) {
        log.info "Running MULTIOME workflow"
        MULTIOME_WORKFLOW(
            sampleChannelBase,
            dropletqc_script_path,
            scdbl_multiome_script_path,
            seurat_multiome_script_path,
            extract_gex_script_path,
            params.cellbender,
            params.gpu,
            params.max_mito,
            params.min_nuclear,
            params.metadata
        )
        seurat_results = MULTIOME_WORKFLOW.out.seurat_results
    } else {
        log.info "Running STANDARD single-cell workflow"
        STANDARD_WORKFLOW(
            sampleChannelBase,
            dropletqc_script_path,
            scdbl_script_path,
            seurat_script_path,
            params.cellbender,
            params.gpu,
            params.max_mito,
            params.min_nuclear,
            params.metadata
        )
        seurat_results = STANDARD_WORKFLOW.out.seurat_results
    }

    // =========================================================================
    // REPORTING
    // =========================================================================
    REPORTING(
        sampleChannelBase,
        seurat_results,
        report_template_path,
        combined_template_path,
        book_template_path,
        params.max_mito,
        params.min_nuclear,
        params.report,
        params.book
    )
}
