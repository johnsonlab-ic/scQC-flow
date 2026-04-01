nextflow.enable.dsl=2

// Import workflows
include { STANDARD_WORKFLOW } from './workflows/workflows'
include { MULTIOME_WORKFLOW } from './workflows/workflows'
include { REPORTING } from './workflows/workflows'
include { ATAC_REPORTING } from './workflows/workflows'
include { RNA_MAPPING_WORKFLOW } from './workflows/workflows'

// =============================================================================
// PARAMETERS
// Runtime params are defined in nextflow.config and can be overridden by CLI.
// =============================================================================

// =============================================================================
// HELP MESSAGE
// =============================================================================
def helpMessage() {
    log.info """
    ===================================
    scQC-flow
    ===================================
    
    Usage:
        nextflow run main.nf -profile offline [params]

        Main params (set in nextflow.config, override on CLI as needed):
            --run_mode         mapping | qc | both
            --raw_data_dir     Directory of FASTQ sample folders (mapping or both mode)
            --mapped_data_dir  Directory of mapped sample folders (qc mode)
            --mapper           Mapper to use (currently: cellranger)
            --cellrangerPath   Path to Cell Ranger binary or install directory
            --transcriptome    Path to Cell Ranger transcriptome reference
            --cellbender       Enable CellBender
            --gpu              Use GPU for CellBender
            --max_mito         Maximum mitochondrial percentage threshold
            --min_nuclear      Minimum nuclear fraction threshold
            --metadata         Optional metadata CSV
            --outputDir        Output directory
            --report           Generate per-sample reports
            --book             Generate combined report book
      --help                Show this help message

        Examples:
            nextflow run main.nf -profile offline --run_mode mapping --raw_data_dir /path/to/raw --cellrangerPath /path/to/cellranger --transcriptome /path/to/ref
            nextflow run main.nf -profile offline --run_mode qc --mapped_data_dir /path/to/mapping
            nextflow run main.nf -profile offline --run_mode both --raw_data_dir /path/to/raw --cellrangerPath /path/to/cellranger --transcriptome /path/to/ref
    """
}

// =============================================================================
// MAIN WORKFLOW
// =============================================================================
workflow {
    if (params.help) {
        helpMessage()
        return
    }

    // Log pipeline parameters
    log.info """
    ===================================
    scQC-flow pipeline
    ===================================
    Output directory    : ${params.outputDir}
    Run mode            : ${params.run_mode}
    Mapper              : ${params.mapper}
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

    def validRunModes = ['mapping', 'qc', 'both']
    if (!validRunModes.contains(params.run_mode)) {
        error "Invalid --run_mode '${params.run_mode}'. Must be one of: mapping, qc, both"
    }

    if (params.run_mode in ['mapping', 'both']) {
        if (!params.raw_data_dir) {
            error "--raw_data_dir is required for run_mode: mapping or both"
        }
        if (!file(params.raw_data_dir).exists()) {
            error "--raw_data_dir does not exist: ${params.raw_data_dir}"
        }
        if (!params.transcriptome) {
            error "--transcriptome is required for run_mode: mapping or both"
        }
        if (!params.cellrangerPath) {
            error "--cellrangerPath is required for run_mode: mapping or both"
        }
    }

    if (params.run_mode == 'qc') {
        if (!params.mapped_data_dir) {
            error "--mapped_data_dir is required for run_mode: qc"
        }
        if (!file(params.mapped_data_dir).exists()) {
            error "--mapped_data_dir does not exist: ${params.mapped_data_dir}"
        }
    }

    // Build sample channels directly from Nextflow params.
    if (params.run_mode == 'qc') {
        channel
            .fromPath("${params.mapped_data_dir}/*", type: 'dir')
            .map { dir -> tuple(dir.name, dir) }
            .set { sampleChannelBase }

    } else {
        // Each raw_data_dir subdirectory = one sample.
        channel
            .fromPath("${params.raw_data_dir}/*", type: 'dir')
            .map { dir -> tuple(dir.name, dir.name, dir.toString()) }
            .set { sampleChannelFastq }

        RNA_MAPPING_WORKFLOW(
            sampleChannelFastq,
            params.mapper,
            params.transcriptome,
            params.cellrangerPath
        )

        sampleChannelBase = RNA_MAPPING_WORKFLOW.out.sample_channel
    }

    // Prepare script/template files as value channels
    dropletqc_script_path = file("${projectDir}/modules/dropletqc/run_dropletqc.R")
    report_template_path = file("${projectDir}/modules/reports/seurat_template.qmd")
    combined_template_path = file("${projectDir}/modules/reports/combined_template.qmd")
    book_template_path = file("${projectDir}/modules/reports/book_template/")

    // Standard module scripts
    scdbl_script_path = file("${projectDir}/modules/scdbl/run_scdbl.R")
    seurat_script_path = file("${projectDir}/modules/seurat/make_seurat.R")

    if (params.run_mode == 'mapping') {
        log.info "Run mode is 'mapping'. Mapping stage completed; skipping QC and reporting."
    } else {
        // =========================================================================
        // STANDARD RNA WORKFLOW
        // =========================================================================
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
        cellbender_comparison_results = STANDARD_WORKFLOW.out.cellbender_comparison_results

        // =========================================================================
        // REPORTING
        // =========================================================================
        REPORTING(
            sampleChannelBase,
            seurat_results,
            cellbender_comparison_results,
            report_template_path,
            combined_template_path,
            book_template_path,
            params.max_mito,
            params.min_nuclear,
            params.report,
            params.book
        )
    }

    // ATAC reporting is intentionally disabled while pipeline scope is RNA-only.
}
