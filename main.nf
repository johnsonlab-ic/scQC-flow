nextflow.enable.dsl=2

// Import workflows
include { STANDARD_WORKFLOW } from './workflows/workflows'
include { MULTIOME_WORKFLOW } from './workflows/workflows'
include { RNA_MAPPING_WORKFLOW } from './workflows/workflows'
include { ALEVINFRY_MAPPING_WORKFLOW } from './workflows/workflows'
include { ALEVINFRY_QC_WORKFLOW } from './workflows/workflows'

// Import processes used standalone (barcode estimation for qc-only alevinfry mode)
include { BARCODE_ESTIMATION; BARCODE_REPORT } from './modules/barcode_estimation/barcode_estimation'

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
            --mapper           Mapper to use (cellranger or alevinfry)
            --cellrangerPath   Path to Cell Ranger binary or install directory
            --transcriptome    Path to Cell Ranger transcriptome reference

        alevin-fry params (when --mapper alevinfry):
            --alevinfry_index       Pre-built simpleaf index directory
            --alevinfry_t2g         Transcript-to-gene map (t2g_3col.tsv)
            --alevinfry_whitelist   Optional barcode permit list override
            --alevinfry_chemistry   Required explicit chemistry (3v2, 5v1, 5v2, 3v3, 5v3, 3v4, 3LT, multiome)
            --alevinfry_orientation Expected read orientation (auto, fw, rc, both; default: auto)
            --alevinfry_fasta       Reference FASTA (for index building)
            --alevinfry_gtf         Reference GTF (for index building)

        QC and downstream params:
            --ambient_method   Ambient RNA removal: 'decontx' (CPU, default) or 'cellbender'
            --gpu              Use GPU for CellBender (only if --ambient_method cellbender)
            --max_mito         Maximum mitochondrial percentage threshold
            --min_nuclear      Minimum nuclear fraction threshold (default: 0.4)
            --max_nuclear      Maximum nuclear fraction threshold (default: 1.0)
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
    Ambient method      : ${params.ambient_method}
    GPU acceleration    : ${params.gpu}
    Generate reports    : ${params.report}
    Generate book       : ${params.book}
    Metadata file       : ${params.metadata}
    QC Thresholds:
      Max mitochondrial % : ${params.max_mito}
      Min nuclear fraction: ${params.min_nuclear}
      Max nuclear fraction: ${params.max_nuclear}
    ===================================
    """

    def validRunModes = ['mapping', 'qc', 'both']
    if (!validRunModes.contains(params.run_mode)) {
        error "Invalid --run_mode '${params.run_mode}'. Must be one of: mapping, qc, both"
    }

    def validMappers = ['cellranger', 'alevinfry']
    if (!validMappers.contains(params.mapper)) {
        error "Invalid --mapper '${params.mapper}'. Must be one of: ${validMappers.join(', ')}"
    }

    if (params.run_mode in ['mapping', 'both']) {
        if (!params.raw_data_dir) {
            error "--raw_data_dir is required for run_mode: mapping or both"
        }
        if (!file(params.raw_data_dir).exists()) {
            error "--raw_data_dir does not exist: ${params.raw_data_dir}"
        }

        if (params.mapper == 'cellranger') {
            if (!params.transcriptome) {
                error "--transcriptome is required for cellranger mapping"
            }
            if (!params.cellrangerPath) {
                error "--cellrangerPath is required for cellranger mapping"
            }
        } else if (params.mapper == 'alevinfry') {
            if (!params.alevinfry_index && (!params.alevinfry_fasta || !params.alevinfry_gtf)) {
                error "Either --alevinfry_index or both --alevinfry_fasta and --alevinfry_gtf are required for alevinfry mapping"
            }
            if (params.alevinfry_index && !params.alevinfry_t2g) {
                error "--alevinfry_t2g is required when using a pre-built --alevinfry_index"
            }
            if (!params.alevinfry_chemistry) {
                error "--alevinfry_chemistry is required for alevinfry mapping"
            }
            def supportedChemistries = ['3v2', '5v1', '5v2', '3v3', '5v3', '3v4', '3LT', 'multiome']
            if (!supportedChemistries.contains(params.alevinfry_chemistry)) {
                if (params.alevinfry_chemistry == '10xv3') {
                    error "--alevinfry_chemistry=10xv3 is ambiguous for whitelist selection. Use explicit chemistry: 3v3 or 5v3"
                }
                error "Unsupported --alevinfry_chemistry '${params.alevinfry_chemistry}'. Supported: ${supportedChemistries.join(', ')}"
            }
            def supportedOrientations = ['auto', 'fw', 'rc', 'both']
            if (!supportedOrientations.contains(params.alevinfry_orientation)) {
                error "Unsupported --alevinfry_orientation '${params.alevinfry_orientation}'. Supported: ${supportedOrientations.join(', ')}"
            }
            if (!params.alevinfry_whitelist && !params.cellrangerPath) {
                error "Either --alevinfry_whitelist must be provided, or provide --cellrangerPath so whitelist can be derived from Cell Ranger barcodes"
            }
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
        if (params.mapper == 'alevinfry') {
            // For alevinfry QC, mapped_data_dir contains per-sample af_quant dirs
            channel
                .fromPath("${params.mapped_data_dir}/*", type: 'dir')
                .map { dir -> tuple(dir.name, dir) }
                .set { afQuantDirs }
        } else {
            channel
                .fromPath("${params.mapped_data_dir}/*", type: 'dir')
                .map { dir -> tuple(dir.name, dir) }
                .set { sampleChannelBase }
        }

    } else {
        // Each raw_data_dir subdirectory = one sample.
        channel
            .fromPath("${params.raw_data_dir}/*", type: 'dir')
            .map { dir -> tuple(dir.name, dir.name, dir.toString()) }
            .set { sampleChannelFastq }

        if (params.mapper == 'cellranger') {
            RNA_MAPPING_WORKFLOW(
                sampleChannelFastq,
                params.mapper,
                params.transcriptome,
                params.cellrangerPath
            )
            sampleChannelBase = RNA_MAPPING_WORKFLOW.out.sample_channel

        } else if (params.mapper == 'alevinfry') {
            ALEVINFRY_MAPPING_WORKFLOW(
                sampleChannelFastq,
                params.cellrangerPath,
                params.alevinfry_index,
                params.alevinfry_t2g,
                params.alevinfry_whitelist,
                params.alevinfry_chemistry,
                params.alevinfry_orientation,
                params.alevinfry_fasta,
                params.alevinfry_gtf
            )
        }
    }

    // Prepare script/template files as value channels
    dropletqc_script_path = file("${projectDir}/modules/dropletqc/run_dropletqc.R")

    // Standard module scripts
    scdbl_script_path = file("${projectDir}/modules/scdbl/run_scdbl.R")
    seurat_script_path = file("${projectDir}/modules/seurat/make_seurat.R")

    if (params.run_mode == 'mapping') {
        log.info "Run mode is 'mapping' only - QC will be skipped."

    } else if (params.mapper == 'alevinfry') {
        // =========================================================================
        // ALEVIN-FRY QC WORKFLOW
        // Barcode estimation ran at the end of ALEVINFRY_MAPPING_WORKFLOW ('both'
        // mode) or is run here from the provided quant dirs ('qc' mode).
        // Chain: ambient correction (decontx|cellbender) → doublets → Seurat
        // =========================================================================
        log.info "Running ALEVIN-FRY QC workflow (ambient_method: ${params.ambient_method})"

        if (params.run_mode == 'qc') {
            // QC-only: user provides quant dirs via --mapped_data_dir.
            // Run barcode estimation now before handing off to the QC workflow.
            barcode_script     = channel.value(file("${projectDir}/modules/barcode_estimation/run_barcode_estimation.R"))
            barcode_report_qmd = channel.value(file("${projectDir}/modules/barcode_estimation/barcode_estimation.qmd"))
            BARCODE_ESTIMATION(afQuantDirs, barcode_script)
            BARCODE_REPORT(
                BARCODE_ESTIMATION.out.knee_data.map { _id, csv -> csv }.collect(),
                BARCODE_ESTIMATION.out.nuclear_fraction.map { _id, csv -> csv }.collect(),
                barcode_report_qmd
            )
            qc_h5_ch             = BARCODE_ESTIMATION.out.h5_files
            qc_knee_data_ch      = BARCODE_ESTIMATION.out.knee_data
            qc_knee_params_ch    = BARCODE_ESTIMATION.out.knee_params
            qc_nuclear_fraction_ch = BARCODE_ESTIMATION.out.nuclear_fraction
            qc_quant_dirs_ch     = afQuantDirs
        } else {
            // 'both' mode: barcode estimation already ran inside ALEVINFRY_MAPPING_WORKFLOW.
            qc_h5_ch             = ALEVINFRY_MAPPING_WORKFLOW.out.h5_files
            qc_knee_data_ch      = ALEVINFRY_MAPPING_WORKFLOW.out.knee_data
            qc_knee_params_ch    = ALEVINFRY_MAPPING_WORKFLOW.out.knee_params
            qc_nuclear_fraction_ch = ALEVINFRY_MAPPING_WORKFLOW.out.nuclear_fraction
            qc_quant_dirs_ch     = ALEVINFRY_MAPPING_WORKFLOW.out.quant_dirs
        }

        ALEVINFRY_QC_WORKFLOW(
            qc_h5_ch,
            qc_knee_data_ch,
            qc_knee_params_ch,
            qc_nuclear_fraction_ch,
            qc_quant_dirs_ch,
            scdbl_script_path,
            seurat_script_path,
            params.ambient_method,
            params.gpu,
            params.max_mito,
            params.min_nuclear,
            params.max_nuclear,
            params.metadata
        )

    } else {
        // =========================================================================
        // STANDARD (CELL RANGER) QC WORKFLOW
        // =========================================================================
        log.info "Running STANDARD Cell Ranger QC workflow"
        STANDARD_WORKFLOW(
            sampleChannelBase,
            dropletqc_script_path,
            scdbl_script_path,
            seurat_script_path,
            params.ambient_method == 'cellbender',
            params.gpu,
            params.max_mito,
            params.min_nuclear,
            params.metadata
        )
        seurat_results = STANDARD_WORKFLOW.out.seurat_results
    }

}
