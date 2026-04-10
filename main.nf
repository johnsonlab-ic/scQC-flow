nextflow.enable.dsl=2

include { MAPPING           } from './workflows/workflows'
include { AMBIENT           } from './workflows/workflows'
include { QC                } from './workflows/workflows'
include { HVG               } from './workflows/workflows'
include { INTEGRATION       } from './workflows/workflows'
include { AMBIENT_DE        } from './modules/ambient_de/ambient_de'
include { STAGE_RAW_H5     } from './modules/ambient_de/ambient_de'
include { INDEX_REPORT      } from './modules/reports/reports'

// =============================================================================
// HELP
// =============================================================================
def helpMessage() {
    log.info """
    ===================================
    scQC-flow
    ===================================

    Usage:
        nextflow run main.nf -profile <profile> [params]

    Required params:
        --raw_data_dir      Parent directory; each subdirectory is one sample
                            (each subdir must contain *_R1_* and *_R2_* FASTQs)
        --cellrangerPath    Path to Cell Ranger install directory
                            (whitelists are extracted from lib/python/cellranger/barcodes/)
        --genome_fasta      Reference genome FASTA (for building simpleaf index)
        --genome_gtf        Reference genome GTF

    Optional:
        --chemistry         10x chemistry string (default: 3v4)
                            Accepted: 3v2, 3v3, 3v4, 3LT, 5v1, 5v2, 5v3, multiome
        --outputDir         Output directory (default: results)
        --run_ambient       Run ambient RNA removal with decontX (default: true)
        --run_qc            Run cell-level QC (default: true; requires --run_ambient)
        --run_hvg           Run HVG selection (default: true; requires --run_qc)
        --run_integration   Run Harmony integration (default: true; requires --run_hvg)
        --metadata_csv      Path to metadata CSV (required when --run_integration)
        --metadata_id_col   Column in metadata CSV that maps to sample IDs (default: sample_id)
        --metadata_vars     Space-separated metadata columns for Harmony correction
                            (e.g. "brainregion condition")
        --hvg_n_hvgs        Number of HVGs to select (default: 4000)
        --help              Show this message

    Profiles:
        offline             Local execution
        imperial            Imperial HPC (PBS)
        slurm               Generic SLURM cluster
        dsi                 DSI cluster

    Example:
        nextflow run main.nf -profile offline \\
            --raw_data_dir    data/mini_landmark_data/raw \\
            --cellrangerPath  /path/to/cellranger-9.0.0 \\
            --genome_fasta    /path/to/refdata-gex-GRCh38/fasta/genome.fa \\
            --genome_gtf      /path/to/refdata-gex-GRCh38/genes/genes.gtf \\
            --chemistry       3v4 \\
            --outputDir       results/
    """
}

// =============================================================================
// MAIN
// =============================================================================
workflow {

    if (params.help) {
        helpMessage()
        return
    }

    // ------------------------------------------------------------------
    // Parameter validation
    // ------------------------------------------------------------------
    if (!params.raw_data_dir) {
        error "--raw_data_dir is required"
    }
    if (!file(params.raw_data_dir).exists()) {
        error "--raw_data_dir does not exist: ${params.raw_data_dir}"
    }
    if (!params.cellrangerPath) {
        error "--cellrangerPath is required (Cell Ranger install directory)"
    }
    if (!params.genome_fasta) {
        error "--genome_fasta is required (reference genome FASTA)"
    }
    if (!file(params.genome_fasta).exists()) {
        error "--genome_fasta does not exist: ${params.genome_fasta}"
    }
    if (!params.genome_gtf) {
        error "--genome_gtf is required (reference genome GTF)"
    }
    if (!file(params.genome_gtf).exists()) {
        error "--genome_gtf does not exist: ${params.genome_gtf}"
    }

    // ------------------------------------------------------------------
    // Log run info
    // ------------------------------------------------------------------
    log.info """
    ===================================
    scQC-flow
    ===================================
    raw_data_dir   : ${params.raw_data_dir}
    cellrangerPath : ${params.cellrangerPath}
    genome_fasta   : ${params.genome_fasta}
    genome_gtf     : ${params.genome_gtf}
    chemistry      : ${params.chemistry}
    run_ambient    : ${params.run_ambient}
    run_qc         : ${params.run_qc}
    run_hvg        : ${params.run_hvg}
    run_integration: ${params.run_integration}
    outputDir      : ${params.outputDir}
    ===================================
    """

    // ------------------------------------------------------------------
    // Sample discovery: each subdir of raw_data_dir is one sample
    // ------------------------------------------------------------------
    samples_ch = channel
        .fromPath("${params.raw_data_dir}/*", type: 'dir')
        .map { dir -> tuple(dir.name, dir.name, dir.toString()) }

    // ------------------------------------------------------------------
    // 1. Mapping (always runs)
    // ------------------------------------------------------------------
    MAPPING(samples_ch)

    // Collect report HTMLs for INDEX_REPORT
    report_htmls = MAPPING.out.report

    // ------------------------------------------------------------------
    // 2. Ambient cleanup (opt-in via --run_ambient)
    // ------------------------------------------------------------------
    if (params.run_ambient) {
        AMBIENT(
            MAPPING.out.h5_files,
            MAPPING.out.knee_data
        )
        report_htmls = report_htmls.mix(AMBIENT.out.report)

        // Ambient DE analysis (always runs if ambient runs)
        // STAGE_RAW_H5 symlinks af_counts_mat.h5 → barcode_matrix_<id>.h5
        // so multiple samples can be collected without filename collision.
        STAGE_RAW_H5(MAPPING.out.h5_files)

        AMBIENT_DE(
            STAGE_RAW_H5.out.h5.collect(),
            MAPPING.out.knee_data.map { _id, csv -> csv }.collect(),
            AMBIENT.out.h5_files.map { _id, h5 -> h5 }.collect(),
            channel.value(file(params.genome_gtf)),
            channel.value(file("${projectDir}/modules/ambient_de/ambient_de.R"))
        )

        // --------------------------------------------------------------
        // 3. Cell-level QC (opt-in via --run_qc, requires ambient)
        // --------------------------------------------------------------
        if (params.run_qc) {
            QC(AMBIENT.out.h5_files)
            report_htmls = report_htmls.mix(QC.out.report)

            // ----------------------------------------------------------
            // 4. HVG selection (opt-in via --run_hvg, requires qc)
            // ----------------------------------------------------------
            if (params.run_hvg) {
                HVG(
                    AMBIENT.out.h5_files,
                    QC.out.qc_metrics,
                    AMBIENT_DE.out.de_table,
                    AMBIENT_DE.out.pb_empties
                )
                report_htmls = report_htmls.mix(HVG.out.report)

                // ------------------------------------------------------
                // 5. Integration (opt-in via --run_integration, requires hvg)
                // ------------------------------------------------------
                if (params.run_integration) {
                    if (!params.metadata_csv) {
                        error "--metadata_csv is required when --run_integration is true"
                    }
                    INTEGRATION(
                        HVG.out.hvg_counts,
                        QC.out.qc_metrics
                    )
                    report_htmls = report_htmls.mix(INTEGRATION.out.report)
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // 6. Index page (standalone — collects all produced report HTMLs)
    // ------------------------------------------------------------------
    INDEX_REPORT(report_htmls.collect())
}
