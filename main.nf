nextflow.enable.dsl=2

include { MAPPING           } from './workflows/workflows'
include { SAMPLE_METADATA   } from './workflows/workflows'
include { AMBIENT           } from './workflows/workflows'
include { AMBIENT_DE_WF     } from './workflows/workflows'
include { QC                } from './workflows/workflows'
include { HVG               } from './workflows/workflows'
include { INTEGRATION       } from './workflows/workflows'
include { ANNOTATION        } from './workflows/workflows'
include { BUILD_SITE        } from './modules/reports/reports'

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
        --siteDir           Optional secondary publish directory for the assembled
                    report site bundle (default: testing/ under the repo)
        --build_site        Build the final assembled report site bundle
                    (default: true)
        --run_ambient       Run ambient RNA removal with decontX (default: true)
        --run_qc            Run cell-level QC (default: true; requires --run_ambient)
        --run_hvg           Run HVG selection (default: true; requires --run_qc)
        --run_integration   Run Harmony integration (default: true; requires --run_hvg)
        --run_annotation    Run pseudobulk marker-gene annotation (default: false;
                    requires --run_integration)
        --metadata_csv      Path to metadata CSV (required when --run_integration)
        --metadata_id_col   Column in metadata CSV that maps to sample IDs (default: sample_id)
        --metadata_vars     Space-separated metadata columns for Harmony correction
                            (e.g. "brainregion condition")
        --hvg_n_hvgs        Number of HVGs to select (default: 2000)
        --annotation_marker_csv   Long-format marker CSV with either label/symbol
                    or cell_type/marker_gene columns
        --annotation_sel_res      Single clustering resolution to carry from
                    integration into annotation (default: 0.2)
        --barcode_v2_splice_context  Splice QC context for Barcode_estimation_v2
                         (snrna, scrna, auto; default: snrna)
        --barcode_v2_ed_fdr          EmptyDrops FDR cutoff for Barcode_estimation_v2
                         (default: 0.001)
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
    if (params.metadata_csv && !file(params.metadata_csv).exists()) {
        error "--metadata_csv does not exist: ${params.metadata_csv}"
    }
    if (params.annotation_marker_csv && !file(params.annotation_marker_csv).exists()) {
        error "--annotation_marker_csv does not exist: ${params.annotation_marker_csv}"
    }
    if (params.run_annotation && !params.annotation_marker_csv) {
        error "--annotation_marker_csv is required when --run_annotation is true"
    }
    if (params.run_annotation && !params.run_integration) {
        error "--run_annotation requires --run_integration"
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
    build_site     : ${params.build_site}
    run_ambient    : ${params.run_ambient}
    run_qc         : ${params.run_qc}
    run_hvg        : ${params.run_hvg}
    run_integration: ${params.run_integration}
    run_annotation : ${params.run_annotation}
    metadata_csv   : ${params.metadata_csv ?: 'not provided'}
    annotation_csv : ${params.annotation_marker_csv ?: 'not provided'}
    outputDir      : ${params.outputDir}
    siteDir        : ${params.siteDir ?: 'disabled'}
    ===================================
    """

    // ------------------------------------------------------------------
    // Sample discovery: each subdir of raw_data_dir is one sample
    // ------------------------------------------------------------------
    samples_ch = channel
        .fromPath("${params.raw_data_dir}/*", type: 'dir')
        .map { dir -> tuple(dir.name, dir.name, dir.toString()) }

    def sample_metadata_ch
    if (params.metadata_csv) {
        SAMPLE_METADATA(samples_ch.map { sampleId, _sampleName, _fastqPath -> sampleId }.collect())
        sample_metadata_ch = SAMPLE_METADATA.out.sample_metadata
    } else {
        sample_metadata_ch = channel.value(file("${projectDir}/templates/NO_FILE"))
    }

    // ------------------------------------------------------------------
    // 1. Mapping (always runs)
    // ------------------------------------------------------------------
    MAPPING(samples_ch)

    // Collect per-report site bundles for BUILD_SITE
    report_sites = MAPPING.out.report

    // ------------------------------------------------------------------
    // 2. Ambient cleanup (opt-in via --run_ambient)
    // ------------------------------------------------------------------
    if (params.run_ambient) {
        AMBIENT(
            MAPPING.out.h5_files,
            MAPPING.out.knee_data
        )
        report_sites = report_sites.mix(AMBIENT.out.report)

        // Ambient DE analysis (always runs if ambient runs)
        AMBIENT_DE_WF(
            MAPPING.out.h5_files,
            MAPPING.out.knee_data,
            AMBIENT.out.h5_files
        )

        // --------------------------------------------------------------
        // 3. Cell-level QC (opt-in via --run_qc, requires ambient)
        // --------------------------------------------------------------
        if (params.run_qc) {
            QC(AMBIENT.out.h5_files, sample_metadata_ch)
            report_sites = report_sites.mix(QC.out.report)

            // ----------------------------------------------------------
            // 4. HVG selection (opt-in via --run_hvg, requires qc)
            // ----------------------------------------------------------
            if (params.run_hvg) {
                HVG(
                    AMBIENT.out.h5_files,
                    QC.out.qc_metrics,
                    AMBIENT_DE_WF.out.de_table,
                    AMBIENT_DE_WF.out.pb_empties
                )
                report_sites = report_sites.mix(HVG.out.report)

                // ------------------------------------------------------
                // 5. Integration (opt-in via --run_integration, requires hvg)
                // ------------------------------------------------------
                if (params.run_integration) {
                    if (!params.metadata_csv) {
                        error "--metadata_csv is required when --run_integration is true"
                    }
                    INTEGRATION(
                        HVG.out.hvg_counts,
                        HVG.out.dbl_hvg_counts,
                        QC.out.qc_metrics
                    )
                    report_sites = report_sites.mix(INTEGRATION.out.report)

                    if (params.run_annotation) {
                        ANNOTATION(
                            AMBIENT.out.h5_files,
                            INTEGRATION.out.integration_dt
                        )
                        report_sites = report_sites.mix(ANNOTATION.out.report)
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // 6. Site bundle (collects all produced report bundles + builds index)
    // ------------------------------------------------------------------
    if (params.build_site) {
        BUILD_SITE(report_sites.collect())
    }
}
