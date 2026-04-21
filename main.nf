nextflow.enable.dsl=2

include { MAPPING           } from './workflows/workflows'
include { SAMPLE_METADATA   } from './workflows/workflows'
include { AMBIENT           } from './workflows/workflows'
include { AMBIENT_DE_WF     } from './workflows/workflows'
include { QC                } from './workflows/workflows'
include { HVG               } from './workflows/workflows'
include { INTEGRATION       } from './workflows/workflows'
include { ANNOTATION        } from './workflows/workflows'
include { ZOOMS             } from './workflows/workflows'
include { REPORT_SITE       } from './modules/reports/reports'
include { EXPORT_SCANPY     } from './modules/export/export_sc'
include { EXPORT_SEURAT     } from './modules/export/export_sc'

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
        --run_annotation    Run pseudobulk marker-gene annotation (default: false;
                    requires --run_integration)
        --export            Export integrated objects (default: none).
                    Accepted: none, anndata, seurat, both
        zoom                Nested zoom configuration map in nextflow.config.
                Each zoom can subset by integration cluster or annotation label,
                then rerun HVG selection, integration, and marker analysis.
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
    if (params.run_qc && !params.run_ambient) {
        error "--run_qc requires --run_ambient"
    }
    if (params.run_hvg && !params.run_qc) {
        error "--run_hvg requires --run_qc"
    }
    if (params.run_integration && !params.run_hvg) {
        error "--run_integration requires --run_hvg"
    }
    if (params.run_annotation && !params.run_integration) {
        error "--run_annotation requires --run_integration"
    }

    def exportMode = (params.export ?: 'none').toString().trim().toLowerCase()
    def validExportModes = ['none', 'anndata', 'seurat', 'both']
    if (!(exportMode in validExportModes)) {
        error "--export must be one of: ${validExportModes.join(', ')}"
    }
    if (exportMode != 'none' && !params.run_integration) {
        error "--export requires --run_integration"
    }

    def normalizedZoomSpecs = []
    def zoomConfig = (params.zoom instanceof Map) ? params.zoom : [:]
    if (zoomConfig.enabled) {
        if (!params.run_integration) {
            error "zoom.enabled requires --run_integration"
        }

        def rawItems = zoomConfig.items ?: []
        if (!(rawItems instanceof List)) {
            error "params.zoom.items must be a list of zoom definitions"
        }

        rawItems.each { rawSpec ->
            if (!(rawSpec instanceof Map)) {
                error "Each zoom definition must be a map"
            }

            def zoomName = rawSpec.name?.toString()?.trim()
            if (!zoomName) {
                error "Each zoom definition requires a non-empty 'name'"
            }

            def zoomSource = rawSpec.source?.toString()?.trim()
            if (!(zoomSource in ['cluster', 'annotation_label'])) {
                error "Zoom '${zoomName}' has invalid source '${zoomSource}'. Use 'cluster' or 'annotation_label'."
            }

            if (zoomSource == 'annotation_label' && !params.run_annotation) {
                error "Zoom '${zoomName}' uses source='annotation_label' and therefore requires --run_annotation"
            }

            def zoomValues = []
            def zoomLabel = null
            if (zoomSource == 'annotation_label') {
                zoomLabel = rawSpec.label?.toString()?.trim()
                if (!zoomLabel && rawSpec.values != null) {
                    def legacyValuesRaw = rawSpec.values
                    def legacyValues = legacyValuesRaw instanceof Collection ? legacyValuesRaw.collect { value -> value.toString() } : [legacyValuesRaw.toString()]
                    legacyValues = legacyValues.findAll { value -> value?.trim() }
                    if (legacyValues.size() > 1) {
                        error "Zoom '${zoomName}' uses source='annotation_label' and must define only one label"
                    }
                    zoomLabel = legacyValues ? legacyValues[0] : null
                }
                if (!zoomLabel) {
                    zoomLabel = zoomName
                }
            } else {
                def zoomValuesRaw = rawSpec.values
                zoomValues = zoomValuesRaw instanceof Collection ? zoomValuesRaw.collect { value -> value.toString() } : (zoomValuesRaw != null ? [zoomValuesRaw.toString()] : [])
                zoomValues = zoomValues.findAll { value -> value?.trim() }
                if (!zoomValues) {
                    error "Zoom '${zoomName}' must define one or more 'values'"
                }
            }

            def metadataVarsRaw = rawSpec.containsKey('metadata_vars') ? rawSpec.metadata_vars : (zoomConfig.default_metadata_vars ?: params.metadata_vars)
            def metadataVars = metadataVarsRaw instanceof Collection ? metadataVarsRaw.collect { value -> value.toString() }.join(' ') : (metadataVarsRaw ?: '').toString()

            def integrationLeidenRaw = rawSpec.integration_leiden_res ?: zoomConfig.default_integration_leiden_res ?: params.integration_leiden_res
            def integrationLeiden = integrationLeidenRaw instanceof Collection ? integrationLeidenRaw.collect { value -> value.toString() }.join(' ') : integrationLeidenRaw.toString()

            def normalizedSpec = [
                name: zoomName,
                source: zoomSource,
                sample_min_cells: (rawSpec.sample_min_cells ?: zoomConfig.default_sample_min_cells ?: 10) as Integer,
                hvg_n_hvgs: (rawSpec.hvg_n_hvgs ?: zoomConfig.default_hvg_n_hvgs ?: params.hvg_n_hvgs) as Integer,
                metadata_vars: metadataVars,
                integration_n_dims: (rawSpec.integration_n_dims ?: zoomConfig.default_integration_n_dims ?: params.integration_n_dims) as Integer,
                integration_theta: (rawSpec.integration_theta ?: zoomConfig.default_integration_theta ?: params.integration_theta) as BigDecimal,
                integration_leiden_res: integrationLeiden,
                integration_dbl_res: (rawSpec.integration_dbl_res ?: zoomConfig.default_integration_dbl_res ?: params.integration_dbl_res) as BigDecimal,
                integration_dbl_cl_prop: (rawSpec.integration_dbl_cl_prop ?: zoomConfig.default_integration_dbl_cl_prop ?: params.integration_dbl_cl_prop) as BigDecimal,
                exclude_mito: rawSpec.containsKey('exclude_mito') ? rawSpec.exclude_mito : (zoomConfig.containsKey('default_exclude_mito') ? zoomConfig.default_exclude_mito : params.exclude_mito),
                marker_sel_res: (rawSpec.marker_sel_res ?: zoomConfig.default_marker_sel_res ?: '0.2').toString(),
                marker_min_cl_size: (rawSpec.marker_min_cl_size ?: zoomConfig.default_marker_min_cl_size ?: 100) as Integer,
                marker_min_cells: (rawSpec.marker_min_cells ?: zoomConfig.default_marker_min_cells ?: 10) as Integer,
                marker_top_n: (rawSpec.marker_top_n ?: zoomConfig.default_marker_top_n ?: 10) as Integer,
                marker_min_cpm: (rawSpec.marker_min_cpm ?: zoomConfig.default_marker_min_cpm ?: params.annotation_min_cpm_mkr) as Integer,
                marker_fdr_cut: (rawSpec.marker_fdr_cut ?: zoomConfig.default_marker_fdr_cut ?: params.annotation_fdr_cut) as BigDecimal,
                marker_max_zero_p: (rawSpec.marker_max_zero_p ?: zoomConfig.default_marker_max_zero_p ?: params.annotation_max_zero_p) as BigDecimal,
            ]

            if (zoomSource == 'cluster') {
                normalizedSpec.values = zoomValues
                normalizedSpec.cluster_res = (rawSpec.cluster_res ?: zoomConfig.default_cluster_res ?: '0.2').toString()
            } else {
                normalizedSpec.label = zoomLabel
            }

            normalizedZoomSpecs << normalizedSpec
        }

        def zoomNames = normalizedZoomSpecs.collect { spec -> spec.name }
        if (zoomNames.toSet().size() != zoomNames.size()) {
            error "Zoom names must be unique: ${zoomNames}"
        }
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
    run_annotation : ${params.run_annotation}
    export         : ${exportMode}
    zooms          : ${normalizedZoomSpecs ? normalizedZoomSpecs.collect { spec -> spec.name }.join(', ') : 'disabled'}
    metadata_csv   : ${params.metadata_csv ?: 'not provided'}
    annotation_csv : ${params.annotation_marker_csv ?: 'not provided'}
    outputDir      : ${params.outputDir}
    ===================================
    """

    def landingPayload = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson([
            pipeline: 'scQC-flow',
            run_name: workflow.runName,
            profile: workflow.profile ?: 'default',
            generated_at: new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()),
            output_dir: params.outputDir,
            reports_dir: "${params.outputDir}/reports",
            pipeline_steps: [
                mapping: true,
                barcode_v2: true,
                ambient: params.run_ambient,
                qc: params.run_ambient && params.run_qc,
                hvg: params.run_ambient && params.run_qc && params.run_hvg,
                integration: params.run_ambient && params.run_qc && params.run_hvg && params.run_integration,
                annotation: params.run_ambient && params.run_qc && params.run_hvg && params.run_integration && params.run_annotation,
                export: params.run_ambient && params.run_qc && params.run_hvg && params.run_integration && exportMode != 'none',
                zooms: normalizedZoomSpecs.size() > 0,
            ],
            params: [
                input_output: [
                    raw_data_dir: params.raw_data_dir,
                    outputDir: params.outputDir,
                    cellrangerPath: params.cellrangerPath,
                ],
                references: [
                    genome_fasta: params.genome_fasta,
                    genome_gtf: params.genome_gtf,
                    chemistry: params.chemistry,
                ],
                workflow_flags: [
                    run_ambient: params.run_ambient,
                    run_qc: params.run_qc,
                    run_hvg: params.run_hvg,
                    run_integration: params.run_integration,
                    run_annotation: params.run_annotation,
                    export: exportMode,
                ],
                barcode_v2: [
                    barcode_v2_min_umis_empty: params.barcode_v2_min_umis_empty,
                    barcode_v2_ed_niters: params.barcode_v2_ed_niters,
                    barcode_v2_ed_fdr: params.barcode_v2_ed_fdr,
                    barcode_v2_splice_context: params.barcode_v2_splice_context,
                    barcode_v2_low_count_strategy: params.barcode_v2_low_count_strategy,
                ],
                metadata: [
                    metadata_csv: params.metadata_csv,
                    metadata_id_col: params.metadata_id_col,
                    metadata_vars: params.metadata_vars,
                ],
                hvg: [
                    hvg_n_hvgs: params.hvg_n_hvgs,
                ],
                integration: [
                    integration_n_dims: params.integration_n_dims,
                    integration_theta: params.integration_theta,
                    integration_leiden_res: params.integration_leiden_res,
                    integration_dbl_res: params.integration_dbl_res,
                    integration_dbl_cl_prop: params.integration_dbl_cl_prop,
                    exclude_mito: params.exclude_mito,
                ],
                annotation: [
                    annotation_marker_csv: params.annotation_marker_csv,
                    annotation_sel_res: params.annotation_sel_res,
                    annotation_min_cl_size: params.annotation_min_cl_size,
                    annotation_min_cells: params.annotation_min_cells,
                    annotation_not_ok_re: params.annotation_not_ok_re,
                    annotation_min_cpm_mkr: params.annotation_min_cpm_mkr,
                    annotation_top_n: params.annotation_top_n,
                    annotation_fdr_cut: params.annotation_fdr_cut,
                    annotation_max_zero_p: params.annotation_max_zero_p,
                ],
                export: [
                    mode: exportMode,
                ],
                qc_thresholds: [
                    qc_hard_min_counts: params.qc_hard_min_counts,
                    qc_hard_min_feats: params.qc_hard_min_feats,
                    qc_hard_max_mito: params.qc_hard_max_mito,
                    qc_min_counts: params.qc_min_counts,
                    qc_min_feats: params.qc_min_feats,
                    qc_max_mito: params.qc_max_mito,
                    qc_min_mito: params.qc_min_mito,
                    qc_max_splice: params.qc_max_splice,
                    qc_min_splice: params.qc_min_splice,
                ],
                zoom: [
                    enabled: zoomConfig.enabled ?: false,
                    items: normalizedZoomSpecs,
                ],
            ],
        ])
    )

    // ------------------------------------------------------------------
    // Sample discovery: each subdir of raw_data_dir is one sample
    // ------------------------------------------------------------------
    samples_ch = channel
        .fromPath("${params.raw_data_dir}/*", type: 'dir')
        .map { dir -> tuple(dir.name, dir.name, dir.toString()) }

    def sample_metadata_ch
    def annotation_cell_labels_ch = channel.value(file("${projectDir}/templates/NO_FILE"))
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
    def report_pages = MAPPING.out.report

    // ------------------------------------------------------------------
    // 2. Ambient cleanup (opt-in via --run_ambient)
    // ------------------------------------------------------------------
    if (params.run_ambient) {
        AMBIENT(
            MAPPING.out.h5_files,
            MAPPING.out.knee_data
        )
        report_pages = report_pages.mix(AMBIENT.out.report)

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
            report_pages = report_pages.mix(QC.out.report)

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
                report_pages = report_pages.mix(HVG.out.report)

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
                    report_pages = report_pages.mix(INTEGRATION.out.report)

                    if (params.run_annotation) {
                        ANNOTATION(
                            AMBIENT.out.h5_files,
                            INTEGRATION.out.integration_dt
                        )
                        annotation_cell_labels_ch = ANNOTATION.out.cell_labels
                        report_pages = report_pages.mix(ANNOTATION.out.report)
                    }

                    if (normalizedZoomSpecs) {
                        ZOOMS(
                            channel.fromList(normalizedZoomSpecs),
                            AMBIENT.out.h5_files,
                            QC.out.qc_metrics,
                            AMBIENT_DE_WF.out.de_table,
                            INTEGRATION.out.integration_dt,
                            annotation_cell_labels_ch
                        )
                        report_pages = report_pages.mix(ZOOMS.out.report)
                    }

                    if (exportMode in ['anndata', 'both']) {
                        EXPORT_SCANPY(
                            AMBIENT.out.h5_files.map { _id, h5 -> h5 }.collect(),
                            QC.out.qc_metrics.map { _id, csv -> csv }.collect(),
                            INTEGRATION.out.integration_dt,
                            annotation_cell_labels_ch,
                            channel.value(file(params.genome_gtf)),
                            channel.value(file("${projectDir}/modules/export/export_anndata.py"))
                        )
                    }

                    if (exportMode in ['seurat', 'both']) {
                        EXPORT_SEURAT(
                            AMBIENT.out.h5_files.map { _id, h5 -> h5 }.collect(),
                            QC.out.qc_metrics.map { _id, csv -> csv }.collect(),
                            INTEGRATION.out.integration_dt,
                            annotation_cell_labels_ch,
                            channel.value(file(params.genome_gtf)),
                            channel.value(file("${projectDir}/modules/export/export_seurat.R")),
                            channel.value(file("${projectDir}/modules/export/export_utils.R"))
                        )
                    }
                }
            }
        }
    }

    REPORT_SITE(
        report_pages.collect(),
        channel.value(file("${projectDir}/modules/reports/landing_page.qmd")),
        landingPayload,
        channel.value(file("${projectDir}/modules/reports/build_report_site.py")),
        channel.value(file("${projectDir}/modules/reports/site.css")),
        channel.value(file("${projectDir}/modules/reports/site.js"))
    )
}
