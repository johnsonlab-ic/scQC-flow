nextflow.enable.dsl=2

include { MAPPING           } from './workflows/workflows'
include { SAMPLE_METADATA   } from './workflows/workflows'
include { AMBIENT           } from './workflows/workflows'
include { QC                } from './workflows/workflows'
include { HVG               } from './workflows/workflows'
include { INTEGRATION       } from './workflows/workflows'
include { ANNOTATION        } from './workflows/workflows'
include { ANNOTATION_METHODS } from './workflows/workflows'
include { ZOOMS             } from './workflows/workflows'
include { REPORT_SITE       } from './modules/reports/reports'
include { EXPORT_SCANPY     } from './modules/export/export_sc'
include { EXPORT_SEURAT     } from './modules/export/export_sc'
include { EXPORT_SCANPY_ZOOM } from './modules/export/export_sc'
include { EXPORT_SEURAT_ZOOM } from './modules/export/export_sc'
include { EXPORT_CELL_METADATA } from './modules/export/export_sc'
include { EXPORT_CELL_METADATA_ZOOM } from './modules/export/export_sc'

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
        --raw_data_dir      Parent directory (or comma-separated list of directories);
                            each subdirectory is one sample
                            (each subdir must contain *_R1_* and *_R2_* FASTQs)
        --cellrangerPath    Path to Cell Ranger install directory
                            (whitelists are extracted from lib/python/cellranger/barcodes/)
        --genome_fasta      Reference genome FASTA (for building simpleaf index)
        --genome_gtf        Reference genome GTF

    Optional:
        --chemistry         10x chemistry string (default: 3v4)
                            Accepted: 3v2, 3v3, 3v4, 3LT, 5v1, 5v2, 5v3, multiome
        --outputDir         Output directory (default: results)
        --publish_mapping_simpleaf
                            Publish large simpleaf mapping directories to
                            outputDir/mapping (default: false)
        --run_ambient       Run ambient RNA cleanup stage (default: true)
        --ambient_method    Ambient method (default: decontx). Accepted: decontx, cellbender
        --cellbender_env_name Preinstalled conda env name used to run CellBender (default: cellbender)
        --cellbender_mig_device MIG UUID for CellBender GPU slice (e.g. MIG-xxxx); empty = auto-detect
        --run_qc            Run cell-level QC (default: true; requires --run_ambient)
        --run_hvg           Run HVG selection (default: true; requires --run_qc)
        --run_integration   Run Harmony integration (default: true; requires --run_hvg)
        --run_annotation    Run pseudobulk marker-gene annotation (default: false;
                    requires --run_integration)
        --export            Export integrated objects (default: none).
                    Accepted: none, anndata, seurat, both
        --export_write_combined
                Also write combined export object(s) in addition to
                per-sample outputs (default: true)
        zoom                Nested zoom configuration map in nextflow.config.
            Each zoom can subset by integration cluster, legacy annotation label,
            or a label from a configured annotation method,
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
        annotation_methods  Optional list of method/reference annotation runs
                defined in config. Each item renders its own report and
                exports namespaced metadata columns.
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
    def rawDataDirs = params.raw_data_dir.toString().split(',').collect { rawDir -> rawDir.trim() }
    rawDataDirs.each { dir ->
        if (!file(dir).exists()) {
            error "--raw_data_dir does not exist: ${dir}"
        }
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
    def ambientMethod = (params.ambient_method ?: 'decontx').toString().trim().toLowerCase()
    if (!(ambientMethod in ['decontx', 'cellbender'])) {
        error "--ambient_method must be one of: decontx, cellbender"
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

    def rawAnnotationMethods = params.annotation_methods ?: []
    if (!(rawAnnotationMethods instanceof List)) {
        error "params.annotation_methods must be a list of annotation method definitions"
    }

    def normalizedAnnotationMethods = rawAnnotationMethods.collect { rawSpec ->
        if (!(rawSpec instanceof Map)) {
            error "Each entry in params.annotation_methods must be a map"
        }

        def methodIdRaw = rawSpec.id?.toString()?.trim()
        if (!methodIdRaw) {
            error "Each annotation method entry must define a non-empty 'id'"
        }
        def methodId = methodIdRaw.replaceAll(/[^A-Za-z0-9_]+/, '_')
        def engine = rawSpec.engine?.toString()?.trim()?.toLowerCase()
        if (!(engine in ['singler', 'xgboost'])) {
            error "Annotation method '${methodIdRaw}' has invalid engine '${engine}'. Use 'singler' or 'xgboost'."
        }

        def clusterCol = rawSpec.cluster_col?.toString()?.trim()
        if (!clusterCol && rawSpec.cluster_res != null) {
            clusterCol = "RNA_snn_res.${rawSpec.cluster_res.toString().trim()}"
        }

        def normalizedSpec = [
            id: methodId,
            engine: engine,
            reference_name: (rawSpec.reference_name ?: methodId).toString(),
            cluster_col: clusterCol ?: '',
        ]

        if (engine == 'singler') {
            if (!rawSpec.reference_rds) {
                error "Annotation method '${methodIdRaw}' requires 'reference_rds'"
            }
            if (!file(rawSpec.reference_rds).exists()) {
                error "Annotation method '${methodIdRaw}' reference_rds does not exist: ${rawSpec.reference_rds}"
            }
            if (!rawSpec.reference_label_col) {
                error "Annotation method '${methodIdRaw}' requires 'reference_label_col'"
            }
            normalizedSpec.reference_rds = rawSpec.reference_rds.toString()
            normalizedSpec.reference_label_col = rawSpec.reference_label_col.toString()
            normalizedSpec.fine_tune = rawSpec.containsKey('fine_tune') ? rawSpec.fine_tune : false
            normalizedSpec.prune = rawSpec.containsKey('prune') ? rawSpec.prune : true
            normalizedSpec.bp_type = (rawSpec.bp_type ?: 'multicore').toString()
        }

        if (engine == 'xgboost') {
            if (!rawSpec.model_rds) {
                error "Annotation method '${methodIdRaw}' requires 'model_rds'"
            }
            if (!file(rawSpec.model_rds).exists()) {
                error "Annotation method '${methodIdRaw}' model_rds does not exist: ${rawSpec.model_rds}"
            }
            if (!rawSpec.class_csv) {
                error "Annotation method '${methodIdRaw}' requires 'class_csv'"
            }
            if (!file(rawSpec.class_csv).exists()) {
                error "Annotation method '${methodIdRaw}' class_csv does not exist: ${rawSpec.class_csv}"
            }
            normalizedSpec.model_rds = rawSpec.model_rds.toString()
            normalizedSpec.class_csv = rawSpec.class_csv.toString()
            normalizedSpec.chunk_size = (rawSpec.chunk_size ?: 10000) as Integer
            normalizedSpec.scale_factor = (rawSpec.scale_factor ?: 10000) as BigDecimal
        }

        normalizedSpec
    }

    def annotationMethodIds = normalizedAnnotationMethods.collect { spec -> spec.id }
    if (annotationMethodIds.toSet().size() != annotationMethodIds.size()) {
        error "Annotation method ids must be unique: ${annotationMethodIds}"
    }
    if (normalizedAnnotationMethods && !params.run_integration) {
        error "params.annotation_methods requires --run_integration"
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
            if (!(zoomSource in ['cluster', 'annotation_label', 'annotation_method_label'])) {
                error "Zoom '${zoomName}' has invalid source '${zoomSource}'. Use 'cluster', 'annotation_label', or 'annotation_method_label'."
            }

            if (zoomSource == 'annotation_label' && !params.run_annotation) {
                error "Zoom '${zoomName}' uses source='annotation_label' and therefore requires --run_annotation"
            }

            def annotationMethodZoomId = null
            if (zoomSource == 'annotation_method_label') {
                annotationMethodZoomId = rawSpec.annotation_method_id?.toString()?.trim()
                if (!annotationMethodZoomId) {
                    error "Zoom '${zoomName}' uses source='annotation_method_label' and must define 'annotation_method_id'"
                }
                if (!(annotationMethodZoomId in annotationMethodIds)) {
                    error "Zoom '${zoomName}' references annotation_method_id='${annotationMethodZoomId}', but that id is not present in params.annotation_methods"
                }
            }

            def zoomValues = []
            def zoomLabel = null
            if (zoomSource in ['annotation_label', 'annotation_method_label']) {
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
            def integrationLeidenVals = integrationLeiden.tokenize(' ').findAll { value -> value?.trim() }
            def finestIntegrationRes = integrationLeidenVals ? integrationLeidenVals.max { a, b -> (new BigDecimal(a)) <=> (new BigDecimal(b)) } : '0.5'

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
                marker_sel_res: (rawSpec.marker_sel_res ?: zoomConfig.default_marker_sel_res ?: finestIntegrationRes).toString(),
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
                normalizedSpec.values = [zoomLabel]
                normalizedSpec.label = zoomLabel
                if (zoomSource == 'annotation_method_label') {
                    normalizedSpec.annotation_method_id = annotationMethodZoomId
                }
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
    ambient_method : ${ambientMethod}
    run_qc         : ${params.run_qc}
    run_hvg        : ${params.run_hvg}
    run_integration: ${params.run_integration}
    run_annotation : ${params.run_annotation}
    annotation_methods: ${normalizedAnnotationMethods ? normalizedAnnotationMethods.collect { spec -> spec.id }.join(', ') : 'none'}
    export         : ${exportMode}
    zooms          : ${normalizedZoomSpecs ? normalizedZoomSpecs.collect { spec -> spec.name }.join(', ') : 'disabled'}
    metadata_csv   : ${params.metadata_csv ?: 'not provided'}
    annotation_csv : ${params.annotation_marker_csv ?: 'not provided'}
    outputDir      : ${params.outputDir}
    ===================================
    """

    def landingStartValue = workflow.start ?: new Date()
    def landingStartedAt
    if (landingStartValue instanceof java.util.Date) {
        landingStartedAt = landingStartValue.toInstant()
            .atZone(java.time.ZoneId.systemDefault())
            .format(java.time.format.DateTimeFormatter.ISO_OFFSET_DATE_TIME)
    } else if (landingStartValue instanceof java.time.Instant) {
        landingStartedAt = landingStartValue
            .atZone(java.time.ZoneId.systemDefault())
            .format(java.time.format.DateTimeFormatter.ISO_OFFSET_DATE_TIME)
    } else if (landingStartValue instanceof java.time.OffsetDateTime) {
        landingStartedAt = landingStartValue.format(java.time.format.DateTimeFormatter.ISO_OFFSET_DATE_TIME)
    } else if (landingStartValue instanceof java.time.ZonedDateTime) {
        landingStartedAt = landingStartValue
            .toOffsetDateTime()
            .format(java.time.format.DateTimeFormatter.ISO_OFFSET_DATE_TIME)
    } else if (landingStartValue instanceof java.time.LocalDateTime) {
        landingStartedAt = landingStartValue
            .atZone(java.time.ZoneId.systemDefault())
            .format(java.time.format.DateTimeFormatter.ISO_OFFSET_DATE_TIME)
    } else {
        landingStartedAt = new Date().toInstant()
            .atZone(java.time.ZoneId.systemDefault())
            .format(java.time.format.DateTimeFormatter.ISO_OFFSET_DATE_TIME)
    }

    def landingPayload = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson([
            pipeline: 'scQC-flow',
            run_name: workflow.runName,
            profile: workflow.profile ?: 'default',
            started_at: landingStartedAt,
            generated_at: null,
            runtime_seconds: null,
            output_dir: params.outputDir,
            reports_dir: "${params.outputDir}/reports",
            pipeline_steps: [
                mapping: true,
                ambient: params.run_ambient,
                qc: params.run_ambient && params.run_qc,
                hvg: params.run_ambient && params.run_qc && params.run_hvg,
                integration: params.run_ambient && params.run_qc && params.run_hvg && params.run_integration,
                annotation: params.run_ambient && params.run_qc && params.run_hvg && params.run_integration && (params.run_annotation || normalizedAnnotationMethods),
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
                    ambient_method: ambientMethod,
                    run_qc: params.run_qc,
                    run_hvg: params.run_hvg,
                    run_integration: params.run_integration,
                    run_annotation: params.run_annotation,
                    annotation_methods: normalizedAnnotationMethods.collect { spec -> spec.id },
                    export: exportMode,
                ],
                cellbender: [
                    env_name: params.cellbender_env_name,
                    mig_device: params.cellbender_mig_device,
                    expected_cells_override: params.cellbender_expected_cells,
                    total_droplets_override: params.cellbender_total_droplets_included,
                    low_count_override: params.cellbender_low_count_threshold,
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
                    annotation_methods: normalizedAnnotationMethods,
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
        .fromPath(rawDataDirs.collect { rawDir -> "${rawDir}/*" }, type: 'dir')
        .map { dir -> tuple(dir.name, dir.name, dir.toString()) }

    def sample_metadata_ch
    def landing_integration_ch = channel.value(file("${projectDir}/templates/NO_FILE"))
    def annotation_cell_labels_ch = channel.value(file("${projectDir}/templates/NO_FILE"))
    def annotation_method_cell_labels_ch = channel.value(file("${projectDir}/templates/NO_FILE"))
    def annotation_export_metadata_ch = channel.empty()
    def hasAnnotationExportMetadata = false
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
                    AMBIENT.out.de_table,
                    AMBIENT.out.pb_empties
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
                    landing_integration_ch = INTEGRATION.out.integration_dt
                    report_pages = report_pages.mix(INTEGRATION.out.report)

                    if (params.run_annotation) {
                        ANNOTATION(
                            AMBIENT.out.h5_files,
                            INTEGRATION.out.integration_dt
                        )
                        annotation_cell_labels_ch = ANNOTATION.out.cell_labels
                        annotation_export_metadata_ch = annotation_export_metadata_ch.mix(ANNOTATION.out.export_metadata)
                        hasAnnotationExportMetadata = true
                        report_pages = report_pages.mix(ANNOTATION.out.report)
                    }

                    if (normalizedAnnotationMethods) {
                        ANNOTATION_METHODS(
                            AMBIENT.out.h5_files,
                            INTEGRATION.out.integration_dt,
                            normalizedAnnotationMethods
                        )
                        annotation_method_cell_labels_ch = ANNOTATION_METHODS.out.cell_labels.collect()
                        annotation_export_metadata_ch = annotation_export_metadata_ch.mix(ANNOTATION_METHODS.out.export_metadata)
                        hasAnnotationExportMetadata = true
                        report_pages = report_pages.mix(ANNOTATION_METHODS.out.report)
                    }

                    def annotation_export_input_ch = hasAnnotationExportMetadata
                        ? annotation_export_metadata_ch.collect()
                        : channel.value(file("${projectDir}/templates/NO_FILE"))

                    if (normalizedZoomSpecs) {
                        ZOOMS(
                            channel.fromList(normalizedZoomSpecs),
                            MAPPING.out.h5_files,
                            MAPPING.out.knee_data,
                            AMBIENT.out.h5_files,
                            QC.out.qc_metrics,
                            INTEGRATION.out.integration_dt,
                            annotation_cell_labels_ch,
                            annotation_method_cell_labels_ch
                        )
                        report_pages = report_pages.mix(ZOOMS.out.report)

                        if (exportMode != 'none') {
                            EXPORT_CELL_METADATA_ZOOM(
                                ZOOMS.out.zoom_int,
                                annotation_export_input_ch,
                                channel.value(file("${projectDir}/modules/export/export_cell_metadata.R")),
                                channel.value(file("${projectDir}/modules/export/export_utils.R"))
                            )
                        }

                        if (exportMode in ['anndata', 'both']) {
                            EXPORT_SCANPY_ZOOM(
                                ZOOMS.out.zoom_int,
                                AMBIENT.out.h5_files.map { _id, h5 -> h5 }.collect(),
                                annotation_export_input_ch,
                                channel.value(file(params.genome_gtf)),
                                channel.value(file("${projectDir}/modules/export/export_anndata.py"))
                            )
                        }

                        if (exportMode in ['seurat', 'both']) {
                            EXPORT_SEURAT_ZOOM(
                                ZOOMS.out.zoom_int,
                                AMBIENT.out.h5_files.map { _id, h5 -> h5 }.collect(),
                                annotation_export_input_ch,
                                channel.value(file(params.genome_gtf)),
                                channel.value(file("${projectDir}/modules/export/export_seurat.R")),
                                channel.value(file("${projectDir}/modules/export/export_utils.R"))
                            )
                        }
                    }

                    if (exportMode in ['anndata', 'both']) {
                        EXPORT_SCANPY(
                            AMBIENT.out.h5_files.map { _id, h5 -> h5 }.collect(),
                            QC.out.qc_metrics.map { _id, csv -> csv }.collect(),
                            INTEGRATION.out.integration_dt,
                            annotation_export_input_ch,
                            channel.value(file(params.genome_gtf)),
                            channel.value(file("${projectDir}/modules/export/export_anndata.py"))
                        )
                    }

                    if (exportMode != 'none') {
                        EXPORT_CELL_METADATA(
                            QC.out.qc_metrics.map { _id, csv -> csv }.collect(),
                            INTEGRATION.out.integration_dt,
                            annotation_export_input_ch,
                            channel.value(file("${projectDir}/modules/export/export_cell_metadata.R")),
                            channel.value(file("${projectDir}/modules/export/export_utils.R"))
                        )
                    }

                    if (exportMode in ['seurat', 'both']) {
                        EXPORT_SEURAT(
                            AMBIENT.out.h5_files.map { _id, h5 -> h5 }.collect(),
                            QC.out.qc_metrics.map { _id, csv -> csv }.collect(),
                            INTEGRATION.out.integration_dt,
                            annotation_export_input_ch,
                            channel.value(file(params.genome_gtf)),
                            channel.value(file("${projectDir}/modules/export/export_seurat.R")),
                            channel.value(file("${projectDir}/modules/export/export_utils.R"))
                        )
                    }
                }
            }
        }
    }

    def existingZoomReportsCh = channel.fromPath("${params.outputDir}/zoom/**/*.html", checkIfExists: false)
    report_pages = report_pages.mix(existingZoomReportsCh)

    def traceFileCh = params.containsKey('trace_report_suffix')
        ? channel.value(file("${params.outputDir}/pipeline_info/execution_trace_${params.trace_report_suffix}.txt"))
        : channel.value(file("${projectDir}/templates/NO_FILE"))

    REPORT_SITE(
        report_pages.collect(),
        channel.value(file("${projectDir}/modules/reports/landing_page.qmd")),
        landingPayload,
        landing_integration_ch,
        channel.value(file("${projectDir}/modules/reports/build_report_site.py")),
        channel.value(file("${projectDir}/modules/reports/site.css")),
        channel.value(file("${projectDir}/modules/reports/site.js")),
        traceFileCh
    )
}
