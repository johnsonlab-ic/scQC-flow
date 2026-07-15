// workflows.nf

include { SIMPLEAF_INDEX    } from '../modules/mapping/mapping.nf'
include { SIMPLEAF_QUANT    } from '../modules/mapping/mapping.nf'
include { BARCODE_ESTIMATION } from '../modules/mapping/mapping.nf'
include { MAPPING_REPORT    } from '../modules/reports/reports.nf'
include { DECONTX           } from '../modules/ambient/ambient.nf'
include { CELLBENDER        } from '../modules/ambient/ambient.nf'
include { AMBIENT_REPORT    } from '../modules/reports/reports.nf'
include { CELLBENDER_REPORT } from '../modules/reports/reports.nf'
include { CELL_CALLING      } from '../modules/cell_calling/cell_calling.nf'
include { CELL_CALLING_REPORT } from '../modules/reports/reports.nf'
include { EMPTYDROPS_CALLING } from '../modules/emptydrops/emptydrops.nf'
include { EMPTYDROPS_REPORT } from '../modules/reports/reports.nf'
include { KNEE_CALLING      } from '../modules/knee/knee.nf'
include { KNEE_REPORT       } from '../modules/reports/reports.nf'
include { CELLSWEEP         } from '../modules/cellsweep/cellsweep.nf'
include { CELLSWEEP_REPORT  } from '../modules/reports/reports.nf'
include { CELLSWEEP_TO_H5   } from '../modules/second_pass/second_pass.nf'
include { HVG_SELECTION_PER_SAMPLE } from '../modules/hvg/hvg.nf'
include { RUN_INTEGRATION_PER_SAMPLE } from '../modules/integration/integration.nf'
include { PER_SAMPLE_ANNOTATION_REPORT } from '../modules/reports/reports.nf'
include { STAGE_RAW_H5      } from '../modules/ambient_de/ambient_de.nf'
include { AMBIENT_DE as AMBIENT_DE_PROC } from '../modules/ambient_de/ambient_de.nf'
include { PREPARE_SAMPLE_METADATA } from '../modules/metadata/metadata.nf'
include { DOUBLET_DETECTION } from '../modules/qc/qc.nf'
include { APPLY_QC          } from '../modules/qc/qc.nf'
include { QC_REPORT         } from '../modules/reports/reports.nf'
include { HVG_SELECTION     } from '../modules/hvg/hvg.nf'
include { HVG_REPORT        } from '../modules/reports/reports.nf'
include { RUN_INTEGRATION   } from '../modules/integration/integration.nf'
include { INTEGRATION_REPORT } from '../modules/reports/reports.nf'
include { RUN_ANNOTATION_MARKERS } from '../modules/annotation/annotation.nf'
include { PREP_REPORT_INPUTS } from '../modules/annotation/annotation.nf'
include { FORMAT_ANNOTATION_EXPORT_METADATA } from '../modules/annotation/annotation.nf'
include { PREPARE_ANNOTATION_QUERY } from '../modules/annotation/annotation.nf'
include { RUN_SINGLER_REFERENCE_ANNOTATION } from '../modules/annotation/annotation.nf'
include { RUN_XGBOOST_REFERENCE_ANNOTATION } from '../modules/annotation/annotation.nf'
include { COMBINE_ANNOTATION_METHOD_OUTPUTS } from '../modules/annotation/annotation.nf'
include { ANNOTATION_REPORT } from '../modules/reports/reports.nf'
include { ANNOTATION_METHOD_REPORT } from '../modules/reports/reports.nf'
include { PREPARE_ZOOM_SUBSET } from '../modules/zoom/zoom.nf'
include { STAGE_ZOOM_RAW_H5 } from '../modules/zoom/zoom.nf'
include { ZOOM_AMBIENT_DE } from '../modules/zoom/zoom.nf'
include { ZOOM_HVG_SELECTION } from '../modules/zoom/zoom.nf'
include { RUN_ZOOM_INTEGRATION } from '../modules/zoom/zoom.nf'
include { RUN_ZOOM_MARKERS } from '../modules/zoom/zoom.nf'
include { ZOOM_REPORT } from '../modules/zoom/zoom.nf'

// =============================================================================
// MAPPING
//
// FASTQ -> simpleaf quant -> barcode estimation -> report
//
// Chemistry is provided by the user (params.chemistry, default '3v4').
// Lookup tables map it to the AF chemistry string, orientation,
// and the Cell Ranger whitelist filename.
//
// Required params (set in nextflow.config or on CLI):
//   --raw_data_dir     parent dir; each subdir is one sample
//   --cellrangerPath   Cell Ranger install directory (for whitelists)
//   --genome_fasta     reference genome FASTA
//   --genome_gtf       reference genome GTF
//   --chemistry        10x chemistry string (default: 3v4)
//   --outputDir
// =============================================================================

workflow MAPPING {

    take:
    samples_ch    // tuple(sampleId, sampleName, fastqPath)

    main:

    // ------------------------------------------------------------------
    // 1. Resolve chemistry params (validated early — fails immediately
    //    if the user passed an unrecognised chemistry string)
    // ------------------------------------------------------------------
    def CHEMISTRY_TO_AF = [
        '3v2': '10xv2', '5v1': '10xv2', '5v2': '10xv2',
        '3v3': '10xv3', '5v3': '10xv3', '3v4': '10xv3',
        '3LT': '10xv3', 'multiome': '10xv3',
    ]

    def CHEMISTRY_TO_ORI = [
        '3v2': 'fw', '5v1': 'rc', '5v2': 'rc',
        '3v3': 'fw', '5v3': 'rc', '3v4': 'fw',
        '3LT': 'fw', 'multiome': 'fw',
    ]

    def CHEMISTRY_TO_WHITELIST = [
        '3v2': '737K-august-2016.txt',
        '5v1': '737K-august-2016.txt',
        '5v2': '737K-august-2016.txt',
        '3v3': '3M-february-2018_TRU.txt.gz',
        '5v3': '3M-5pgex-jan-2023.txt.gz',
        '3v4': '3M-3pgex-may-2023_TRU.txt.gz',
        '3LT': '9K-LT-march-2021.txt.gz',
        'multiome': '737K-arc-v1.txt.gz',
    ]

    if (!CHEMISTRY_TO_AF.containsKey(params.chemistry)) {
        error "Unknown chemistry '${params.chemistry}'. Accepted values: ${CHEMISTRY_TO_AF.keySet().join(', ')}"
    }
    def af_chemistry       = CHEMISTRY_TO_AF[params.chemistry]
    def orientation        = CHEMISTRY_TO_ORI[params.chemistry]
    def whitelist_filename = CHEMISTRY_TO_WHITELIST[params.chemistry]

    // ------------------------------------------------------------------
    // 2. Build simpleaf index from genome FASTA + GTF
    // ------------------------------------------------------------------
    SIMPLEAF_INDEX(
        channel.value(file(params.genome_fasta)),
        channel.value(file(params.genome_gtf))
    )

    // ------------------------------------------------------------------
    // 3. simpleaf quant (per sample)
    // ------------------------------------------------------------------
    SIMPLEAF_QUANT(
        samples_ch,
        af_chemistry,
        orientation,
        whitelist_filename,
        params.cellrangerPath,
        SIMPLEAF_INDEX.out.index_dir,
        SIMPLEAF_INDEX.out.t2g
    )

    // ------------------------------------------------------------------
    // 4. Barcode estimation (per sample)
    //    Produces H5 + knee CSV + ambient params
    // ------------------------------------------------------------------
    BARCODE_ESTIMATION(
        SIMPLEAF_QUANT.out.quant_dirs,
        channel.value(file("${projectDir}/modules/mapping/barcode_estimation.R"))
    )

    // ------------------------------------------------------------------
    // 5. Mapping report
    // ------------------------------------------------------------------
    MAPPING_REPORT(
        BARCODE_ESTIMATION.out.knee_data.map { _id, csv -> csv }.collect(),
        channel.value(file("${projectDir}/modules/reports/mapping_report.qmd"))
    )

    emit:
    h5_files        = BARCODE_ESTIMATION.out.h5_files
    knee_data       = BARCODE_ESTIMATION.out.knee_data
    ambient_params  = BARCODE_ESTIMATION.out.ambient_params
    alevinfry_stats = BARCODE_ESTIMATION.out.alevinfry_stats
    report          = MAPPING_REPORT.out.html
}

// =============================================================================
// AMBIENT  (dev-cellsweep: cell calling, no ambient removal)
//
// H5 + knee CSV -> CELL_CALLING (splice-aware GMM + guardrails) -> report
//
// decontX/CellBender are dropped on this branch: cells flow downstream
// UNCORRECTED (CellSweep decontamination is evaluated post-hoc). CELL_CALLING
// emits the is_cell matrix under the same filt_counts_<id>.h5 filename contract
// the QC/HVG stages already expect, plus the is_empty barcode set used to build
// the ambient profile for AMBIENT_DE (and, later, CellSweep).
// =============================================================================

workflow AMBIENT {

    take:
    h5_ch             // tuple(sampleId, h5_file)      from MAPPING.h5_files
    knee_ch           // tuple(sampleId, knee_csv)      from MAPPING.knee_data
    alevinfry_stats_ch // tuple(sampleId, stats_csv)    from MAPPING.alevinfry_stats

    main:

    joined_input = h5_ch.join(knee_ch)   // tuple(sampleId, h5, knee_csv)

    if (params.cell_caller == 'emptydrops') {
        EMPTYDROPS_CALLING(
            h5_ch,
            channel.value(file("${projectDir}/modules/emptydrops/emptydrops_calling.R"))
        )
        EMPTYDROPS_REPORT(
            EMPTYDROPS_CALLING.out.summaries.map { _id, csv -> csv }.collect(),
            EMPTYDROPS_CALLING.out.labels.map    { _id, csv -> csv }.collect(),
            channel.value(file("${projectDir}/modules/reports/emptydrops_report.qmd"))
        )
        cc_h5       = EMPTYDROPS_CALLING.out.h5_files
        cc_barcodes = EMPTYDROPS_CALLING.out.barcodes
        cc_empties  = EMPTYDROPS_CALLING.out.empties
        cc_labels   = EMPTYDROPS_CALLING.out.labels
        cc_report   = EMPTYDROPS_REPORT.out.html
    } else if (params.cell_caller == 'knee') {
        KNEE_CALLING(
            joined_input,
            channel.value(file("${projectDir}/modules/knee/knee_calling.R"))
        )
        KNEE_REPORT(
            knee_ch.map { _id, csv -> csv }.collect(),
            alevinfry_stats_ch.map { _id, csv -> csv }.collect(),
            KNEE_CALLING.out.summaries.map { _id, csv -> csv }.collect(),
            KNEE_CALLING.out.labels.map    { _id, csv -> csv }.collect(),
            channel.value(file("${projectDir}/modules/reports/knee_report.qmd"))
        )
        cc_h5       = KNEE_CALLING.out.h5_files
        cc_barcodes = KNEE_CALLING.out.barcodes
        cc_empties  = KNEE_CALLING.out.empties
        cc_labels   = KNEE_CALLING.out.labels
        cc_report   = KNEE_REPORT.out.html
    } else {
        CELL_CALLING(
            joined_input,
            channel.value(file("${projectDir}/modules/cell_calling/cell_calling.R"))
        )
        CELL_CALLING_REPORT(
            alevinfry_stats_ch.map { _id, csv -> csv }.collect(),
            CELL_CALLING.out.summaries.map { _id, csv -> csv }.collect(),
            CELL_CALLING.out.labels.map    { _id, csv -> csv }.collect(),
            CELL_CALLING.out.gmm.map       { _id, rds -> rds }.collect(),
            channel.value(file("${projectDir}/modules/reports/cell_calling_report.qmd"))
        )
        cc_h5       = CELL_CALLING.out.h5_files
        cc_barcodes = CELL_CALLING.out.barcodes
        cc_empties  = CELL_CALLING.out.empties
        cc_labels   = CELL_CALLING.out.labels
        cc_report   = CELL_CALLING_REPORT.out.html
    }

    STAGE_RAW_H5(h5_ch)

    AMBIENT_DE_PROC(
        STAGE_RAW_H5.out.h5.collect(),
        knee_ch.map { _id, csv -> csv }.collect(),
        cc_h5.map      { _id, h5 -> h5 }.collect(),
        cc_empties.map { _id, csv -> csv }.collect(),
        channel.value(file(params.genome_gtf)),
        channel.value(file("${projectDir}/modules/ambient_de/ambient_de.R"))
    )

    emit:
    h5_files   = cc_h5
    barcodes   = cc_barcodes
    empties    = cc_empties
    labels     = cc_labels
    report     = cc_report
    de_table   = AMBIENT_DE_PROC.out.de_table
    pb_empties = AMBIENT_DE_PROC.out.pb_empties
}

// =============================================================================
// PER_SAMPLE_ANNOTATION  (workflow_improvement.md step 4)
//
// Per sample, ahead of any cross-sample pooling: normalize -> HVG -> PCA ->
// Leiden(cellsweep_celltype_col resolution) -> XGBoost annotation. Feeds
// CellSweep's per-cell "celltype" grouping. Runs only when
// params.ambient_method == 'cellsweep'.
// =============================================================================

workflow PER_SAMPLE_ANNOTATION_WF {

    take:
    h5_ch         // tuple(sampleId, filt_counts.h5)   from AMBIENT.h5_files
    qc_metrics_ch // tuple(sampleId, qc_metrics.csv.gz) from QC.qc_metrics
    de_table      // path edger_dt.csv.gz               from AMBIENT.de_table
    xgb_spec      // normalized xgboost annotation method spec (cluster_col already set)

    main:
    hvg_input_ch = h5_ch.join(qc_metrics_ch)   // tuple(sid, h5, qc_csv)

    HVG_SELECTION_PER_SAMPLE(
        hvg_input_ch,
        channel.value(file(params.genome_gtf)),
        de_table,
        channel.value(file("${projectDir}/modules/hvg/hvg_selection.py"))
    )

    integration_input_ch = HVG_SELECTION_PER_SAMPLE.out.hvg_counts
        .join(HVG_SELECTION_PER_SAMPLE.out.dbl_hvg_counts)
        .join(qc_metrics_ch)   // tuple(sid, hvg_counts, dbl_hvg_counts, qc_csv)

    RUN_INTEGRATION_PER_SAMPLE(
        integration_input_ch,
        channel.value(file("${projectDir}/modules/integration/run_integration.py"))
    )

    query_input_ch = h5_ch.join(RUN_INTEGRATION_PER_SAMPLE.out.integration_dt)   // tuple(sid, h5, integration_dt)

    PREPARE_ANNOTATION_QUERY(
        query_input_ch,
        channel.value(file("${projectDir}/modules/annotation/prepare_annotation_query.R")),
        channel.value(file("${projectDir}/modules/export/export_utils.R"))
    )

    def spec_b64 = groovy.json.JsonOutput.toJson(xgb_spec).bytes.encodeBase64().toString()
    xgb_run_inputs_ch = PREPARE_ANNOTATION_QUERY.out.query.map { sample_id, query_rds ->
        tuple(xgb_spec.id.toString(), spec_b64, sample_id, query_rds, file(xgb_spec.model_rds.toString()), file(xgb_spec.class_csv.toString()))
    }

    RUN_XGBOOST_REFERENCE_ANNOTATION(
        xgb_run_inputs_ch,
        channel.value(file("${projectDir}/modules/annotation/run_xgboost_reference_annotation.R"))
    )

    PER_SAMPLE_ANNOTATION_REPORT(
        RUN_XGBOOST_REFERENCE_ANNOTATION.out.result.map { _mid, _sb64, _sid, cells_csv, _cl, _exp -> cells_csv }.collect(),
        RUN_XGBOOST_REFERENCE_ANNOTATION.out.result.map { _mid, _sb64, _sid, _cells, cluster_csv, _exp -> cluster_csv }.collect(),
        qc_metrics_ch.map { _id, csv -> csv }.collect(),
        channel.value(file("${projectDir}/modules/reports/per_sample_annotation_report.qmd")),
        channel.value(file("${projectDir}/modules/integration/integration_plots.R")),
        channel.value(file("${projectDir}/modules/qc/qc_plots.R")),
        xgb_spec.id.toString(),
        xgb_spec.cluster_col.toString()
    )

    emit:
    integration_dt   = RUN_INTEGRATION_PER_SAMPLE.out.integration_dt
    cluster_summary  = RUN_XGBOOST_REFERENCE_ANNOTATION.out.result.map { _mid, _sb64, sample_id, _cells, cluster_csv, _exp -> tuple(sample_id, cluster_csv) }
    annotation_cells = RUN_XGBOOST_REFERENCE_ANNOTATION.out.result.map { _mid, _sb64, sample_id, cells_csv, _cl, _exp -> tuple(sample_id, cells_csv) }
    report           = PER_SAMPLE_ANNOTATION_REPORT.out.html
}

// =============================================================================
// CELLSWEEP  (per-sample ambient decontamination)
//
// Per sample: is_cell counts + raw H5 (for the is_empty ambient profile) + that
// sample's own PER_SAMPLE_ANNOTATION_WF cluster labels -> CellSweep EM ->
// decontaminated counts + alpha_hat -> converted back to the filt_counts.h5 /
// qc_metrics.csv.gz contract the shared HVG/INTEGRATION/ANNOTATION* workflows
// expect. Fully self-contained: callers just swap AMBIENT.h5_files/QC.qc_metrics
// for CELLSWEEP_WF.h5_files/qc_metrics downstream.
// =============================================================================

workflow CELLSWEEP_WF {

    take:
    filt_ch           // tuple(sampleId, filt_counts.h5)    from AMBIENT.h5_files
    raw_ch            // tuple(sampleId, af_counts_mat.h5)  from MAPPING.h5_files
    empties_ch        // tuple(sampleId, empty_barcodes.csv) from AMBIENT.empties
    integration_dt_ch // tuple(sampleId, integration_dt.csv.gz)   from PER_SAMPLE_ANNOTATION_WF.integration_dt
    cluster_ann_ch    // tuple(sampleId, cluster_summary.csv.gz)  from PER_SAMPLE_ANNOTATION_WF.cluster_summary
    annotation_cells_ch // tuple(sampleId, annotation_cells.csv.gz) from PER_SAMPLE_ANNOTATION_WF.annotation_cells
    labels_ch         // tuple(sampleId, <caller>_labels.csv.gz) from AMBIENT.labels
    qc_metrics_ch     // tuple(sampleId, qc_metrics.csv.gz) from QC.qc_metrics
    xgb_method_id     // val — id of the xgboost annotation method used by PER_SAMPLE_ANNOTATION_WF
    xgb_cluster_col   // val — its cluster_col (params.cellsweep_celltype_col)

    main:
    in_ch = filt_ch.join(raw_ch).join(empties_ch).join(integration_dt_ch).join(cluster_ann_ch)

    CELLSWEEP(
        in_ch,
        channel.value(file("${projectDir}/modules/cellsweep/cellsweep_run_sample.py"))
    )

    CELLSWEEP_TO_H5(
        CELLSWEEP.out.h5ad.join(qc_metrics_ch),
        channel.value(file("${projectDir}/modules/second_pass/cellsweep_to_h5.py"))
    )

    CELLSWEEP_REPORT(
        labels_ch.map            { _id, csv -> csv }.collect(),
        qc_metrics_ch.map        { _id, csv -> csv }.collect(),
        CELLSWEEP.out.alpha.map  { _id, csv -> csv }.collect(),
        annotation_cells_ch.map  { _id, csv -> csv }.collect(),
        CELLSWEEP_TO_H5.out.qc_metrics.map { _id, csv -> csv }.collect(),
        channel.value(file("${projectDir}/modules/reports/cellsweep_report.qmd")),
        channel.value(file("${projectDir}/modules/integration/integration_plots.R")),
        xgb_method_id,
        xgb_cluster_col
    )

    emit:
    alpha      = CELLSWEEP.out.alpha
    h5_files   = CELLSWEEP_TO_H5.out.h5_files    // tuple(sampleId, filt_counts.h5) — corrected counts
    qc_metrics = CELLSWEEP_TO_H5.out.qc_metrics  // tuple(sampleId, qc_metrics.csv.gz) — carries alpha_hat
    report     = CELLSWEEP_REPORT.out.html
}

// =============================================================================
// SAMPLE_METADATA
//
// Normalize and validate sample metadata once, keyed by sample_id.
// =============================================================================

workflow SAMPLE_METADATA {

    take:
    sample_ids_ch

    main:

    PREPARE_SAMPLE_METADATA(
        sample_ids_ch,
        channel.value(file(params.metadata_csv)),
        channel.value(file("${projectDir}/modules/metadata/prepare_metadata.py"))
    )

    emit:
    PREPARE_SAMPLE_METADATA.out.sample_metadata
}

// =============================================================================
// QC
//
// Filtered H5 -> per-sample QC (metrics + doublet detection) -> QC report
//
// Takes ambient-cleaned filtered H5 files from AMBIENT and the genome GTF.
// Runs SAMPLE_QC per sample (scDblFinder, hard + soft thresholds),
// then renders a combined QC report.
// =============================================================================

workflow QC {

    take:
    h5_ch        // tuple(sampleId, h5_file)  from AMBIENT.h5_files
    sample_meta_ch // path sample_metadata.csv.gz or NO_FILE placeholder

    main:

    // ------------------------------------------------------------------
    // 1. Doublet detection (per sample) — cached unless H5 changes
    // ------------------------------------------------------------------
    DOUBLET_DETECTION(
        h5_ch,
        channel.value(file(params.genome_gtf)),
        channel.value(file("${projectDir}/modules/qc/doublet_detection.R"))
    )

    // ------------------------------------------------------------------
    // 2. QC metrics + threshold filtering (per sample)
    //    Re-runs on threshold changes; doublet detection is unaffected.
    // ------------------------------------------------------------------
    APPLY_QC(
        h5_ch.join(DOUBLET_DETECTION.out.dbl_results),
        sample_meta_ch,
        channel.value(file(params.genome_gtf)),
        channel.value(file("${projectDir}/modules/qc/apply_qc.R"))
    )

    // ------------------------------------------------------------------
    // 3. QC report (one HTML across all samples)
    // ------------------------------------------------------------------
    QC_REPORT(
        APPLY_QC.out.qc_metrics.map { _id, csv -> csv }.collect(),
        channel.value(file("${projectDir}/modules/qc/qc_plots.R")),
        channel.value(file("${projectDir}/modules/reports/qc_report.qmd"))
    )

    emit:
    qc_metrics  = APPLY_QC.out.qc_metrics
    qc_summary  = APPLY_QC.out.qc_summary
    dbl_results = DOUBLET_DETECTION.out.dbl_results
    report      = QC_REPORT.out.html
}

// =============================================================================
// HVG_SELECTION
//
// Filtered H5 files + QC metrics -> HVG selection -> HVG report
//
// Collects all per-sample filtered H5 files from AMBIENT and QC metrics from
// QC. Sums S+U+A counts, filters to QC-passing singlets, runs scanpy HVG
// selection (seurat_v3, batch_key=sample_id), and saves:
//   - hvg_stats.csv.gz    (per-gene mean/variance/HVG flag)
//   - hvg_counts.h5       (raw SUA-summed counts for HVG genes only)
// =============================================================================

workflow HVG {

    take:
    h5_ch         // tuple(sampleId, h5_file)      from AMBIENT.h5_files
    qc_metrics_ch // tuple(sampleId, csv_gz)        from QC.qc_metrics
    de_table      // path edger_dt.csv.gz           from AMBIENT.de_table
    pb_empties    // path pb_empties.rds            from AMBIENT.pb_empties

    main:

    HVG_SELECTION(
        h5_ch.map { _id, h5 -> h5 }.collect(),
        qc_metrics_ch.map { _id, csv -> csv }.collect(),
        channel.value(file(params.genome_gtf)),
        de_table,
        channel.value(file("${projectDir}/modules/hvg/hvg_selection.py"))
    )

    HVG_REPORT(
        HVG_SELECTION.out.hvg_stats,
        de_table,
        pb_empties,
        channel.value(file("${projectDir}/modules/reports/hvg_report.qmd")),
        channel.value(file("${projectDir}/modules/hvg/hvg_plots.R"))
    )

    emit:
    hvg_counts     = HVG_SELECTION.out.hvg_counts
    dbl_hvg_counts = HVG_SELECTION.out.dbl_hvg_counts
    hvg_stats      = HVG_SELECTION.out.hvg_stats
    report         = HVG_REPORT.out.html
}

// =============================================================================
// INTEGRATION
//
// HVG count matrix + QC metadata -> Harmony integration -> integration report
//
// Loads the HVG count matrix from HVG, reads metadata already attached to
// QC metrics, normalises (CPM10k + log1p),
// scales, runs PCA -> Harmony -> neighbors -> Leiden (multiple resolutions)
// -> UMAP, and saves integration_dt.csv.gz.
// =============================================================================

workflow INTEGRATION {

    take:
    hvg_counts_ch      // path  hvg_counts.h5 (singlets) from HVG
    dbl_hvg_counts_ch  // path  dbl_hvg_counts.h5 (doublets) from HVG
    qc_metrics_ch      // tuple(sampleId, csv_gz) from QC.qc_metrics

    main:

    RUN_INTEGRATION(
        hvg_counts_ch,
        dbl_hvg_counts_ch,
        qc_metrics_ch.map { _id, csv -> csv }.collect(),
        channel.value(file("${projectDir}/modules/integration/run_integration.py"))
    )

    INTEGRATION_REPORT(
        RUN_INTEGRATION.out.integration_dt,
        RUN_INTEGRATION.out.dbl_sweep,
        qc_metrics_ch.map { _id, csv -> csv }.collect(),
        channel.value(file("${projectDir}/modules/reports/integration_report.qmd")),
        channel.value(file("${projectDir}/modules/integration/integration_plots.R"))
    )

    emit:
    integration_dt = RUN_INTEGRATION.out.integration_dt
    report         = INTEGRATION_REPORT.out.html
}

// =============================================================================
// ANNOTATION
//
// Filtered H5 files + integration output + user marker panel -> pseudobulk
// marker DE -> annotation report.
//
// Uses a single selected clustering resolution from the integration output,
// aggregates counts per sample within each cluster, runs one-vs-rest edgeR
// marker testing, and renders heatmaps for the provided marker panel together
// with top marker plots for each cluster.
// =============================================================================

workflow ANNOTATION {

    take:
    h5_ch             // tuple(sampleId, h5_file) from AMBIENT.h5_files
    integration_dt_ch // path integration_dt.csv.gz from INTEGRATION

    main:

    RUN_ANNOTATION_MARKERS(
        h5_ch.map { _id, h5 -> h5 }.collect(),
        integration_dt_ch,
        channel.value(file(params.genome_gtf)),
        channel.value(file(params.annotation_marker_csv)),
        channel.value(file("${projectDir}/modules/annotation/marker_genes.R")),
        channel.value(file("${projectDir}/modules/annotation/annotation_utils.R"))
    )

    PREP_REPORT_INPUTS(
        h5_ch.map { _id, h5 -> h5 }.collect(),
        integration_dt_ch,
        RUN_ANNOTATION_MARKERS.out.markers,
        channel.value(file("${projectDir}/modules/annotation/prep_report_inputs.R")),
        channel.value(file("${projectDir}/modules/annotation/annotation_utils.R"))
    )

    FORMAT_ANNOTATION_EXPORT_METADATA(
        RUN_ANNOTATION_MARKERS.out.cell_labels,
        channel.value(file("${projectDir}/modules/annotation/format_annotation_export_metadata.R"))
    )

    ANNOTATION_REPORT(
        integration_dt_ch,
        RUN_ANNOTATION_MARKERS.out.markers,
        RUN_ANNOTATION_MARKERS.out.logcpms,
        RUN_ANNOTATION_MARKERS.out.marker_panel,
        RUN_ANNOTATION_MARKERS.out.marker_expr,
        PREP_REPORT_INPUTS.out.top_marker_expr,
        RUN_ANNOTATION_MARKERS.out.cell_labels,
        channel.value(file("${projectDir}/modules/reports/annotation_report.qmd")),
        channel.value(file("${projectDir}/modules/annotation/annotation_utils.R")),
        channel.value(file("${projectDir}/modules/integration/integration_plots.R"))
    )

    emit:
    markers = RUN_ANNOTATION_MARKERS.out.markers
    logcpms = RUN_ANNOTATION_MARKERS.out.logcpms
    marker_panel = RUN_ANNOTATION_MARKERS.out.marker_panel
    marker_expr = RUN_ANNOTATION_MARKERS.out.marker_expr
    top_marker_expr = PREP_REPORT_INPUTS.out.top_marker_expr
    cell_labels = RUN_ANNOTATION_MARKERS.out.cell_labels
    export_metadata = FORMAT_ANNOTATION_EXPORT_METADATA.out.export_metadata
    report = ANNOTATION_REPORT.out.html
}

// =============================================================================
// ANNOTATION_METHODS
//
// Per-sample prepared query SCEs -> method/reference-specific annotation runs
// -> combined method outputs -> one HTML report per configured annotation method.
// =============================================================================

workflow ANNOTATION_METHODS {

    take:
    h5_ch
    integration_dt_ch
    annotation_methods

    main:

    PREPARE_ANNOTATION_QUERY(
        h5_ch.combine(integration_dt_ch),
        channel.value(file("${projectDir}/modules/annotation/prepare_annotation_query.R")),
        channel.value(file("${projectDir}/modules/export/export_utils.R"))
    )

    def singler_methods = annotation_methods.findAll { spec -> spec.engine == 'singler' }
    def xgboost_methods = annotation_methods.findAll { spec -> spec.engine == 'xgboost' }

    def annotation_sample_results_ch = channel.empty()

    if (singler_methods) {
        def singler_specs_ch = channel.fromList(singler_methods)
            .map { spec ->
                def spec_json = groovy.json.JsonOutput.toJson(spec)
                tuple(spec.id.toString(), spec_json.bytes.encodeBase64().toString(), file(spec.reference_rds.toString()))
            }

        def singler_run_inputs_ch = singler_specs_ch
            .combine(PREPARE_ANNOTATION_QUERY.out.query)
            .map { joined -> tuple(joined[0], joined[1], joined[3], joined[4], joined[2]) }

        RUN_SINGLER_REFERENCE_ANNOTATION(
            singler_run_inputs_ch,
            channel.value(file("${projectDir}/modules/annotation/run_singler_reference_annotation.R"))
        )

        annotation_sample_results_ch = annotation_sample_results_ch.mix(RUN_SINGLER_REFERENCE_ANNOTATION.out.result)
    }

    if (xgboost_methods) {
        def xgboost_specs_ch = channel.fromList(xgboost_methods)
            .map { spec ->
                def spec_json = groovy.json.JsonOutput.toJson(spec)
                tuple(
                    spec.id.toString(),
                    spec_json.bytes.encodeBase64().toString(),
                    file(spec.model_rds.toString()),
                    file(spec.class_csv.toString())
                )
            }

        def xgboost_run_inputs_ch = xgboost_specs_ch
            .combine(PREPARE_ANNOTATION_QUERY.out.query)
            .map { joined -> tuple(joined[0], joined[1], joined[4], joined[5], joined[2], joined[3]) }

        RUN_XGBOOST_REFERENCE_ANNOTATION(
            xgboost_run_inputs_ch,
            channel.value(file("${projectDir}/modules/annotation/run_xgboost_reference_annotation.R"))
        )

        annotation_sample_results_ch = annotation_sample_results_ch.mix(RUN_XGBOOST_REFERENCE_ANNOTATION.out.result)
    }

    def annotation_method_inputs_ch = annotation_sample_results_ch
        .map { method_id, spec_b64, _sample_id, cells_csv, _cluster_csv, export_csv ->
            tuple(method_id, spec_b64, cells_csv, export_csv)
        }
        .groupTuple(by: [0, 1])

    COMBINE_ANNOTATION_METHOD_OUTPUTS(
        annotation_method_inputs_ch,
        channel.value(file("${projectDir}/modules/annotation/combine_annotation_method_outputs.R"))
    )

    ANNOTATION_METHOD_REPORT(
        COMBINE_ANNOTATION_METHOD_OUTPUTS.out.result,
        channel.value(file("${projectDir}/modules/reports/annotation_method_report.qmd"))
    )

    emit:
    cell_labels = COMBINE_ANNOTATION_METHOD_OUTPUTS.out.result.map { _id, _spec_b64, cells_csv, _cluster_csv, _export_csv -> cells_csv }
    export_metadata = COMBINE_ANNOTATION_METHOD_OUTPUTS.out.result.map { _id, _spec_b64, _cells_csv, _cluster_csv, export_csv -> export_csv }
    cluster_summary = COMBINE_ANNOTATION_METHOD_OUTPUTS.out.result.map { _id, _spec_b64, _cells_csv, cluster_csv, _export_csv -> cluster_csv }
    report = ANNOTATION_METHOD_REPORT.out.html
}

// =============================================================================
// ZOOMS
//
// Post-integration subset reruns driven either by cluster identity or by
// broad annotation labels. Each zoom subsets QC rows, reruns HVG selection and
// integration on that subset, then computes top marker genes per cluster and
// renders a dedicated report.
// =============================================================================

workflow ZOOMS {

    take:
    zoom_specs_ch
    raw_h5_ch
    knee_ch
    h5_ch
    empty_bc_ch
    qc_metrics_ch
    integration_dt_ch
    annotation_cell_labels_ch
    annotation_method_cell_labels_ch

    main:

    def encoded_zoom_specs_ch = zoom_specs_ch.map { spec ->
        def spec_json = groovy.json.JsonOutput.toJson(spec)
        tuple(spec.name.toString(), spec_json.bytes.encodeBase64().toString())
    }

    PREPARE_ZOOM_SUBSET(
        encoded_zoom_specs_ch,
        integration_dt_ch,
        qc_metrics_ch.map { _id, csv -> csv }.collect(),
        annotation_cell_labels_ch,
        annotation_method_cell_labels_ch,
        channel.value(file("${projectDir}/modules/zoom/prepare_zoom_subset.py"))
    )

    STAGE_ZOOM_RAW_H5(raw_h5_ch)

    ZOOM_AMBIENT_DE(
        PREPARE_ZOOM_SUBSET.out.zoom_inputs,
        STAGE_ZOOM_RAW_H5.out.h5.collect(),
        knee_ch.map { _id, csv -> csv }.collect(),
        h5_ch.map { _id, h5 -> h5 }.collect(),
        empty_bc_ch.map { _id, csv -> csv }.collect(),
        channel.value(file(params.genome_gtf)),
        channel.value(file("${projectDir}/modules/ambient_de/ambient_de.R"))
    )

    ZOOM_HVG_SELECTION(
        ZOOM_AMBIENT_DE.out.zoom_ambient,
        h5_ch.map { _id, h5 -> h5 }.collect(),
        channel.value(file(params.genome_gtf)),
        channel.value(file("${projectDir}/modules/hvg/hvg_selection.py"))
    )

    RUN_ZOOM_INTEGRATION(
        ZOOM_HVG_SELECTION.out.zoom_hvg,
        channel.value(file("${projectDir}/modules/integration/run_integration.py"))
    )

    RUN_ZOOM_MARKERS(
        RUN_ZOOM_INTEGRATION.out.zoom_int,
        h5_ch.map { _id, h5 -> h5 }.collect(),
        channel.value(file(params.genome_gtf)),
        channel.value(file("${projectDir}/modules/zoom/zoom_markers.R")),
        channel.value(file("${projectDir}/modules/annotation/annotation_utils.R"))
    )

    ZOOM_REPORT(
        RUN_ZOOM_MARKERS.out.zoom_markers,
        integration_dt_ch,
        channel.value(file("${projectDir}/modules/reports/zoom_report.qmd")),
        channel.value(file("${projectDir}/modules/integration/integration_plots.R")),
        channel.value(file("${projectDir}/modules/annotation/annotation_utils.R"))
    )

    emit:
    report = ZOOM_REPORT.out.html
    zoom_int = RUN_ZOOM_INTEGRATION.out.zoom_int
}
