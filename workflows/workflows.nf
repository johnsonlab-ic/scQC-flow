// workflows.nf

include { SIMPLEAF_INDEX    } from '../modules/mapping/mapping.nf'
include { SIMPLEAF_QUANT    } from '../modules/mapping/mapping.nf'
include { BARCODE_ESTIMATION } from '../modules/mapping/mapping.nf'
include { MAPPING_REPORT    } from '../modules/reports/reports.nf'
include { DECONTX           } from '../modules/ambient/ambient.nf'
include { CELLBENDER        } from '../modules/ambient/ambient.nf'
include { AMBIENT_REPORT    } from '../modules/reports/reports.nf'
include { CELLBENDER_REPORT } from '../modules/reports/reports.nf'
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
include { ANNOTATION_REPORT } from '../modules/reports/reports.nf'
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
    h5_files       = BARCODE_ESTIMATION.out.h5_files
    knee_data      = BARCODE_ESTIMATION.out.knee_data
    ambient_params = BARCODE_ESTIMATION.out.ambient_params
    report         = MAPPING_REPORT.out.html
}

// =============================================================================
// AMBIENT
//
// H5 + knee CSV -> {decontX | CellBender} -> ambient report
//
// Both methods consume the same raw H5 + knee CSV inputs from MAPPING and
// emit filtered H5 matrices under the same filename contract.
// =============================================================================

workflow AMBIENT {

    take:
    h5_ch       // tuple(sampleId, h5_file)   from MAPPING.h5_files
    knee_ch     // tuple(sampleId, knee_csv)  from MAPPING.knee_data

    main:

    def ambientMethod = (params.ambient_method ?: 'decontx').toString().trim().toLowerCase()
    if (!(ambientMethod in ['decontx', 'cellbender'])) {
        error "Unsupported --ambient_method '${ambientMethod}'. Use one of: decontx, cellbender"
    }

    joined_input = h5_ch.join(knee_ch)   // tuple(sampleId, h5, knee_csv)

    def ambient_h5_files
    def ambient_barcodes
    def ambient_report

    if (ambientMethod == 'cellbender') {
        CELLBENDER(
            joined_input,
            channel.value(file("${projectDir}/modules/ambient/cellbender_postprocess.py"))
        )

        CELLBENDER_REPORT(
            CELLBENDER.out.summaries.map { _id, csv -> csv }.collect(),
            CELLBENDER.out.labels.map { _id, csv -> csv }.collect(),
            channel.value(file("${projectDir}/modules/reports/cellbender_report.qmd"))
        )

        ambient_h5_files = CELLBENDER.out.h5_files
        ambient_barcodes = CELLBENDER.out.barcodes
        ambient_report = CELLBENDER_REPORT.out.html
    } else {
        DECONTX(
            joined_input,
            channel.value(file("${projectDir}/modules/ambient/decontx.R"))
        )

        AMBIENT_REPORT(
            DECONTX.out.qc_metrics.map { _id, csv -> csv }.collect(),
            DECONTX.out.barcodes.map   { _id, csv -> csv }.collect(),
            DECONTX.out.dcx_params.map { _id, csv -> csv }.collect(),
            DECONTX.out.summaries.map  { _id, csv -> csv }.collect(),
            channel.value(file("${projectDir}/modules/reports/ambient_report.qmd"))
        )

        ambient_h5_files = DECONTX.out.h5_files
        ambient_barcodes = DECONTX.out.barcodes
        ambient_report = AMBIENT_REPORT.out.html
    }

    STAGE_RAW_H5(h5_ch)

    AMBIENT_DE_PROC(
        STAGE_RAW_H5.out.h5.collect(),
        knee_ch.map { _id, csv -> csv }.collect(),
        ambient_h5_files.map { _id, h5 -> h5 }.collect(),
        channel.value(file(params.genome_gtf)),
        channel.value(file("${projectDir}/modules/ambient_de/ambient_de.R"))
    )

    emit:
    h5_files   = ambient_h5_files
    barcodes   = ambient_barcodes
    report     = ambient_report
    de_table   = AMBIENT_DE_PROC.out.de_table
    pb_empties = AMBIENT_DE_PROC.out.pb_empties
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
    sample_metadata = PREPARE_SAMPLE_METADATA.out.sample_metadata
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
    report = ANNOTATION_REPORT.out.html
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
    qc_metrics_ch
    integration_dt_ch
    annotation_cell_labels_ch

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
        channel.value(file("${projectDir}/modules/zoom/prepare_zoom_subset.py"))
    )

    STAGE_ZOOM_RAW_H5(raw_h5_ch)

    ZOOM_AMBIENT_DE(
        PREPARE_ZOOM_SUBSET.out.zoom_inputs,
        STAGE_ZOOM_RAW_H5.out.h5.collect(),
        knee_ch.map { _id, csv -> csv }.collect(),
        h5_ch.map { _id, h5 -> h5 }.collect(),
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
