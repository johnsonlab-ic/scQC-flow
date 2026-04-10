// workflows.nf

include { SIMPLEAF_INDEX    } from '../modules/mapping/mapping.nf'
include { SIMPLEAF_QUANT    } from '../modules/mapping/mapping.nf'
include { BARCODE_ESTIMATION } from '../modules/mapping/mapping.nf'
include { MAPPING_REPORT    } from '../modules/reports/reports.nf'
include { DECONTX           } from '../modules/ambient/ambient.nf'
include { AMBIENT_REPORT    } from '../modules/reports/reports.nf'
include { DOUBLET_DETECTION } from '../modules/qc/qc.nf'
include { APPLY_QC          } from '../modules/qc/qc.nf'
include { QC_REPORT         } from '../modules/reports/reports.nf'
include { HVG_SELECTION     } from '../modules/hvg/hvg.nf'
include { HVG_REPORT        } from '../modules/reports/reports.nf'
include { RUN_INTEGRATION   } from '../modules/integration/integration.nf'
include { INTEGRATION_REPORT } from '../modules/reports/reports.nf'
include { INDEX_REPORT      } from '../modules/reports/reports.nf'

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
// H5 + knee CSV -> decontX -> ambient report
//
// Takes the H5 and knee CSV from MAPPING (v1 barcode estimation).
// Runs decontX per sample using the in_empty_plateau flag from the knee CSV
// to define the background, and expected_cells to define the cell population.
// =============================================================================

workflow AMBIENT {

    take:
    h5_ch       // tuple(sampleId, h5_file)   from MAPPING.h5_files
    knee_ch     // tuple(sampleId, knee_csv)  from MAPPING.knee_data

    main:

    // ------------------------------------------------------------------
    // 1. Join H5 + knee CSV per sample, run decontX
    // ------------------------------------------------------------------
    dcx_input = h5_ch.join(knee_ch)   // tuple(sampleId, h5, knee_csv)

    DECONTX(
        dcx_input,
        channel.value(file("${projectDir}/modules/ambient/decontx.R"))
    )

    // ------------------------------------------------------------------
    // 2. Ambient report (one HTML across all samples)
    // ------------------------------------------------------------------
    AMBIENT_REPORT(
        DECONTX.out.qc_metrics.map { _id, csv -> csv }.collect(),
        DECONTX.out.barcodes.map   { _id, csv -> csv }.collect(),
        DECONTX.out.dcx_params.map { _id, csv -> csv }.collect(),
        DECONTX.out.summaries.map  { _id, csv -> csv }.collect(),
        channel.value(file("${projectDir}/modules/reports/ambient_report.qmd"))
    )

    emit:
    h5_files = DECONTX.out.h5_files
    barcodes = DECONTX.out.barcodes
    report   = AMBIENT_REPORT.out.html
}

// =============================================================================
// AMBIENT_DE
//
// Raw H5 + knee CSV + filtered H5 -> EdgeR comparison -> DE table + pseudobulk
//
// Takes raw H5 files (all droplets) from MAPPING, knee CSVs from MAPPING
// (to identify empty vs. cell-ranked droplets), and filtered H5 files
// (decontX-processed cells) from AMBIENT. Constructs pseudobulk matrices:
// - Empty: sum S+U+A counts across empty-ranked barcodes per sample
// - Cells: sum S+U+A counts across decontX-identified cells per sample
// Then runs EdgeR comparison to identify ambient genes (FDR<0.05, logFC>0).
//
// Outputs:
//   - edger_dt.csv.gz      differential expression results
//   - pb_empties.rds       SummarizedExperiment with empty pseudobulk counts
// =============================================================================

// =============================================================================
// QC
//
// Filtered H5 -> per-sample QC (metrics + doublet detection) -> QC report
//
// Takes decontX-filtered H5 files from AMBIENT and the genome GTF.
// Runs SAMPLE_QC per sample (scDblFinder, hard + soft thresholds),
// then renders a combined QC report.
// =============================================================================

workflow QC {

    take:
    h5_ch        // tuple(sampleId, h5_file)  from AMBIENT.h5_files

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
    de_table      // path edger_dt.csv.gz           from AMBIENT_DE.de_table
    pb_empties    // path pb_empties.rds            from AMBIENT_DE.pb_empties

    main:

    HVG_SELECTION(
        h5_ch.map { _id, h5 -> h5 }.collect(),
        qc_metrics_ch.map { _id, csv -> csv }.collect(),
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
    hvg_counts = HVG_SELECTION.out.hvg_counts
    hvg_stats  = HVG_SELECTION.out.hvg_stats
    report     = HVG_REPORT.out.html
}

// =============================================================================
// INTEGRATION
//
// HVG count matrix + metadata CSV -> Harmony integration -> integration report
//
// Loads the HVG count matrix from HVG, joins with user-provided
// metadata CSV (on metadata_id_col -> sample_id), normalises (CPM10k + log1p),
// scales, runs PCA -> Harmony -> neighbors -> Leiden (multiple resolutions)
// -> UMAP, and saves integration_dt.csv.gz.
// =============================================================================

workflow INTEGRATION {

    take:
    hvg_counts_ch  // path  hvg_counts.h5 from HVG
    qc_metrics_ch  // tuple(sampleId, csv_gz) from QC.qc_metrics

    main:

    meta_csv_ch = params.metadata_csv
        ? channel.value(file(params.metadata_csv))
        : channel.value(file("NO_FILE"))

    RUN_INTEGRATION(
        hvg_counts_ch,
        meta_csv_ch,
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
