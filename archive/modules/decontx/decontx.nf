// DECONTX
//
// CPU-native ambient RNA removal using decontX (barcodeRanks cell calling).
// GPU-free alternative to CellBender.
//
// Cell calling mirrors scprocess/scripts/ambient.R:
//   - Cells  : top expected_cells barcodes from knee_data (sorted by rank)
//   - Empties: barcodes with in_empty_plateau == TRUE in knee_data
//   No emptyDrops is run — it is too slow on unfiltered matrices.
//
// Inputs:
//   h5_file       — unfiltered total (S+U+A) count matrix from BARCODE_ESTIMATION
//   knee_data_csv — per-barcode CSV from BARCODE_ESTIMATION
//
// Outputs (match CellBender interface for unified downstream use):
//   filtered_h5    — ambient-corrected, cell-only count matrix (.h5)
//   barcodes       — retained barcode list CSV (no header)
//   contamination  — per-barcode contamination fraction + decontX cluster
//   usa_metrics    — per-barcode pre/post S+U+A counts for top total_droplets_included
//                    barcodes; input for DECONTX_REPORT

process DECONTX {
    label "process_high"
    tag  "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/decontx", mode: 'copy', overwrite: true

    input:
    tuple val(sampleId),
          path(h5_file),
          path(knee_data_csv)
    path(script)

    output:
    tuple val(sampleId), path("${sampleId}_decontx_filtered.h5"),       emit: filtered_h5
    tuple val(sampleId), path("${sampleId}_decontx_barcodes.csv"),      emit: barcodes
    tuple val(sampleId), path("${sampleId}_decontx_contamination.csv"), emit: contamination
    tuple val(sampleId), path("${sampleId}_decontx_usa_metrics.csv"),   emit: usa_metrics

    script:
    """
    Rscript "${script}" \\
        "${sampleId}" \\
        "${h5_file}" \\
        "${knee_data_csv}" \\
        "${task.cpus}" \\
        "${sampleId}_decontx_filtered.h5" \\
        "${sampleId}_decontx_barcodes.csv" \\
        "${sampleId}_decontx_contamination.csv" \\
        "${sampleId}_decontx_usa_metrics.csv"
    """
}

// =============================================================================
// DECONTX_REPORT
//
// Collects per-sample decontX outputs and renders a single ambient QC report.
// Mirrors scprocess ambient.Rmd — three sections:
//   1. Spliced % vs UMIs (tabset per sample, pre/post decontX)
//   2. Reads removed as ambient (violin per sample, chunked)
//   3. Droplets kept table (one row per sample)
// =============================================================================
process DECONTX_REPORT {
    label "process_reports"
    tag  "decontx_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path(usa_metrics_csvs)   // collected list of *_decontx_usa_metrics.csv
    path(barcodes_csvs)      // collected list of *_decontx_barcodes.csv
    path(knee_data_csvs)     // collected list of *_knee_data.csv (for total_droplets_included)
    path(report_qmd)

    output:
    path("decontx_report.html"), emit: html_report

    script:
    """
    export HOME="\$PWD"
    quarto render "${report_qmd}" --output decontx_report.html
    """
}
