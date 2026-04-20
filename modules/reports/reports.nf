// reports.nf — report rendering processes

// ---------------------------------------------------------------------------
// Mapping report (one HTML across all samples)
// ---------------------------------------------------------------------------

process MAPPING_REPORT {
    label     "process_reports"
    tag       "mapping_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
  publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path knee_data_csvs   // collected knee_plot_data_*.csv files — v1 (DropletUtils 1.30)
    path report_qmd       // mapping_report.qmd

    output:
    path "mapping_report.html", emit: html

    script:
    """
    export HOME="\$PWD"
    quarto render "${report_qmd}" --output mapping_report.html
    """
}

// ---------------------------------------------------------------------------
// Barcode caller v2 report (one HTML across all samples)
// ---------------------------------------------------------------------------

process BARCODE_REPORT_V2 {
    label     "process_reports"
    tag       "barcode_report_v2"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
  publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path audit_csvs     // barcode_audit_v2_*.csv.gz
    path summary_csvs   // barcode_summary_v2_*.csv
    path report_qmd     // barcode_report_v2.qmd
    path plots_r        // barcode_v2_plots.R

    output:
    path "barcode_report_v2.html", emit: html

    script:
    """
    export HOME="\$PWD"
    export BARCODE_V2_SPLICE_CONTEXT="${params.barcode_v2_splice_context}"
    export BARCODE_V2_ED_FDR="${params.barcode_v2_ed_fdr}"
    quarto render "${report_qmd}" --output barcode_report_v2.html
    """
}

// ---------------------------------------------------------------------------
// Ambient report (one HTML across all samples)
// ---------------------------------------------------------------------------

process AMBIENT_REPORT {
    label     "process_reports"
    tag       "ambient_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
  publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path qc_metrics_csvs   // barcodes_qc_metrics_*.csv.gz  — pre/post S/U/A per barcode
    path barcodes_csvs     // cell_barcodes_*.csv            — accepted cell barcodes
    path dcx_params_csvs   // dcx_params_*.csv               — contamination + cluster
    path summaries_csvs    // dcx_summary_*.csv              — per-sample summary stats
    path report_qmd        // ambient_report.qmd

    output:
    path "ambient_report.html", emit: html

    script:
    """
    export HOME="\$PWD"
    quarto render "${report_qmd}" --output ambient_report.html
    """
}

// ---------------------------------------------------------------------------
// QC report (one HTML across all samples)
// ---------------------------------------------------------------------------

process QC_REPORT {
    label     "process_reports"
    tag       "qc_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
  publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path qc_metrics_csvs   // collected qc_metrics_*.csv.gz files
    path qc_plots_r        // qc_plots.R helper script
    path report_qmd        // qc_report.qmd

    output:
    path "qc_report.html", emit: html

    script:
    """
    export HOME="\$PWD"
    export HARD_MIN_COUNTS="${params.qc_hard_min_counts}"
    export HARD_MIN_FEATS="${params.qc_hard_min_feats}"
    export HARD_MAX_MITO="${params.qc_hard_max_mito}"
    export MIN_COUNTS="${params.qc_min_counts}"
    export MIN_FEATS="${params.qc_min_feats}"
    export MAX_MITO="${params.qc_max_mito}"
    export MIN_MITO="${params.qc_min_mito}"
    export MAX_SPLICE="${params.qc_max_splice}"
    export MIN_SPLICE="${params.qc_min_splice}"
    quarto render "${report_qmd}" --output qc_report.html
    """
}

// ---------------------------------------------------------------------------
// HVG report (one HTML across all samples)
// ---------------------------------------------------------------------------

process HVG_REPORT {
    label     "process_reports"
    tag       "hvg_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
  publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path hvg_stats_csv     // hvg_stats.csv.gz — per-gene HVG stats
    path de_table_gz       // edger_dt.csv.gz — ambient DE results
    path pb_empties_rds    // pb_empties.rds — pseudobulk empty SummarizedExperiment
    path report_qmd        // hvg_report.qmd
    path plots_r           // hvg_plots.R — plotting helpers sourced by the report

    output:
    path "hvg_report.html", emit: html

    script:
    """
    export HOME="\$PWD"
    export N_HVG="${params.hvg_n_hvgs}"
    quarto render "${report_qmd}" --output hvg_report.html
    """
}

// ---------------------------------------------------------------------------
// Integration report (one HTML across all samples)
// ---------------------------------------------------------------------------

process INTEGRATION_REPORT {
    label     "process_reports"
    tag       "integration_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
  publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path integration_csv   // integration_dt.csv.gz — UMAP + cluster assignments
    path qc_metrics_csvs   // collected qc_metrics_*.csv.gz files for cluster QC summaries
    path report_qmd        // integration_report.qmd
    path plots_r           // integration_plots.R — plotting helpers sourced by the report

    output:
    path "integration_report.html", emit: html

    script:
    """
    export HOME="\$PWD"
    export METADATA_VARS="${params.metadata_vars}"
    export LEIDEN_RES="${params.integration_leiden_res}"
    quarto render "${report_qmd}" --output integration_report.html
    """
}

  // ---------------------------------------------------------------------------
  // Annotation report (one HTML across all samples)
  // ---------------------------------------------------------------------------

process ANNOTATION_REPORT {
    label     "process_reports"
    tag       "annotation_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path integration_csv
    path marker_stats_csv
    path logcpms_csv
    path marker_panel_csv
    path marker_expr_rds
    path report_qmd
    path utils_r
    path plots_r

    output:
    path "annotation_report.html", emit: html

    script:
    """
    export HOME="\$PWD"
    export ANNOTATION_SEL_RES="${params.annotation_sel_res}"
    export ANNOTATION_MIN_CPM_MKR="${params.annotation_min_cpm_mkr}"
    export ANNOTATION_NOT_OK_RE="${params.annotation_not_ok_re}"
    export ANNOTATION_TOP_N="${params.annotation_top_n}"
    export ANNOTATION_FDR_CUT="${params.annotation_fdr_cut}"
    export ANNOTATION_MAX_ZERO_P="${params.annotation_max_zero_p}"
    quarto render "${report_qmd}" --output annotation_report.html
    """
}

process REPORT_SITE {
    label     "process_reports"
    tag       "report_site"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true, saveAs: { filename -> filename.replaceFirst('^site/', '') }

    input:
    path report_htmls
    path landing_qmd
    val landing_payload_json
    path builder_script
    path site_css
    path site_js

    output:
    path "site/*", emit: site

    script:
    def payloadB64 = landing_payload_json.toString().bytes.encodeBase64().toString()
    """
    set -euo pipefail
    export HOME="\$PWD"

    printf '%s' '${payloadB64}' | base64 --decode > landing_page_payload.json

    quarto render "${landing_qmd}" --output index.html

    python3 "${builder_script}" \
        --payload landing_page_payload.json \
        --outdir site \
        --css "${site_css}" \
        --js "${site_js}" \
        *.html
    """
}
