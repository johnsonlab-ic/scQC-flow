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
// Index page linking all produced reports (dynamic)
// ---------------------------------------------------------------------------

process INDEX_REPORT {
    label     "process_reports"
    tag       "index"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path report_htmls   // collected *.html files from all workflows

    output:
    path "index.html"

    script:
    // Build nav links dynamically from the collected HTML files
    """
    set -euo pipefail
    python3 <<'PY'
import glob
import html

report_labels = {
    'mapping_report.html':    'Mapping Report',
  'barcode_report_v2.html': 'Barcode Caller V2 Report',
    'ambient_report.html':    'Ambient Report',
    'qc_report.html':         'QC Report',
    'hvg_report.html':        'HVG Report',
    'integration_report.html': 'Integration Report',
}

htmls = sorted(glob.glob('*.html'))
# Build link entries in a fixed order
ordered = [h for h in report_labels if h in htmls]
# Add any unexpected HTMLs at the end
ordered += [h for h in htmls if h not in report_labels and h != 'index.html']

links = []
for h in ordered:
    label = report_labels.get(h, h.replace('_', ' ').replace('.html', '').title())
    links.append(f'    <a href="{html.escape(h)}" target="report">{html.escape(label)}</a>')

first_src = ordered[0] if ordered else ''
nav_links = chr(10).join(links)

page = f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>scQC-flow Reports</title>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ display: flex; flex-direction: column; height: 100vh; font-family: sans-serif; }}
    nav {{
      background: #1e1e2e; color: #cdd6f4; padding: 10px 16px;
      display: flex; align-items: center; gap: 24px; flex-shrink: 0;
    }}
    nav span {{ font-weight: 600; font-size: 1rem; }}
    nav a {{
      color: #89b4fa; text-decoration: none; font-size: 0.9rem;
      padding: 4px 10px; border-radius: 4px; border: 1px solid #45475a;
    }}
    nav a:hover {{ background: #313244; }}
    iframe {{ flex: 1; border: none; width: 100%; }}
  </style>
</head>
<body>
  <nav>
    <span>scQC-flow</span>
{nav_links}
  </nav>
  <iframe name="report" src="{first_src}"></iframe>
</body>
</html>'''

with open('index.html', 'w') as f:
    f.write(page)
PY
    """
}
