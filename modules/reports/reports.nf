// reports.nf — report rendering processes

// ---------------------------------------------------------------------------
// Mapping report (one HTML across all samples)
// ---------------------------------------------------------------------------

process MAPPING_REPORT {
    label     "process_reports"
    tag       "mapping_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/mapping", mode: params.publish_mode_reports, overwrite: true

    input:
    path knee_data_csvs   // collected knee_plot_data_*.csv files — v1 (DropletUtils 1.30)
    path report_qmd       // mapping_report.qmd

    output:
    path "mapping_report.html", emit: html
    path "plots/**",             optional: true

    script:
    """
    export HOME="\$PWD"
    quarto render "${report_qmd}" --output mapping_report.html
    """
}

// ---------------------------------------------------------------------------
// CellBender report (one HTML across all samples)
// ---------------------------------------------------------------------------

process CELLBENDER_REPORT {
    label     "process_reports"
    tag       "cellbender_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/ambient", mode: params.publish_mode_reports, overwrite: true

    input:
    path summary_csvs   // cb_summary_*.csv
    path labels_csvs    // cb_barcode_labels_*.csv.gz
    path report_qmd     // cellbender_report.qmd

    output:
    path "cellbender_report.html", emit: html

    script:
    """
    export HOME="\$PWD"
    quarto render "${report_qmd}" --output cellbender_report.html
    """
}

// ---------------------------------------------------------------------------
// Ambient report (one HTML across all samples)
// ---------------------------------------------------------------------------

process AMBIENT_REPORT {
    label     "process_reports"
    tag       "ambient_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/ambient", mode: params.publish_mode_reports, overwrite: true

    input:
    path qc_metrics_csvs   // barcodes_qc_metrics_*.csv.gz  — pre/post S/U/A per barcode
    path barcodes_csvs     // cell_barcodes_*.csv            — accepted cell barcodes
    path dcx_params_csvs   // dcx_params_*.csv               — contamination + cluster
    path summaries_csvs    // dcx_summary_*.csv              — per-sample summary stats
    path report_qmd        // ambient_report.qmd

    output:
    path "ambient_report.html", emit: html
    path "plots/**",             optional: true

    script:
    """
    export HOME="\$PWD"
    quarto render "${report_qmd}" --output ambient_report.html
    """
}

process CELL_CALLING_REPORT {
    label     "process_reports"
    tag       "cell_calling_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/cell_calling", mode: params.publish_mode_reports, overwrite: true

    input:
    path alevinfry_stats_csvs // alevinfry_stats_*.csv         — read mapping + per-barcode stats from piscem/alevin-fry
    path summary_csvs         // cell_calling_summary_*.csv     — per-sample counts + cuts
    path label_csvs           // cell_calling_labels_*.csv.gz   — per-barcode posteriors + population
    path gmm_rds              // cell_calling_gmm_*.rds         — fitted GMM params + cuts
    path report_qmd           // cell_calling_report.qmd

    output:
    path "cell_calling_report.html", emit: html
    path "plots/**",                  optional: true

    script:
    """
    export HOME="\$PWD"
    export CC_BHATTACHARYYA_WARN="${params.cc_bhattacharyya_warn}"
    export CC_BHATTACHARYYA_FAIL="${params.cc_bhattacharyya_fail}"
    quarto render "${report_qmd}" --output cell_calling_report.html
    """
}

process KNEE_REPORT {
    label     "process_reports"
    tag       "knee_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/cell_calling", mode: params.publish_mode_reports, overwrite: true

    input:
    path knee_data_csvs      // knee_plot_data_*.csv
    path alevinfry_stats_csvs // alevinfry_stats_*.csv — read mapping + per-barcode stats from piscem/alevin-fry
    path summary_csvs        // knee_summary_*.csv
    path label_csvs          // knee_labels_*.csv.gz
    path report_qmd          // knee_report.qmd

    output:
    path "knee_report.html", emit: html
    path "plots/**",         optional: true

    script:
    """
    export HOME="\$PWD"
    quarto render "${report_qmd}" --output knee_report.html
    """
}

process EMPTYDROPS_REPORT {
    label     "process_reports"
    tag       "emptydrops_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/cell_calling", mode: params.publish_mode_reports, overwrite: true

    input:
    path summary_csvs   // emptydrops_summary_*.csv
    path label_csvs     // emptydrops_labels_*.csv.gz
    path report_qmd     // emptydrops_report.qmd

    output:
    path "emptydrops_report.html", emit: html
    path "plots/**",               optional: true

    script:
    """
    export HOME="\$PWD"
    quarto render "${report_qmd}" --output emptydrops_report.html
    """
}

process CELLSWEEP_REPORT {
    label     "process_reports"
    tag       "cellsweep_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/cellsweep", mode: params.publish_mode_reports, overwrite: true

    input:
    path label_csvs                            // <caller>_labels_*.csv.gz from AMBIENT.labels
    path qc_metrics_csvs, stageAs: 'pre_qc/*'   // qc_metrics_*.csv.gz (pre-CellSweep) from QC.qc_metrics
    path alpha_csvs                             // <sid>_alpha_hat.csv.gz from CELLSWEEP
    path annotation_cells_csvs                  // annotation_cells_<method_id>_<sid>.csv.gz from PER_SAMPLE_ANNOTATION_WF
    path corrected_qc_csvs, stageAs: 'post_qc/*' // qc_metrics_*.csv.gz (post-CellSweep, has alpha_hat) from CELLSWEEP_TO_H5
    path report_qmd                             // cellsweep_report.qmd
    path plots_r                                // integration_plots.R
    val  method_id
    val  cluster_col

    output:
    path "cellsweep_report.html", emit: html
    path "plots/**",               optional: true

    script:
    """
    export HOME="\$PWD"
    export ANNOTATION_METHOD_ID='${method_id}'
    export ANNOTATION_CLUSTER_COL='${cluster_col}'
    export CELLSWEEP_REPORT_N_WORST='${params.cellsweep_report_n_worst}'
    export CELLSWEEP_REPORT_N_BEST='${params.cellsweep_report_n_best}'
    export CELLSWEEP_REPORT_FOCUS_SAMPLES='${params.cellsweep_report_focus_samples}'
    quarto render "${report_qmd}" --output cellsweep_report.html
    """
}

// ---------------------------------------------------------------------------
// QC report (one HTML across all samples)
// ---------------------------------------------------------------------------

process QC_REPORT {
    label     "process_reports"
    tag       "qc_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/qc", mode: params.publish_mode_reports, overwrite: true

    input:
    path qc_metrics_csvs   // collected qc_metrics_*.csv.gz files
    path qc_plots_r        // qc_plots.R helper script
    path report_qmd        // qc_report.qmd

    output:
    path "qc_report.html", emit: html
    path "plots/**",        optional: true

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
    export QC_MAD_FILTER="${params.qc_mad_filter}"
    export QC_MAD_NMADS="${params.qc_mad_nmads}"
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
    publishDir "${params.outputDir}/hvg", mode: params.publish_mode_reports, overwrite: true

    input:
    path hvg_stats_csv     // hvg_stats.csv.gz — per-gene HVG stats
    path de_table_gz, stageAs: 'de_in/*'    // edger_dt.csv.gz (or NO_FILE placeholder) — ambient DE results
    path pb_empties_rds, stageAs: 'pb_in/*' // pb_empties.rds (or NO_FILE placeholder) — pseudobulk empty SummarizedExperiment
    path report_qmd        // hvg_report.qmd
    path plots_r           // hvg_plots.R — plotting helpers sourced by the report

    output:
    path "hvg_report.html", emit: html
    path "plots/**",         optional: true

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
    publishDir "${params.outputDir}/integration", mode: params.publish_mode_reports, overwrite: true

    input:
    path integration_csv   // integration_dt.csv.gz — UMAP + cluster assignments
    path dbl_sweep_csv     // dbl_sweep.csv.gz — doublet proportion sweep
    path qc_metrics_csvs   // collected qc_metrics_*.csv.gz files for cluster QC summaries
    path report_qmd        // integration_report.qmd
    path plots_r           // integration_plots.R — plotting helpers sourced by the report

    output:
    path "integration_report.html", emit: html
    path "plots/**",                 optional: true

    script:
    """
    export HOME="\$PWD"
    export METADATA_VARS="${params.metadata_vars}"
    export LEIDEN_RES="${params.integration_leiden_res}"
    export DBL_CL_PROP="${params.integration_dbl_cl_prop}"
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
    publishDir "${params.outputDir}/annotation", mode: params.publish_mode_reports, overwrite: true

    input:
    path integration_csv
    path marker_stats_csv
    path logcpms_csv
    path marker_panel_csv
    path marker_expr_rds
    path top_marker_expr_rds
    path cell_labels_csv
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

process ANNOTATION_METHOD_REPORT {
    label     "process_reports"
    tag       "annotation_method_${method_id}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/annotation", mode: params.publish_mode_reports, overwrite: true, saveAs: { filename -> "${method_id}/${filename}" }

    input:
    tuple val(method_id), val(spec_b64), path(cells_csv), path(cluster_csv), path(export_csv)
    path report_qmd

    output:
    path "annotation_report_${method_id}.html", emit: html

    script:
    def spec = new groovy.json.JsonSlurper().parseText(new String(spec_b64.decodeBase64()))
    def engine = spec.engine.toString().replace("'", "'\"'\"'")
    def referenceName = spec.reference_name.toString().replace("'", "'\"'\"'")
    def referenceLabelCol = (spec.reference_label_col ?: '').toString().replace("'", "'\"'\"'")
    """
    export HOME="\$PWD"
    export ANNOTATION_METHOD_ID='${method_id}'
    export ANNOTATION_ENGINE='${engine}'
    export ANNOTATION_REFERENCE_NAME='${referenceName}'
    export ANNOTATION_REFERENCE_LABEL_COL='${referenceLabelCol}'
    export ANNOTATION_CELLS_CSV='${cells_csv.name}'
    export ANNOTATION_CLUSTER_CSV='${cluster_csv.name}'
    export ANNOTATION_FEATURE_RES='${params.annotation_sel_res}'

    quarto render ${report_qmd} --output annotation_report_${method_id}.html
    """
}

process PER_SAMPLE_ANNOTATION_REPORT {
    label     "process_reports"
    tag       "per_sample_annotation_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/per_sample_annotation", mode: params.publish_mode_reports, overwrite: true

    input:
    path cells_csvs   // annotation_cells_<method_id>_<sid>.csv.gz, one per sample
    path summary_csvs // annotation_cluster_summary_<method_id>_<sid>.csv.gz, one per sample
    path qc_metrics_csvs // qc_metrics_<sid>.csv.gz (pre-CellSweep, all cells) from QC.qc_metrics
    path report_qmd   // per_sample_annotation_report.qmd
    path plots_r      // integration_plots.R
    path qc_plots_r   // qc_plots.R
    val  method_id
    val  cluster_col

    output:
    path "per_sample_annotation_report.html", emit: html
    path "plots/**",                          optional: true

    script:
    """
    export HOME="\$PWD"
    export ANNOTATION_METHOD_ID='${method_id}'
    export ANNOTATION_CLUSTER_COL='${cluster_col}'
    export HARD_MIN_COUNTS="${params.qc_hard_min_counts}"
    export HARD_MIN_FEATS="${params.qc_hard_min_feats}"
    export HARD_MAX_MITO="${params.qc_hard_max_mito}"
    export MIN_COUNTS="${params.qc_min_counts}"
    export MIN_FEATS="${params.qc_min_feats}"
    export MAX_MITO="${params.qc_max_mito}"
    export MIN_MITO="${params.qc_min_mito}"
    export MAX_SPLICE="${params.qc_max_splice}"
    export MIN_SPLICE="${params.qc_min_splice}"
    quarto render "${report_qmd}" --output per_sample_annotation_report.html
    """
}

process REPORT_SITE {
    label     "process_reports"
    tag       "report_site"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/reports", mode: params.publish_mode_reports, overwrite: true, saveAs: { filename -> filename.replaceFirst('^site/', '') }

    input:
    path report_htmls
    path landing_qmd
    val landing_payload_json
    path integration_csv
    path builder_script
    path trace_file

    output:
    path "site/*", emit: site

    script:
    def payloadB64 = landing_payload_json.toString().bytes.encodeBase64().toString()
    """
    set -euo pipefail
    export HOME="\$PWD"

    printf '%s' '${payloadB64}' | base64 --decode > landing_page_payload.json

    python3 - <<'PY'
import csv
import datetime as dt
import json
import re
from pathlib import Path


def parse_nf_duration_ms(s):
    if not s or s.strip() in ('', '-'):
        return None
    total = 0.0
    for val, unit in re.findall(r'([0-9.]+)(ms|s|m|h|d)', s):
        v = float(val)
        if unit == 'ms':   total += v
        elif unit == 's':  total += v * 1_000
        elif unit == 'm':  total += v * 60_000
        elif unit == 'h':  total += v * 3_600_000
        elif unit == 'd':  total += v * 86_400_000
    return total or None


def runtime_from_trace(trace_path):
    min_submit = None
    max_complete = None
    try:
        with open(trace_path, newline='', encoding='utf-8') as fh:
            for row in csv.DictReader(fh, delimiter=chr(9)):
                if row.get('status') != 'COMPLETED':
                    continue
                submit_str = (row.get('submit') or '').strip()
                dur_str    = (row.get('duration') or '').strip()
                if not submit_str or submit_str == '-':
                    continue
                submit_dt = None
                for fmt in ('%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S'):
                    try:
                        submit_dt = dt.datetime.strptime(submit_str, fmt)
                        break
                    except ValueError:
                        continue
                if submit_dt is None:
                    continue
                if min_submit is None or submit_dt < min_submit:
                    min_submit = submit_dt
                dur_ms = parse_nf_duration_ms(dur_str)
                if dur_ms is not None:
                    complete_dt = submit_dt + dt.timedelta(milliseconds=dur_ms)
                    if max_complete is None or complete_dt > max_complete:
                        max_complete = complete_dt
    except Exception:
        pass
    if min_submit is not None and max_complete is not None:
        return max(0, int((max_complete - min_submit).total_seconds()))
    return None


payload_path = Path("landing_page_payload.json")
payload = json.loads(payload_path.read_text(encoding="utf-8"))

now = dt.datetime.now().astimezone()
payload["generated_at"] = now.strftime("%Y-%m-%d %H:%M:%S %Z")

trace_path = "${trace_file}"
runtime_seconds = None
if trace_path and not trace_path.endswith("NO_FILE"):
    runtime_seconds = runtime_from_trace(trace_path)

if runtime_seconds is None:
    started_at = payload.get("started_at")
    if started_at:
        try:
            runtime_seconds = max(0, int((now - dt.datetime.fromisoformat(started_at)).total_seconds()))
        except ValueError:
            pass

payload["runtime_seconds"] = runtime_seconds
payload_path.write_text(json.dumps(payload, indent=2) + chr(10), encoding="utf-8")
PY

    quarto render "${landing_qmd}" --output index.html

    mapfile -t report_html_files < <(find -L . -type f -name '*.html' ! -path './site/*' -printf '%P\n' | sort)

    python3 "${builder_script}" \
        --payload landing_page_payload.json \
        --outdir site \
        "\${report_html_files[@]}"
    """
}
