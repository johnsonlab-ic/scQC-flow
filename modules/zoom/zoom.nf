// zoom.nf — post-integration zoom workflows

process STAGE_ZOOM_RAW_H5 {
    label     "process_low"
    tag       "$sampleId"

    input:
    tuple val(sampleId), path(h5_file)

    output:
    path "barcode_matrix_${sampleId}.h5", emit: h5

    script:
    """
    ln -s ${h5_file} barcode_matrix_${sampleId}.h5
    """
}

process PREPARE_ZOOM_SUBSET {
    label     "process_low"
    tag       "prepare_zoom_${zoom_name}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    publishDir "${params.outputDir}/zoom", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${zoom_name}/${filename}" }

    input:
    tuple val(zoom_name), val(spec_b64)
    path integration_csv
    path qc_metrics_csvs
    path annotation_cell_labels_csv
    path script

    output:
    tuple val(zoom_name), val(spec_b64), path("zoom_qc_metrics.csv.gz"), path("zoom_selection.csv.gz"), emit: zoom_inputs

    script:
    def annotation_arg = annotation_cell_labels_csv.name == 'NO_FILE' ? 'NO_FILE' : annotation_cell_labels_csv.name
    """
    set -euo pipefail
    export HOME="\$PWD"

    python3 ${script} \
        --integration_csv ${integration_csv} \
        --annotation_labels_csv ${annotation_arg} \
        --qc_pattern 'qc_metrics_*.csv.gz' \
        --zoom_spec_b64 ${spec_b64} \
        --out_qc_csv zoom_qc_metrics.csv.gz \
        --out_selection_csv zoom_selection.csv.gz
    """
}


process ZOOM_AMBIENT_DE {
    label     "process_high"
    tag       "zoom_ambient_${zoom_name}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    publishDir "${params.outputDir}/zoom", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${zoom_name}/${filename}" }

    input:
    tuple val(zoom_name), val(spec_b64), path(zoom_qc_metrics), path(zoom_selection)
    path raw_h5_files
    path knee_csvs
    path filt_h5_files
    path genome_gtf
    path script

    output:
    tuple val(zoom_name), val(spec_b64), path("zoom_qc_metrics.csv.gz"), path("zoom_selection.csv.gz"), path("zoom_edger_dt.csv.gz"), emit: zoom_ambient

    script:
    """
    set -euo pipefail

    Rscript ${script} ${genome_gtf} zoom_qc_metrics.csv.gz
    mv edger_dt.csv.gz zoom_edger_dt.csv.gz
    """
}


process ZOOM_HVG_SELECTION {
    label     "process_high"
    tag       "zoom_hvg_${zoom_name}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    publishDir "${params.outputDir}/zoom", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${zoom_name}/${filename}" }

    input:
    tuple val(zoom_name), val(spec_b64), path(zoom_qc_metrics), path(zoom_selection), path(zoom_edger_csv)
    path h5_files
    path genome_gtf
    path script

    output:
    tuple val(zoom_name), val(spec_b64), path("zoom_qc_metrics.csv.gz"), path("zoom_selection.csv.gz"), path("zoom_hvg_stats.csv.gz"), path("zoom_hvg_counts.h5"), path("zoom_dbl_hvg_counts.h5"), emit: zoom_hvg

    script:
    def spec = new groovy.json.JsonSlurper().parseText(new String(spec_b64.decodeBase64()))
    def nTopGenes = (spec.hvg_n_hvgs ?: params.hvg_n_hvgs) as Integer
    """
    set -euo pipefail
    export HOME="\$PWD"
    export MPLCONFIGDIR="\$PWD"
    export NUMBA_CACHE_DIR="\$PWD"

    python3 ${script} \
        --h5_pattern  'filt_counts_*.h5' \
        --qc_pattern  'zoom_qc_metrics.csv.gz' \
        --n_top_genes ${nTopGenes} \
        --out_stats   zoom_hvg_stats.csv.gz \
        --out_h5      zoom_hvg_counts.h5 \
        --out_dbl_h5  zoom_dbl_hvg_counts.h5 \
        --gtf         ${genome_gtf} \
        --edger_csv   ${zoom_edger_csv}
    """
}


process RUN_ZOOM_INTEGRATION {
    label     "process_high"
    tag       "zoom_integration_${zoom_name}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    publishDir "${params.outputDir}/zoom", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${zoom_name}/${filename}" }

    input:
    tuple val(zoom_name), val(spec_b64), path(zoom_qc_metrics), path(zoom_selection), path(zoom_hvg_stats), path(zoom_hvg_counts), path(zoom_dbl_hvg_counts)
    path script

    output:
    tuple val(zoom_name), val(spec_b64), path("zoom_qc_metrics.csv.gz"), path("zoom_selection.csv.gz"), path("zoom_integration_dt.csv.gz"), emit: zoom_int

    script:
    def spec = new groovy.json.JsonSlurper().parseText(new String(spec_b64.decodeBase64()))
    def metadataVars = (spec.metadata_vars ?: params.metadata_vars ?: '').toString()
    def metaArg = metadataVars ? "--metadata_vars '${metadataVars.replace("'", "'\"'\"'")}'" : ''
    def excludeMito = (spec.containsKey('exclude_mito') ? spec.exclude_mito : params.exclude_mito) ? '--exclude_mito' : ''
    def nDims = (spec.integration_n_dims ?: params.integration_n_dims) as Integer
    def dblRes = (spec.integration_dbl_res ?: params.integration_dbl_res) as BigDecimal
    def dblClProp = (spec.integration_dbl_cl_prop ?: params.integration_dbl_cl_prop) as BigDecimal
    def theta = (spec.integration_theta ?: params.integration_theta) as BigDecimal
    def leidenRes = (spec.integration_leiden_res ?: params.integration_leiden_res).toString().replace("'", "'\"'\"'")
    """
    set -euo pipefail
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export NUMBA_CACHE_DIR="\$PWD/.numba"
    export PYTHONUNBUFFERED=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    export NUMEXPR_MAX_THREADS=${task.cpus}

    python3 -u ${script} \
        --hvg_h5          ${zoom_hvg_counts} \
        --singlet_only \
        ${metaArg} \
        ${excludeMito} \
        --n_dims          ${nDims} \
        --dbl_res         ${dblRes} \
        --dbl_cl_prop     ${dblClProp} \
        --theta           ${theta} \
        --leiden_res      '${leidenRes}' \
        --qc_pattern      'zoom_qc_metrics.csv.gz' \
        --out_csv         zoom_integration_dt.csv.gz
    """
}


process RUN_ZOOM_MARKERS {
    label     "process_high"
    tag       "zoom_markers_${zoom_name}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    publishDir "${params.outputDir}/zoom", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${zoom_name}/${filename}" }

    input:
    tuple val(zoom_name), val(spec_b64), path(zoom_qc_metrics), path(zoom_selection), path(zoom_integration_dt)
    path h5_files
    path genome_gtf
    path script
    path utils_r

    output:
    tuple val(zoom_name), val(spec_b64), path("zoom_qc_metrics.csv.gz"), path("zoom_selection.csv.gz"), path("zoom_integration_dt.csv.gz"), path("zoom_marker_stats.csv.gz"), path("zoom_marker_logcpms.csv.gz"), emit: zoom_markers

    script:
    def spec = new groovy.json.JsonSlurper().parseText(new String(spec_b64.decodeBase64()))
    def markerSelRes = (spec.marker_sel_res ?: '0.2').toString().replace("'", "'\"'\"'")
    def markerMinClSize = (spec.marker_min_cl_size ?: 100) as Integer
    def markerMinCells = (spec.marker_min_cells ?: 10) as Integer
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        zoom_integration_dt.csv.gz \
        ${genome_gtf} \
        '${markerSelRes}' \
        ${markerMinClSize} \
        ${markerMinCells} \
        zoom_marker_stats.csv.gz \
        zoom_marker_logcpms.csv.gz \
        'filt_counts_*.h5'
    """
}


process ZOOM_REPORT {
    label     "process_reports"
    tag       "zoom_report_${zoom_name}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    publishDir "${params.outputDir}/zoom", mode: params.publish_mode_reports, overwrite: true

    input:
    tuple val(zoom_name), val(spec_b64), path(zoom_qc_metrics), path(zoom_selection), path(zoom_integration_dt), path(zoom_marker_stats), path(zoom_marker_logcpms)
    path report_qmd
    path integration_plots_r
    path annotation_utils_r

    output:
    path "zoom_${zoom_name}_report.html", emit: html

    script:
    def spec = new groovy.json.JsonSlurper().parseText(new String(spec_b64.decodeBase64()))
    def zoomSource = spec.source.toString().replace("'", "'\"'\"'")
    def zoomValues = ((spec.values ?: (spec.label ? [spec.label] : [])) as List).collect { value -> value.toString() }.join(',').replace("'", "'\"'\"'")
    def markerSelRes = (spec.marker_sel_res ?: '0.2').toString().replace("'", "'\"'\"'")
    def markerTopN = (spec.marker_top_n ?: 10) as Integer
    def markerMinCpm = (spec.marker_min_cpm ?: params.annotation_min_cpm_mkr) as Integer
    def markerFdrCut = (spec.marker_fdr_cut ?: params.annotation_fdr_cut) as BigDecimal
    def markerMaxZero = (spec.marker_max_zero_p ?: params.annotation_max_zero_p) as BigDecimal
    def metadataVars = (spec.metadata_vars ?: params.metadata_vars ?: '').toString().replace("'", "'\"'\"'")
    """
    export HOME="\$PWD"
    export ZOOM_NAME='${zoom_name}'
    export ZOOM_SOURCE='${zoomSource}'
    export ZOOM_VALUES='${zoomValues}'
    export ZOOM_MARKER_SEL_RES='${markerSelRes}'
    export ZOOM_MARKER_TOP_N='${markerTopN}'
    export ZOOM_MARKER_MIN_CPM='${markerMinCpm}'
    export ZOOM_MARKER_FDR_CUT='${markerFdrCut}'
    export ZOOM_MARKER_MAX_ZERO_P='${markerMaxZero}'
    export METADATA_VARS='${metadataVars}'

    quarto render ${report_qmd} --output zoom_${zoom_name}_report.html
    """
}