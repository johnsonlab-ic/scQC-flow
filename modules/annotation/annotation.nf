process RUN_ANNOTATION_MARKERS {
    label     "process_high"
    tag       "annotation"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/annotation", mode: params.publish_mode_nonreport, overwrite: true

    input:
    path h5_files
    path integration_csv
    path genome_gtf
    path marker_csv
    path script
    path utils_r

    output:
    path "annotation_markers.csv.gz", emit: markers
    path "annotation_logcpms.csv.gz", emit: logcpms
    path "annotation_pseudobulk.rds", emit: pseudobulk
    path "annotation_marker_panel.csv.gz", emit: marker_panel
    path "annotation_marker_expression.rds", emit: marker_expr
    path "annotation_cell_labels.csv.gz", emit: cell_labels

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        ${integration_csv} \
        ${genome_gtf} \
        ${marker_csv} \
        '${params.annotation_sel_res}' \
        ${params.annotation_min_cl_size} \
        ${params.annotation_min_cells} \
        annotation_markers.csv.gz \
        annotation_logcpms.csv.gz \
        annotation_marker_panel.csv.gz \
        annotation_marker_expression.rds \
        annotation_cell_labels.csv.gz \
        annotation_pseudobulk.rds \
        ${task.cpus} \
        'filt_counts_*.h5'
    """
}

process PREP_REPORT_INPUTS {
    label     "process_high"
    tag       "prep_report_inputs"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/annotation", mode: params.publish_mode_nonreport, overwrite: true

    input:
    path h5_files
    path integration_csv
    path marker_stats_csv
    path script
    path utils_r

    output:
    path "annotation_top_marker_expression.rds", emit: top_marker_expr

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"
    export ANNOTATION_SEL_RES="${params.annotation_sel_res}"
    export ANNOTATION_MIN_CPM_MKR="${params.annotation_min_cpm_mkr}"
    export ANNOTATION_NOT_OK_RE="${params.annotation_not_ok_re}"
    export ANNOTATION_TOP_N="${params.annotation_top_n}"
    export ANNOTATION_FDR_CUT="${params.annotation_fdr_cut}"
    export ANNOTATION_MAX_ZERO_P="${params.annotation_max_zero_p}"

    Rscript ${script} \
        ${integration_csv} \
        ${marker_stats_csv} \
        annotation_top_marker_expression.rds \
        'filt_counts_*.h5'
    """
}

process FORMAT_ANNOTATION_EXPORT_METADATA {
    label     "process_low"
    tag       "annotation_export_metadata"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/annotation", mode: params.publish_mode_nonreport, overwrite: true

    input:
    path cell_labels_csv
    path script

    output:
    path "annotation_export_marker_panel.csv.gz", emit: export_metadata

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --input_csv ${cell_labels_csv} \
        --prefix marker_panel \
        --output_csv annotation_export_marker_panel.csv.gz
    """
}

process PREPARE_ANNOTATION_QUERY {
    label     "process_medium"
    tag       "prepare_annotation_query_${sample_id}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"

    input:
    tuple val(sample_id), path(h5_file), path(integration_csv)
    path script
    path utils_r

    output:
    tuple val(sample_id), path("annotation_query_${sample_id}.rds"), emit: query

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --sample_id '${sample_id}' \
        --integration_csv ${integration_csv} \
        --utils_r ${utils_r} \
        --h5_file ${h5_file} \
        --output_rds annotation_query_${sample_id}.rds
    """
}

process RUN_SINGLER_REFERENCE_ANNOTATION {
    label     "process_high"
    tag       "${method_id}_${sample_id}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/annotation", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${method_id}/${filename}" }

    input:
    tuple val(method_id), val(spec_b64), val(sample_id), path(query_rds), path(reference_rds)
    path script

    output:
    tuple val(method_id), val(spec_b64), val(sample_id), path("annotation_cells_${method_id}_${sample_id}.csv.gz"), path("annotation_cluster_summary_${method_id}_${sample_id}.csv.gz"), path("annotation_export_${method_id}_${sample_id}.csv.gz"), emit: result

    script:
    def spec = new groovy.json.JsonSlurper().parseText(new String(spec_b64.decodeBase64()))
    def labelCol = spec.reference_label_col.toString().replace("'", "'\"'\"'")
    def referenceName = spec.reference_name.toString().replace("'", "'\"'\"'")
    def fineTune = spec.containsKey('fine_tune') ? spec.fine_tune.toString() : 'false'
    def prune = spec.containsKey('prune') ? spec.prune.toString() : 'true'
    def bpType = (spec.bp_type ?: 'multicore').toString().replace("'", "'\"'\"'")
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --query_rds ${query_rds} \
        --reference_rds ${reference_rds} \
        --reference_label_col '${labelCol}' \
        --method_id '${method_id}' \
        --reference_name '${referenceName}' \
        --output_cells_csv annotation_cells_${method_id}_${sample_id}.csv.gz \
        --output_export_csv annotation_export_${method_id}_${sample_id}.csv.gz \
        --output_cluster_csv annotation_cluster_summary_${method_id}_${sample_id}.csv.gz \
        --ncores ${task.cpus} \
        --bp_type '${bpType}' \
        --fine_tune ${fineTune} \
        --prune ${prune}
    """
}

process RUN_XGBOOST_REFERENCE_ANNOTATION {
    label     "process_high"
    tag       "${method_id}_${sample_id}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/annotation", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${method_id}/${filename}" }

    input:
    tuple val(method_id), val(spec_b64), val(sample_id), path(query_rds), path(model_rds), path(class_csv)
    path script

    output:
    tuple val(method_id), val(spec_b64), val(sample_id), path("annotation_cells_${method_id}_${sample_id}.csv.gz"), path("annotation_cluster_summary_${method_id}_${sample_id}.csv.gz"), path("annotation_export_${method_id}_${sample_id}.csv.gz"), emit: result

    script:
    def spec = new groovy.json.JsonSlurper().parseText(new String(spec_b64.decodeBase64()))
    def clusterArg = spec.cluster_col ? "--cluster_col '${spec.cluster_col.toString().replace("'", "'\"'\"'")}'" : ''
    def referenceName = spec.reference_name.toString().replace("'", "'\"'\"'")
    def chunkSize = (spec.chunk_size ?: 10000) as Integer
    def scaleFactor = (spec.scale_factor ?: 10000) as BigDecimal
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --query_rds ${query_rds} \
        --model_rds ${model_rds} \
        --class_csv ${class_csv} \
        ${clusterArg} \
        --method_id '${method_id}' \
        --reference_name '${referenceName}' \
        --output_cells_csv annotation_cells_${method_id}_${sample_id}.csv.gz \
        --output_export_csv annotation_export_${method_id}_${sample_id}.csv.gz \
        --output_cluster_csv annotation_cluster_summary_${method_id}_${sample_id}.csv.gz \
        --chunk_size ${chunkSize} \
        --scale_factor ${scaleFactor}
    """
}

process COMBINE_ANNOTATION_METHOD_OUTPUTS {
    label     "process_medium"
    tag       "combine_${method_id}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/annotation", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${method_id}/${filename}" }

    input:
    tuple val(method_id), val(spec_b64), path(cells_csvs), path(export_csvs)
    path script

    output:
    tuple val(method_id), val(spec_b64), path("annotation_cells_${method_id}.csv.gz"), path("annotation_cluster_summary_${method_id}.csv.gz"), path("annotation_export_${method_id}.csv.gz"), emit: result

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --method_id '${method_id}' \
        --cells_pattern 'annotation_cells_*.csv.gz' \
        --export_pattern 'annotation_export_*.csv.gz' \
        --output_cells_csv annotation_cells_${method_id}.csv.gz \
        --output_cluster_csv annotation_cluster_summary_${method_id}.csv.gz \
        --output_export_csv annotation_export_${method_id}.csv.gz
    """
}