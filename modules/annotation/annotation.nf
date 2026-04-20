process RUN_ANNOTATION_MARKERS {
    label     "process_high"
    tag       "annotation"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/annotation", mode: 'copy', overwrite: true

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
        'filt_counts_*.h5'
    """
}