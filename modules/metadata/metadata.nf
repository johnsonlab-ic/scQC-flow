process PREPARE_SAMPLE_METADATA {
    label     "process_medium"
    tag       "sample_metadata"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/metadata", mode: 'copy', overwrite: true

    input:
    val sample_ids
    path metadata_csv
    path script

    output:
    path "sample_metadata.csv.gz", emit: sample_metadata

    script:
    def sample_ids_arg = sample_ids.join(',')
    def metadata_vars_arg = params.metadata_vars ? "--metadata_vars '${params.metadata_vars}'" : ""
    """
    set -euo pipefail

    python3 ${script} \
        --metadata_csv    ${metadata_csv} \
        --metadata_id_col '${params.metadata_id_col}' \
        --sample_ids      '${sample_ids_arg}' \
        ${metadata_vars_arg} \
        --out_csv         sample_metadata.csv.gz
    """
}