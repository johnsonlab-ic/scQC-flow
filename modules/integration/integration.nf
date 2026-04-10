// integration.nf — Harmony integration process

process RUN_INTEGRATION {
    label     "process_high"
    tag       "integration"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/integration", mode: 'copy', overwrite: true

    input:
    path hvg_counts    // hvg_counts.h5 from HVG_SELECTION
    path metadata_csv  // user-supplied metadata CSV
    path script        // run_integration.py

    output:
    path "integration_dt.csv.gz", emit: integration_dt

    script:
    def meta_vars = params.metadata_vars ? "--metadata_vars '${params.metadata_vars}'" : ""
    def meta_csv  = metadata_csv.name != 'NO_FILE' ? "--metadata_csv '${metadata_csv}'" : ""
    """
    set -euo pipefail
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export NUMBA_CACHE_DIR="\$PWD/.numba"
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    python3 ${script} \
        --hvg_h5          ${hvg_counts} \
        ${meta_csv} \
        --metadata_id_col '${params.metadata_id_col}' \
        ${meta_vars} \
        --n_dims          ${params.integration_n_dims} \
        --theta           ${params.integration_theta} \
        --leiden_res      '${params.integration_leiden_res}' \
        --out_csv         integration_dt.csv.gz
    """
}
