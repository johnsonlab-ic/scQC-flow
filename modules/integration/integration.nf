// integration.nf — Two-pass Harmony integration process

process RUN_INTEGRATION {
    label     "process_high"
    tag       "integration"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/integration", mode: 'copy', overwrite: true

    input:
    path hvg_counts      // hvg_counts.h5 (singlets) from HVG_SELECTION
    path dbl_hvg_counts  // dbl_hvg_counts.h5 (doublets) from HVG_SELECTION
    path qc_csvs         // collected qc_metrics_*.csv.gz from QC
    path script          // run_integration.py

    output:
    path "integration_dt.csv.gz", emit: integration_dt

    script:
    def meta_vars    = params.metadata_vars ? "--metadata_vars '${params.metadata_vars}'" : ""
    def dbl_h5_arg   = dbl_hvg_counts.name != 'NO_FILE' ? "--dbl_hvg_h5 '${dbl_hvg_counts}'" : ""
    def excl_mito    = params.exclude_mito ? "--exclude_mito" : ""
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
        --hvg_h5          ${hvg_counts} \
        ${dbl_h5_arg} \
        ${meta_vars} \
        ${excl_mito} \
        --n_dims          ${params.integration_n_dims} \
        --dbl_res         ${params.integration_dbl_res} \
        --dbl_cl_prop     ${params.integration_dbl_cl_prop} \
        --theta           ${params.integration_theta} \
        --leiden_res      '${params.integration_leiden_res}' \
        --qc_pattern      'qc_metrics_*.csv.gz' \
        --out_csv         integration_dt.csv.gz
    """
}
