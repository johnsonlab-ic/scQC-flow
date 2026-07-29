// integration.nf — Two-pass Harmony integration process

process RUN_INTEGRATION {
    label     "process_high"
    tag       "integration"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/integration", mode: params.publish_mode_nonreport, overwrite: true

    input:
    path hvg_counts      // hvg_counts.h5 (singlets) from HVG_SELECTION
    path dbl_hvg_counts  // dbl_hvg_counts.h5 (doublets) from HVG_SELECTION
    path qc_csvs         // collected qc_metrics_*.csv.gz from QC
    path script          // run_integration.py
    path dbl_qc_csvs, stageAs: 'dbl_qc/*'  // apply_qc QC (scdbl_class) for doublet rows (or NO_FILE)

    output:
    path "integration_dt.csv.gz", emit: integration_dt
    path "dbl_sweep.csv.gz",      emit: dbl_sweep

    script:
    def meta_vars    = params.metadata_vars ? "--metadata_vars '${params.metadata_vars}'" : ""
    def dbl_h5_arg   = dbl_hvg_counts.name != 'NO_FILE' ? "--dbl_hvg_h5 '${dbl_hvg_counts}'" : ""
    def has_dbl_qc   = dbl_qc_csvs instanceof List || dbl_qc_csvs.name != 'NO_FILE'
    def dbl_qc_arg   = has_dbl_qc ? "--dbl_qc_pattern 'dbl_qc/qc_metrics_*.csv.gz'" : ""
    def excl_mito    = params.exclude_mito ? "--exclude_mito" : ""
    def paga_flag    = params.integration_use_paga ? "--use_paga" : ""
    def chunk_arg    = params.integration_chunk_size > 0 ? "--chunk_size ${params.integration_chunk_size}" : ""
    def harmony_arg  = (params.harmony_impl ?: 'harmonypy') == 'harmony2' ? "--harmony_impl harmony2 --harmony2_script ${projectDir}/modules/integration/harmony2_correct.R --harmony2_lib ${params.harmony2_lib}" : ""
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
        ${paga_flag} \
        --n_dims          ${params.integration_n_dims} \
        --cluster_seed    ${params.integration_cluster_seed} \
        --dbl_res         ${params.integration_dbl_res} \
        --dbl_cl_prop     ${params.integration_dbl_cl_prop} \
        --theta           ${params.integration_theta} \
        --leiden_res      '${params.integration_leiden_res}' \
        --n_neighbors     ${params.integration_n_neighbors} \
        ${chunk_arg} \
        ${harmony_arg} \
        ${dbl_qc_arg} \
        --qc_pattern      'qc_metrics_*.csv.gz' \
        --out_csv         integration_dt.csv.gz \
        --dbl_sweep_csv   dbl_sweep.csv.gz
    """
}

// Per-sample flavor for PER_SAMPLE_ANNOTATION_WF (workflow_improvement.md): single-batch
// PCA -> Leiden(0.2) -> UMAP ahead of any cross-sample pooling. --singlet_only skips Harmony
// (nothing to correct for with n_batches=1).
process RUN_INTEGRATION_PER_SAMPLE {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/per_sample_annotation/${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(hvg_counts), path(dbl_hvg_counts), path(qc_csv)
    path script

    output:
    tuple val(sampleId), path("integration_dt_${sampleId}.csv.gz"), emit: integration_dt

    script:
    def dbl_h5_arg = dbl_hvg_counts.name != 'NO_FILE' ? "--dbl_hvg_h5 '${dbl_hvg_counts}'" : ""
    def excl_mito  = params.exclude_mito ? "--exclude_mito" : ""
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
        ${excl_mito} \
        --n_dims          ${params.integration_n_dims} \
        --cluster_seed    ${params.integration_cluster_seed} \
        --dbl_res         ${params.integration_dbl_res} \
        --dbl_cl_prop     ${params.integration_dbl_cl_prop} \
        --theta           ${params.integration_theta} \
        --leiden_res      '${params.cellsweep_celltype_col.replaceFirst(/^RNA_snn_res\./, '')}' \
        --n_neighbors     ${params.integration_n_neighbors} \
        --singlet_only \
        --qc_pattern      '${qc_csv}' \
        --out_csv         integration_dt_${sampleId}.csv.gz \
        --dbl_sweep_csv   dbl_sweep_${sampleId}.csv.gz
    """
}
