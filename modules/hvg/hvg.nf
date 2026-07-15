// hvg.nf — HVG selection process

process HVG_SELECTION {
    label     "process_high"
    tag       "hvg_selection"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/hvg", mode: params.publish_mode_nonreport, overwrite: true

    input:
    path h5_files        // collected filt_counts_*.h5 from DECONTX
    path qc_metrics_csvs // collected qc_metrics_*.csv.gz from APPLY_QC
    path genome_gtf      // reference GTF (for SYMBOL_ENSEMBL gene IDs)
    path edger_csv       // edger_dt.csv.gz from AMBIENT_DE (ambient gene exclusion)
    path script          // hvg_selection.py

    output:
    path "hvg_stats.csv.gz",    emit: hvg_stats
    path "hvg_counts.h5",       emit: hvg_counts
    path "dbl_hvg_counts.h5",   emit: dbl_hvg_counts

    script:
    def edger_arg = edger_csv.name != 'NO_FILE' ? "--edger_csv '${edger_csv}'" : ""
    """
    set -euo pipefail
    export HOME="\$PWD"
    export MPLCONFIGDIR="\$PWD"
    export NUMBA_CACHE_DIR="\$PWD"

    python3 ${script} \
        --h5_pattern  'filt_counts_*.h5' \
        --qc_pattern  'qc_metrics_*.csv.gz' \
        --n_top_genes ${params.hvg_n_hvgs} \
        --out_stats   hvg_stats.csv.gz \
        --out_h5      hvg_counts.h5 \
        --out_dbl_h5  dbl_hvg_counts.h5 \
        --gtf         ${genome_gtf} \
        ${edger_arg}
    """
}

// Per-sample flavor for PER_SAMPLE_ANNOTATION_WF (workflow_improvement.md): same script,
// but tagged/scoped to one sample instead of the whole collected cohort, so it can run as
// an independent parallel task per sample ahead of any cross-sample pooling.
process HVG_SELECTION_PER_SAMPLE {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/per_sample_annotation/${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5_file), path(qc_metrics_csv)
    path genome_gtf
    path edger_csv
    path script

    output:
    tuple val(sampleId), path("hvg_stats_${sampleId}.csv.gz"),  emit: hvg_stats
    tuple val(sampleId), path("hvg_counts_${sampleId}.h5"),     emit: hvg_counts
    tuple val(sampleId), path("dbl_hvg_counts_${sampleId}.h5"), emit: dbl_hvg_counts

    script:
    def edger_arg = edger_csv.name != 'NO_FILE' ? "--edger_csv '${edger_csv}'" : ""
    """
    set -euo pipefail
    export HOME="\$PWD"
    export MPLCONFIGDIR="\$PWD"
    export NUMBA_CACHE_DIR="\$PWD"

    python3 ${script} \
        --h5_pattern  '${h5_file}' \
        --qc_pattern  '${qc_metrics_csv}' \
        --n_top_genes ${params.hvg_n_hvgs} \
        --out_stats   hvg_stats_${sampleId}.csv.gz \
        --out_h5      hvg_counts_${sampleId}.h5 \
        --out_dbl_h5  dbl_hvg_counts_${sampleId}.h5 \
        --gtf         ${genome_gtf} \
        ${edger_arg}
    """
}
