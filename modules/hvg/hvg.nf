// hvg.nf — HVG selection process

process HVG_SELECTION {
    label     "process_high"
    tag       "hvg_selection"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/hvg", mode: 'copy', overwrite: true

    input:
    path h5_files        // collected filt_counts_*.h5 from DECONTX
    path qc_metrics_csvs // collected qc_metrics_*.csv.gz from SAMPLE_QC
    path script          // hvg_selection.py

    output:
    path "hvg_stats.csv.gz",  emit: hvg_stats
    path "hvg_counts.h5",     emit: hvg_counts

    script:
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
        --out_h5      hvg_counts.h5
    """
}
