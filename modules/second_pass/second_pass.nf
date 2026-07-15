// second_pass.nf — adapter turning CellSweep-corrected h5ads back into the
// per-sample filt_counts/qc_metrics contract the downstream chain consumes.

process CELLSWEEP_TO_H5 {
    label     "process_medium"
    tag       "$sampleId"
    cache     'lenient'   // inputs are pass-1 published symlinks; ignore timestamps so -resume hits
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/corrected_counts", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5ad), path(pass1_qc, stageAs: 'pass1_qc_in/*')
    path  script

    output:
    tuple val(sampleId), path("filt_counts_${sampleId}.h5"),       emit: h5_files
    tuple val(sampleId), path("qc_metrics_${sampleId}.csv.gz"),    emit: qc_metrics

    script:
    """
    export HOME="\$PWD"
    export NUMBA_CACHE_DIR="\$PWD/.numba"
    export MPLCONFIGDIR="\$PWD/.mpl"
    python3 ${script} \\
      --sample_id ${sampleId} \\
      --h5ad ${h5ad} \\
      --pass1_qc ${pass1_qc} \\
      --alpha_max ${params.cell_calling_alpha_max} \\
      --out_h5 filt_counts_${sampleId}.h5 \\
      --out_qc qc_metrics_${sampleId}.csv.gz
    """
}
