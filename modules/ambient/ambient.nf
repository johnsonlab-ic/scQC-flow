// ambient.nf — decontX ambient RNA removal

process DECONTX {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/ambient/dcx_${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5_file), path(knee_csv)
    path  script

    output:
    tuple val(sampleId), path("filt_counts_${sampleId}.h5"),                emit: h5_files
    tuple val(sampleId), path("cell_barcodes_${sampleId}.csv"),              emit: barcodes
    tuple val(sampleId), path("dcx_params_${sampleId}.csv"),                 emit: dcx_params
    tuple val(sampleId), path("barcodes_qc_metrics_${sampleId}.csv.gz"),     emit: qc_metrics
    tuple val(sampleId), path("dcx_summary_${sampleId}.csv"),                emit: summaries
    tuple val(sampleId),
          env('DCX_EXPECTED_CELLS'),
          env('DCX_TOTAL_INCLUDED'),
          env('DCX_MEAN_CONTAMINATION'),
          emit: ambient_params

    script:
    """
    set -euo pipefail

    Rscript ${script} "${sampleId}" "${h5_file}" "${knee_csv}"

    source "ambient_params_${sampleId}.env"
    export DCX_EXPECTED_CELLS
    export DCX_TOTAL_INCLUDED
    export DCX_MEAN_CONTAMINATION
    """
}
