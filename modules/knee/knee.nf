// knee.nf — "classic" knee/shin barcode-rank cell calling (the selection that
// originally fed decontX; decontX itself is dropped). Drop-in alternative to
// CELL_CALLING / EMPTYDROPS_CALLING with the same downstream contract.

process KNEE_CALLING {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/cell_calling", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5_file), path(knee_csv)
    path  script

    output:
    tuple val(sampleId), path("filt_counts_${sampleId}.h5"),         emit: h5_files
    tuple val(sampleId), path("cell_barcodes_${sampleId}.csv"),      emit: barcodes
    tuple val(sampleId), path("empty_barcodes_${sampleId}.csv"),     emit: empties
    tuple val(sampleId), path("knee_labels_${sampleId}.csv.gz"),     emit: labels
    tuple val(sampleId), path("knee_summary_${sampleId}.csv"),       emit: summaries
    tuple val(sampleId), path("cell_calling_${sampleId}.env"),       emit: params_env

    script:
    """
    Rscript ${script} "${sampleId}" "${h5_file}" "${knee_csv}"
    """
}
