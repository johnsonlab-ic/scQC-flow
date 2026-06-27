// emptydrops.nf — EmptyDrops cell calling (CellSweep-recommended; no GMM)
//
// Drop-in alternative to CELL_CALLING: identifies cells with DropletUtils::emptyDrops
// and emits the same filt_counts / cell_barcodes / empty_barcodes contract so the
// downstream (QC / HVG / integration / AMBIENT_DE / CellSweep) is unchanged.

process EMPTYDROPS_CALLING {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/cell_calling", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5_file)
    path  script

    output:
    tuple val(sampleId), path("filt_counts_${sampleId}.h5"),             emit: h5_files
    tuple val(sampleId), path("cell_barcodes_${sampleId}.csv"),          emit: barcodes
    tuple val(sampleId), path("empty_barcodes_${sampleId}.csv"),         emit: empties
    tuple val(sampleId), path("emptydrops_labels_${sampleId}.csv.gz"),   emit: labels
    tuple val(sampleId), path("emptydrops_summary_${sampleId}.csv"),     emit: summaries
    tuple val(sampleId), path("cell_calling_${sampleId}.env"),           emit: params_env

    script:
    """
    Rscript ${script} "${sampleId}" "${h5_file}"
    """
}
