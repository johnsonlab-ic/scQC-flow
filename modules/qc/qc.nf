// qc.nf — cell-level quality control processes

// ---------------------------------------------------------------------------
// DOUBLET_DETECTION: per-sample scDblFinder (slow — cached by Nextflow)
//
// Takes the decontX-filtered H5 and GTF; excludes rRNA genes and runs
// scDblFinder on the combined S+U+A counts matrix.  The GTF is used only
// to strip high-expression rRNA genes before doublet scoring.
//
// This process is intentionally separated from APPLY_QC so that changing
// QC thresholds does not re-trigger doublet detection.
// ---------------------------------------------------------------------------

process DOUBLET_DETECTION {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/qc/doublet_detection/${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5_file)
    path  genome_gtf
    path  script

    output:
    tuple val(sampleId), path("scdblfinder_${sampleId}.csv.gz"), emit: dbl_results

    script:
    """
    set -euo pipefail

    Rscript ${script} \
        "${sampleId}" \
        "${h5_file}" \
        "${genome_gtf}"
    """
}

// ---------------------------------------------------------------------------
// APPLY_QC: per-sample metric calculation + threshold filtering (fast)
//
// Takes the decontX H5, GTF, and pre-computed doublet results, then
// calculates per-cell QC metrics and applies hard + soft thresholds.
// Re-runs whenever threshold params change; doublet detection is unchanged.
// ---------------------------------------------------------------------------

process APPLY_QC {
    label     "process_medium"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/qc/apply_qc/${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5_file), path(dbl_csv)
    path  metadata_csv
    path  genome_gtf
    path  script

    output:
    tuple val(sampleId), path("qc_metrics_${sampleId}.csv.gz"), emit: qc_metrics
    tuple val(sampleId), path("qc_summary_${sampleId}.csv"),    emit: qc_summary

    script:
    """
    set -euo pipefail

    Rscript ${script} \
        "${sampleId}" \
        "${h5_file}" \
        "${genome_gtf}" \
        "${dbl_csv}" \
        "${metadata_csv}" \
        "${params.qc_hard_min_counts}" \
        "${params.qc_hard_min_feats}" \
        "${params.qc_hard_max_mito}" \
        "${params.qc_min_counts}" \
        "${params.qc_min_feats}" \
        "${params.qc_max_mito}" \
        "${params.qc_min_mito}" \
        "${params.qc_max_splice}" \
        "${params.qc_min_splice}"
    """
}
