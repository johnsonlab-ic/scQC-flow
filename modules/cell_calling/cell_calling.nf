// cell_calling.nf — splice-aware per-sample cell selection
//
// Replaces decontX's cell-calling role. Classifies barcodes with an anchored
// per-sample GMM in (log10 UMI, logit splice%) plus hard biological guardrails,
// and emits a drop-in filtered cell matrix (UNCORRECTED) for the downstream
// QC/HVG contract, plus the is_empty barcode set for ambient-profile estimation.

process CELL_CALLING {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/cell_calling", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5_file), path(knee_csv)
    path  script

    output:
    tuple val(sampleId), path("filt_counts_${sampleId}.h5"),             emit: h5_files
    tuple val(sampleId), path("cell_barcodes_${sampleId}.csv"),          emit: barcodes
    tuple val(sampleId), path("empty_barcodes_${sampleId}.csv"),         emit: empties
    tuple val(sampleId), path("cell_calling_labels_${sampleId}.csv.gz"), emit: labels
    tuple val(sampleId), path("cell_calling_summary_${sampleId}.csv"),   emit: summaries
    tuple val(sampleId), path("cell_calling_gmm_${sampleId}.rds"),       emit: gmm
    tuple val(sampleId), path("cell_calling_${sampleId}.env"),           emit: params_env

    script:
    def incl = ((params.cell_calling_include_damaged ?: false).toString().toBoolean()) ? 'TRUE' : 'FALSE'
    def alpha_dir = (params.cell_calling_alpha_dir ?: '').toString()
    def alpha_csv = alpha_dir ? "${alpha_dir}/${sampleId}_alpha_hat.csv.gz" : ''
    def alpha_max = params.cell_calling_alpha_max ?: 0.5
    """
    Rscript ${script} "${sampleId}" "${h5_file}" "${knee_csv}" ${incl} "${alpha_csv}" ${alpha_max}
    """
}
