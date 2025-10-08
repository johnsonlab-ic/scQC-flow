// DropletQC module for nuclear fraction analysis
// This module provides dropletQC analysis for Cell Ranger outputs

process DROPLETQC {
    label "process_dropletqc"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true

    input:
    tuple val(sampleName), path(bamFile), path(barcodesFile), path(run_dropletqc_R)

    output:
    tuple val(sampleName), path("${sampleName}_dropletqc_metrics.csv"), emit: metrics
    tuple val(sampleName), path("${sampleName}_dropletqc_summary.txt"), emit: summary

    script:
    """
    echo "Running DropletQC analysis for sample: ${sampleName}"
    echo "BAM file: ${bamFile}"
    echo "Barcodes file: ${barcodesFile}"

    # Run the R script with arguments
    echo "Executing DropletQC analysis with 20 cores..."
    Rscript ${run_dropletqc_R} --bam_file ${bamFile} --barcodes_file ${barcodesFile} --sample_name ${sampleName} --cores 20

    echo "DropletQC analysis completed for ${sampleName}"
    """
}
