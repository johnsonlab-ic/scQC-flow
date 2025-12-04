// DropletQC module for nuclear fraction analysis
// This module provides dropletQC analysis for Cell Ranger outputs

process DROPLETQC {
    label "process_dropletqc"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    // No publishDir - intermediate files not needed in final output

    input:
    tuple val(sampleName), path(bamFile), path(bamIndex), path(barcodesFile), path(run_dropletqc_R)

    output:
    tuple val(sampleName), path("${sampleName}_dropletqc_metrics.csv"), emit: metrics
    tuple val(sampleName), path("${sampleName}_dropletqc_summary.txt"), emit: summary

    script:
    """
    echo "Running DropletQC analysis for sample: ${sampleName}"
    echo "BAM file: ${bamFile}"
    echo "BAM index: ${bamIndex}"
    echo "Barcodes file: ${barcodesFile}"

    # Copy BAM index into workdir (some compute environments don't stage sidecar files in the same dir)
    echo "Copying BAM index into workdir"
    touch ${bamIndex} .

    # Run the R script with arguments
    echo "Executing DropletQC analysis with 20 cores..."
    Rscript ${run_dropletqc_R} --bam_file ${bamFile} --barcodes_file ${barcodesFile} --sample_name ${sampleName} --cores 20

    echo "DropletQC analysis completed for ${sampleName}"
    """
}
