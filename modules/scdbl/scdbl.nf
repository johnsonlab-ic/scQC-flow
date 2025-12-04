// scDblFinder module for doublet detection
// This module provides scDblFinder analysis for Cell Ranger outputs

process SCDBL {
    label "process_dropletqc"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    // No publishDir - intermediate files not needed in final output

    input:
    tuple val(sampleName), path(h5File), path(run_scdbl_R)

    output:
    tuple val(sampleName), path("${sampleName}_scdbl_metrics.csv"), emit: metrics
    tuple val(sampleName), path("${sampleName}_scdbl_summary.txt"), emit: summary

    script:
    """
    echo "Running scDblFinder analysis for sample: ${sampleName}"
    echo "H5 file: ${h5File}"

    # Run the R script with arguments
    echo "Executing scDblFinder doublet detection..."
    Rscript ${run_scdbl_R} --h5_file ${h5File} --sample_name ${sampleName}

    echo "scDblFinder analysis completed for ${sampleName}"
    """
}
