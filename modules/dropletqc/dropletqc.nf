// DropletQC module for nuclear fraction analysis
// This module provides dropletQC analysis for Cell Ranger outputs

process DROPLETQC {
    label "process_dropletqc"
    tag { sampleName }
    container "ah3918/pilot-analyses:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir)

    output:
    tuple val(sampleName), path("${sampleName}_dropletqc_metrics.csv"), emit: metrics
    tuple val(sampleName), path("${sampleName}_dropletqc_summary.txt"), emit: summary

    script:
    """
    echo "Running DropletQC analysis for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"

    # Copy the R script and execute it with parameters
    cp ${projectDir}/modules/dropletqc/run_dropletqc.R .

    # Run the R script with arguments
    echo "Executing DropletQC analysis with 10 cores..."
    Rscript run_dropletqc.R --mapping_dir ${mappingDir} --sample_name ${sampleName} --cores 10

    echo "DropletQC analysis completed for ${sampleName}"
    """
}
