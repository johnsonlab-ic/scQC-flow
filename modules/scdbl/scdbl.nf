// scDblFinder module for doublet detection
// This module provides scDblFinder analysis for Cell Ranger outputs

process SCDBL {
    label "process_dropletqc"
    tag { sampleName }
    container "ah3918/pilot-analyses:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir)

    output:
    tuple val(sampleName), path("${sampleName}_scdbl_metrics.csv"), emit: metrics
    tuple val(sampleName), path("${sampleName}_scdbl_summary.txt"), emit: summary

    script:
    """
    echo "Running scDblFinder analysis for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"

    # Copy the R script and execute it with parameters
    cp ${projectDir}/modules/scdbl/run_scdbl.R .

    # Run the R script with arguments
    echo "Executing scDblFinder doublet detection..."
    Rscript run_scdbl.R --mapping_dir ${mappingDir} --sample_name ${sampleName}

    echo "scDblFinder analysis completed for ${sampleName}"
    """
}
