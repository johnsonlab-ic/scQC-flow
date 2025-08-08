// Reports module for generating QC reports from Cell Ranger outputs
// This module provides Quarto-based report generation

process GENERATE_REPORTS {
    label "process_reports"
    tag { sampleName }
    container "ah3918/pilot-analyses:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir)

    output:
    tuple val(sampleName), path("${sampleName}_qc_report.html"), emit: html_report
    tuple val(sampleName), path("${sampleName}_qc_report.qmd"), emit: qmd_source

    script:
    """
    echo "Generating QC report for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"

    # Copy the template and replace placeholders
    cp ${projectDir}/templates/qc_report_template.qmd ${sampleName}_qc_report.qmd
    
    # Replace placeholders with actual values
    sed -i "s|SAMPLE_NAME_PLACEHOLDER|${sampleName}|g" ${sampleName}_qc_report.qmd
    sed -i "s|DATA_PATH_PLACEHOLDER|${mappingDir}|g" ${sampleName}_qc_report.qmd

    # Render the report
    echo "Rendering Quarto report..."
    quarto render ${sampleName}_qc_report.qmd

    echo "QC report completed for ${sampleName}"
    """
}
