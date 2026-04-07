// =============================================================================
// ATAC QC REPORT FOR MULTIOME DATA
// =============================================================================
process GENERATE_ATAC_REPORT {
    label "process_reports"
    tag "$sampleName"
    container "ghcr.io/johnsonlab-ic/landmark-singlecell-atac:with-radian"
    publishDir "${params.outputDir}/reports/atac", mode: 'copy', overwrite: true

    input:
    tuple val(sampleName), path(seurat_rds), path(fragment_file), path(peak_file), path(template_qmd)

    output:
    tuple val(sampleName), path("${sampleName}_atac_qc_report.html"), emit: html_report
    tuple val(sampleName), path("${sampleName}_atac_qc_report.qmd"), emit: qmd_source

    script:
    """
    echo "Generating ATAC QC report for sample: ${sampleName}"
    echo "Seurat RDS: ${seurat_rds}"
    echo "Fragment file: ${fragment_file}"
    echo "Peak file: ${peak_file}"

    # Copy the template from input path to a new file
    cp ${template_qmd} ${sampleName}_atac_qc_report.qmd

    # Get absolute paths for replacement
    SEURAT_PATH=\$(realpath ${seurat_rds})
    FRAGMENT_PATH=\$(realpath ${fragment_file})
    PEAK_PATH=\$(realpath ${peak_file})

    # Replace placeholders with actual values
    sed -i "s|SAMPLE_NAME_PLACEHOLDER|${sampleName}|g" ${sampleName}_atac_qc_report.qmd
    sed -i "s|SEURAT_PATH_PLACEHOLDER|\$SEURAT_PATH|g" ${sampleName}_atac_qc_report.qmd
    sed -i "s|FRAGMENT_PATH_PLACEHOLDER|\$FRAGMENT_PATH|g" ${sampleName}_atac_qc_report.qmd
    sed -i "s|PEAK_PATH_PLACEHOLDER|\$PEAK_PATH|g" ${sampleName}_atac_qc_report.qmd

    # Render the report
    echo "Rendering ATAC QC Quarto report..."
    export HOME="\$PWD"
    quarto render ${sampleName}_atac_qc_report.qmd

    echo "ATAC QC report completed for ${sampleName}"
    """
}
