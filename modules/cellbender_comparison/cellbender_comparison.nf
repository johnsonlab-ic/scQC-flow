// CellBender vs Cell Ranger Comparison Module
// Compares droplet calling between Cell Ranger and CellBender
// Generates metrics CSV and knee-plot visualization

process CELLBENDER_COMPARISON {
    label "process_medium"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    publishDir "${params.outputDir}/${sampleName}/cellbender_comparison", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir), path(cellbender_h5)

    output:
    tuple val(sampleName), path("${sampleName}_cellbender_comparison_metrics.csv"), path("${sampleName}_cellbender_comparison_kneeplot.png"), emit: comparison_results

    script:
    """
    echo "Running CellBender vs Cell Ranger comparison for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"
    echo "CellBender H5: ${cellbender_h5}"

    Rscript ${projectDir}/modules/cellbender_comparison/cellbender_comparison.R \\
        --raw_h5 "${mappingDir}/outs/raw_feature_bc_matrix.h5" \\
        --filtered_h5 "${mappingDir}/outs/filtered_feature_bc_matrix.h5" \\
        --cellbender_h5 "${cellbender_h5}" \\
        --output_metrics "${sampleName}_cellbender_comparison_metrics.csv" \\
        --output_plot "${sampleName}_cellbender_comparison_kneeplot.png" \\
        --sample_name "${sampleName}"

    echo "CellBender comparison completed for ${sampleName}"
    """
}

process CELLBENDER_COMPARISON_STATS_ONLY {
    label "process_medium"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    publishDir "${params.outputDir}/${sampleName}/cellbender_comparison", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir)

    output:
    tuple val(sampleName), path("${sampleName}_cellbender_comparison_metrics.csv"), path("${sampleName}_cellbender_comparison_kneeplot.png"), emit: comparison_results

    script:
    """
    echo "Running droplet calling metrics for sample: ${sampleName} (Cell Ranger only)"
    echo "Mapping directory: ${mappingDir}"

    Rscript ${projectDir}/modules/cellbender_comparison/cellbender_comparison.R \\
        --raw_h5 "${mappingDir}/outs/raw_feature_bc_matrix.h5" \\
        --filtered_h5 "${mappingDir}/outs/filtered_feature_bc_matrix.h5" \\
        --output_metrics "${sampleName}_cellbender_comparison_metrics.csv" \\
        --output_plot "${sampleName}_cellbender_comparison_kneeplot.png" \\
        --sample_name "${sampleName}"

    echo "Droplet calling metrics completed for ${sampleName}"
    """
}
