// CellBender module for ambient RNA removal and empty droplet filtering
// This module provides both CPU and GPU variants of CellBender

process CELLBENDER {
    label "process_cellbender"
    tag { sampleName }
    container "us.gcr.io/broad-dsde-methods/cellbender:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir)

    output:
    tuple val(sampleName), path("${sampleName}_cellbender_output"), emit: cellbender_output
    path "${sampleName}_cellbender_output/cellbender_out.h5", emit: h5_file
    path "${sampleName}_cellbender_output/summary.txt", emit: summary

    script:
    """
    echo "Running CellBender (CPU) for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"

    # Extract estimated number of cells from metrics file
    ESTIMATED_CELLS=\$(awk -F'"' 'NR==2 {gsub(/,/, "", \$2); print \$2}' ${mappingDir}/outs/metrics_summary.csv)
    echo "Estimated cells from metrics: \$ESTIMATED_CELLS"

    mkdir -p ${sampleName}_cellbender_output

    cellbender remove-background \\
                 --input ${mappingDir}/outs/raw_feature_bc_matrix.h5 \\
                 --output ${sampleName}_cellbender_output/cellbender_out.h5 \\
                 --expected-cells \$ESTIMATED_CELLS

    echo "CellBender processing completed" > ${sampleName}_cellbender_output/summary.txt
    echo "CellBender (CPU) completed for ${sampleName}"
    """
}

process CELLBENDER_GPU {
    label "process_gpu"
    tag { sampleName }
    container "us.gcr.io/broad-dsde-methods/cellbender:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir)

    output:
    tuple val(sampleName), path("${sampleName}_cellbender_output"), emit: cellbender_output
    path "${sampleName}_cellbender_output/cellbender_out.h5", emit: h5_file
    path "${sampleName}_cellbender_output/summary.txt", emit: summary

    script:
    """
    echo "Running CellBender (GPU) for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"

    # Extract estimated number of cells from metrics file
    ESTIMATED_CELLS=\$(awk -F'"' 'NR==2 {gsub(/,/, "", \$2); print \$2}' ${mappingDir}/outs/metrics_summary.csv)
    echo "Estimated cells from metrics: \$ESTIMATED_CELLS"

    mkdir -p ${sampleName}_cellbender_output

    cellbender remove-background \\
                 --input ${mappingDir}/outs/raw_feature_bc_matrix.h5 \\
                 --output ${sampleName}_cellbender_output/cellbender_out.h5 \\
                 --cuda \\
                 --expected-cells \$ESTIMATED_CELLS

    echo "CellBender processing completed" > ${sampleName}_cellbender_output/summary.txt
    echo "CellBender (GPU) completed for ${sampleName}"
    """
}

process CELLBENDER_H5_CONVERT {
    label "process_low"
    tag { sampleName }
    container "ah3918/pilot-analyses:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(cellbender_output)

    output:
    tuple val(sampleName), path("${sampleName}_cellbender_output_seurat.h5"), emit: seurat_h5

    script:
    """
    echo "Converting CellBender H5 to Seurat-compatible format for: ${sampleName}"
    
    # Use ptrepack to create Seurat-compatible H5 file, overwriting nodes if needed
    ptrepack --complevel 5 ${cellbender_output}/cellbender_out_filtered.h5:/matrix ${sampleName}_cellbender_output_seurat.h5:/matrix
    
    echo "H5 conversion completed for ${sampleName}"
    """
}
