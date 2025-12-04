// CellBender module for multiome data
// Uses extracted Gene Expression H5 instead of raw multiome H5

process CELLBENDER_MULTIOME {
    label "process_cellbender"
    tag { sampleName }
    container "us.gcr.io/broad-dsde-methods/cellbender:latest"
    // No publishDir - intermediate files not needed in final output
    
    input:
    tuple val(sampleName), path(gex_h5)

    output:
    tuple val(sampleName), path("${sampleName}_cellbender_output"), emit: cellbender_output
    path "${sampleName}_cellbender_output/cellbender_out.h5", emit: h5_file
    path "${sampleName}_cellbender_output/summary.txt", emit: summary

    script:
    """
    echo "Running CellBender (CPU) for multiome sample: ${sampleName}"
    echo "Input GEX H5: ${gex_h5}"

    mkdir -p ${sampleName}_cellbender_output

    cellbender remove-background \\
                 --input ${gex_h5} \\
                 --output ${sampleName}_cellbender_output/cellbender_out.h5

    echo "CellBender processing completed" > ${sampleName}_cellbender_output/summary.txt
    echo "CellBender (CPU) completed for ${sampleName}"
    """
}

process CELLBENDER_MULTIOME_GPU {
    label "process_gpu"
    tag { sampleName }
    container "us.gcr.io/broad-dsde-methods/cellbender:latest"
    // No publishDir - intermediate files not needed in final output
    
    input:
    tuple val(sampleName), path(gex_h5)

    output:
    tuple val(sampleName), path("${sampleName}_cellbender_output"), emit: cellbender_output
    path "${sampleName}_cellbender_output/cellbender_out.h5", emit: h5_file
    path "${sampleName}_cellbender_output/summary.txt", emit: summary

    script:
    """
    echo "Running CellBender (GPU) for multiome sample: ${sampleName}"
    echo "Input GEX H5: ${gex_h5}"

    mkdir -p ${sampleName}_cellbender_output

    cellbender remove-background \\
                 --input ${gex_h5} \\
                 --output ${sampleName}_cellbender_output/cellbender_out.h5 \\
                 --cuda

    echo "CellBender processing completed" > ${sampleName}_cellbender_output/summary.txt
    echo "CellBender (GPU) completed for ${sampleName}"
    """
}
