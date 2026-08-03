// CellBender module for ambient RNA removal and empty droplet filtering
// This module provides both CPU and GPU variants of CellBender

process CELLBENDER {
    label "process_cellbender"
    tag "$sampleName"
    container "us.gcr.io/broad-dsde-methods/cellbender:latest"
    publishDir "${params.outputDir}/cellbender", mode: 'copy', overwrite: true
    
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


    mkdir -p ${sampleName}_cellbender_output

    cellbender remove-background \\
                 --input ${mappingDir}/outs/raw_feature_bc_matrix.h5 \\
                 --output ${sampleName}_cellbender_output/cellbender_out.h5

    echo "CellBender processing completed" > ${sampleName}_cellbender_output/summary.txt
    echo "CellBender (CPU) completed for ${sampleName}"
    """
}

process CELLBENDER_GPU {
    label "process_gpu"
    tag "$sampleName"
    container "us.gcr.io/broad-dsde-methods/cellbender:latest"
    publishDir "${params.outputDir}/cellbender", mode: 'copy', overwrite: true
    
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

    mkdir -p ${sampleName}_cellbender_output

    cellbender remove-background \\
                 --input ${mappingDir}/outs/raw_feature_bc_matrix.h5 \\
                 --output ${sampleName}_cellbender_output/cellbender_out.h5 \\
                 --cuda

    echo "CellBender processing completed" > ${sampleName}_cellbender_output/summary.txt
    echo "CellBender (GPU) completed for ${sampleName}"
    """
}

// CellBender variant for alevin-fry: takes a pre-built H5 and knee parameters
// from BARCODE_ESTIMATION instead of a Cell Ranger mapping directory.
process CELLBENDER_FROM_H5 {
    label "process_cellbender"
    tag "$sampleId"
    container "us.gcr.io/broad-dsde-methods/cellbender:latest"
    publishDir "${params.outputDir}/cellbender", mode: 'copy', overwrite: true

    input:
    tuple val(sampleId),
          path(h5_file),
          val(expected_cells),
          val(total_droplets),
          val(low_count_thr)

    output:
    tuple val(sampleId), path("${sampleId}_cellbender_filtered.h5"), emit: filtered_h5
    tuple val(sampleId), path("${sampleId}_cellbender_barcodes.csv"), emit: barcodes

    script:
    def ec_arg  = (expected_cells as int)  > 0 ? "--expected-cells ${expected_cells}"    : ""
    def td_arg  = (total_droplets as int)  > 0 ? "--total-droplets-included ${total_droplets}" : ""
    def lct_arg = (low_count_thr  as int)  > 0 ? "--low-count-threshold ${low_count_thr}" : ""
    """
    cellbender remove-background \\
        --input  "${h5_file}" \\
        --output "${sampleId}_cellbender_out.h5" \\
        ${ec_arg} ${td_arg} ${lct_arg}

    # Rename the filtered output to our standard name
    cp "${sampleId}_cellbender_out_filtered.h5" "${sampleId}_cellbender_filtered.h5"

    # Extract barcodes from the filtered H5
    python3 -c "
import h5py, csv
with h5py.File('${sampleId}_cellbender_filtered.h5', 'r') as f:
    bcs = [b.decode() for b in f['matrix/barcodes'][:]]
with open('${sampleId}_cellbender_barcodes.csv', 'w', newline='') as fh:
    for bc in bcs:
        fh.write(bc + '\\n')
"
    """
}

process CELLBENDER_H5_CONVERT {
    label "process_low"
    tag "$sampleName"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    publishDir "${params.outputDir}/cellbender", mode: 'copy', overwrite: true
    
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
