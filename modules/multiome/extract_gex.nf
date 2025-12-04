// Extract Gene Expression modality from multiome H5 files
// This is required before running CellBender on multiome data

process EXTRACT_GEX_H5 {
    label "process_low"
    tag { sampleName }
    container "ah3918/pilot-analyses:latest"
    // No publishDir - intermediate files not needed in final output
    
    input:
    tuple val(sampleName), path(mappingDir), path(extract_script)

    output:
    tuple val(sampleName), path("${sampleName}_gex_raw.h5"), emit: gex_h5
    tuple val(sampleName), path(mappingDir), emit: mapping_dir

    script:
    """
    echo "Extracting Gene Expression from multiome H5 for sample: ${sampleName}"
    
    Rscript ${extract_script} \\
        --input_h5 ${mappingDir}/outs/raw_feature_bc_matrix.h5 \\
        --output_h5 ${sampleName}_gex_raw.h5 \\
        --sample_name ${sampleName}
    
    echo "Gene Expression extraction completed for ${sampleName}"
    """
}
