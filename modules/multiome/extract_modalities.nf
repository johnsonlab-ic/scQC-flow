// Extract Gene Expression and ATAC modalities from multiome H5 files
// GEX H5 is required before running CellBender on multiome data
// ATAC H5 is passed to Seurat for creating the ChromatinAssay

process EXTRACT_MODALITIES {
    label "process_low"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image:latest"
    // No publishDir - intermediate files not needed in final output
    
    input:
    tuple val(sampleName), path(mappingDir), path(extract_script)

    output:
    tuple val(sampleName), path("${sampleName}_gex_raw.h5"), emit: gex_h5
    tuple val(sampleName), path("${sampleName}_atac_raw.h5"), emit: atac_h5
    tuple val(sampleName), path(mappingDir), emit: mapping_dir

    script:
    """
    echo "Extracting modalities from multiome H5 for sample: ${sampleName}"
    
    Rscript ${extract_script} \\
        --input_h5 ${mappingDir}/outs/raw_feature_bc_matrix.h5 \\
        --output_gex_h5 ${sampleName}_gex_raw.h5 \\
        --output_atac_h5 ${sampleName}_atac_raw.h5 \\
        --sample_name ${sampleName}
    
    echo "Modality extraction completed for ${sampleName}"
    """
}
