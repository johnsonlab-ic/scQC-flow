// Seurat module for Multiome data: create Seurat objects pre- and post-QC
// Extracts Gene Expression modality from multiome H5 files
// Also copies ATAC files for downstream analysis

process CREATE_SEURAT_MULTIOME {
    label "process_seurat"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true

  input:
  tuple val(sampleName), path(mappingDir), path(dropletqc_metrics), path(scdbl_metrics), path(seurat_script), val(max_mito), val(min_nuclear), val(metadata_file), path(h5_path)

  output:
  tuple val(sampleName), path("${sampleName}_seurat_object.rds"), path("${sampleName}_seurat_object_postqc.rds")
  path "atac/", emit: atac_files

  script:
  """
  echo "Creating Seurat objects for multiome sample: ${sampleName}"
  echo "Mapping dir: ${mappingDir}"
  echo "DropletQC metrics: ${dropletqc_metrics}"
  echo "scDbl metrics: ${scdbl_metrics}"
  echo "QC thresholds: max_mito=${max_mito}, min_nuclear=${min_nuclear}"
  echo "H5 counts file: ${h5_path}"

  # Run the external R script with QC parameters, H5 path, and optional metadata
  if [ ! -z "${metadata_file}" ] && [ "${metadata_file}" != "null" ]; then
    echo "Using metadata file: ${metadata_file}"
    Rscript ${seurat_script} "${sampleName}" "${mappingDir}" "${dropletqc_metrics}" "${scdbl_metrics}" \
      --max_mito ${max_mito} --min_nuclear ${min_nuclear} --metadata "${metadata_file}" --h5_path "${h5_path}"
  else
    Rscript ${seurat_script} "${sampleName}" "${mappingDir}" "${dropletqc_metrics}" "${scdbl_metrics}" \
      --max_mito ${max_mito} --min_nuclear ${min_nuclear} --h5_path "${h5_path}"
  fi

  echo "Seurat objects created for ${sampleName}"

  # Copy ATAC files for downstream multiome analysis
  echo "Copying ATAC files..."
  mkdir -p atac
  cp ${mappingDir}/outs/atac_* atac/ 2>/dev/null || echo "No atac_* files found"
  
  # List what was copied
  echo "ATAC files copied:"
  ls -la atac/

  echo "Multiome Seurat processing completed for ${sampleName}"
  """
}
