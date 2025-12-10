// Seurat module for Multiome data: create Seurat objects pre- and post-QC
// Creates RNA assay from CellBender/10X H5, adds ATAC ChromatinAssay from ATAC H5
// Filters ATAC cells to match CellBender-filtered cells

process CREATE_SEURAT_MULTIOME {
    label "process_seurat"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true

  input:
  tuple val(sampleName), path(mappingDir), path(dropletqc_metrics), path(scdbl_metrics), path(seurat_script), val(max_mito), val(min_nuclear), val(metadata_file), path(h5_path), path(atac_h5_path)

  output:
  tuple val(sampleName), path("${sampleName}_seurat_object.rds"), path("${sampleName}_seurat_object_postqc.rds"), emit: seurat_rds
  path "atac/", emit: atac_files

  script:
  """
  echo "Creating Seurat objects for multiome sample: ${sampleName}"
  echo "Mapping dir: ${mappingDir}"
  echo "DropletQC metrics: ${dropletqc_metrics}"
  echo "scDbl metrics: ${scdbl_metrics}"
  echo "QC thresholds: max_mito=${max_mito}, min_nuclear=${min_nuclear}"
  echo "H5 counts file (GEX): ${h5_path}"
  echo "ATAC H5 file: ${atac_h5_path}"

  # Run the external R script with QC parameters, H5 paths, and optional metadata
  if [ ! -z "${metadata_file}" ] && [ "${metadata_file}" != "null" ]; then
    echo "Using metadata file: ${metadata_file}"
    Rscript ${seurat_script} "${sampleName}" "${mappingDir}" "${dropletqc_metrics}" "${scdbl_metrics}" \\
      --max_mito ${max_mito} --min_nuclear ${min_nuclear} --metadata "${metadata_file}" \\
      --h5_path "${h5_path}" --atac_h5_path "${atac_h5_path}"
  else
    Rscript ${seurat_script} "${sampleName}" "${mappingDir}" "${dropletqc_metrics}" "${scdbl_metrics}" \\
      --max_mito ${max_mito} --min_nuclear ${min_nuclear} \\
      --h5_path "${h5_path}" --atac_h5_path "${atac_h5_path}"
  fi

  echo "Seurat objects created for ${sampleName}"

  # Copy ATAC files for downstream multiome analysis (excluding BAM files to save space)
  echo "Copying ATAC files (excluding BAM files)..."
  mkdir -p atac
  for f in ${mappingDir}/outs/atac_*; do
    if [[ ! "\$f" == *bam* ]]; then
      cp "\$f" atac/ 2>/dev/null || true
    fi
  done
  
  # List what was copied
  echo "ATAC files copied:"
  ls -la atac/

  echo "Multiome Seurat processing completed for ${sampleName}"
  """
}
