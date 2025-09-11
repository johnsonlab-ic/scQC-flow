// Seurat module: create Seurat objects pre- and post-QC and add DropletQC/scDbl metadata

process CREATE_SEURAT {
    label "process_seurat"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true

  input:
  tuple val(sampleName), path(mappingDir), path(dropletqc_metrics), path(scdbl_metrics), path(seurat_script), val(max_mito), val(min_nuclear), val(metadata_file)

  output:
  tuple val(sampleName), path("${sampleName}_seurat_object.rds"), path("${sampleName}_seurat_object_postqc.rds")

  script:
  """
  echo "Creating Seurat objects for sample: ${sampleName}"
  echo "Mapping dir: ${mappingDir}"
  echo "DropletQC metrics: ${dropletqc_metrics}"
  echo "scDbl metrics: ${scdbl_metrics}"
  echo "QC thresholds: max_mito=${max_mito}, min_nuclear=${min_nuclear}"
  if [ ! -z "${metadata_file}" ] && [ "${metadata_file}" != "null" ]; then
    echo "Using metadata file: ${metadata_file}"
    METADATA_ARG="--metadata ${metadata_file}"
  else
    METADATA_ARG=""
  fi

  # Run the external R script with QC parameters and optional metadata
  Rscript ${seurat_script} "${sampleName}" "${mappingDir}" "${dropletqc_metrics}" "${scdbl_metrics}" \
    --max_mito ${max_mito} --min_nuclear ${min_nuclear} $METADATA_ARG

  echo "Seurat objects created for ${sampleName}"
  """
}
