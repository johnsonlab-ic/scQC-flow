// Seurat module: create Seurat objects pre- and post-QC and add DropletQC/scDbl metadata

process CREATE_SEURAT {
    label "process_seurat"
    tag { sampleName }
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true

  input:
  tuple val(sampleName), path(mappingDir), path(dropletqc_metrics), path(scdbl_metrics), path(seurat_script)

  output:
  tuple val(sampleName), path(mappingDir), path("${sampleName}_seurat_object.rds"), path("${sampleName}_seurat_object_postqc.rds")

  script:
  """
  echo "Creating Seurat objects for sample: ${sampleName}"
  echo "Mapping dir: ${mappingDir}"
  echo "DropletQC metrics: ${dropletqc_metrics}"
  echo "scDbl metrics: ${scdbl_metrics}"

  # Run the external R script with (sample, mappingDir, dropletqc, scdbl)
  Rscript ${seurat_script} "${sampleName}" "${mappingDir}" "${dropletqc_metrics}" "${scdbl_metrics}"

  echo "Seurat objects created for ${sampleName}"
  """
}
