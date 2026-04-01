// RNA mapping module
// Produces Cell Ranger-style mapped directories for downstream QC.

process MAP_CELLRANGER {
    label "process_higher_memory"
    tag "$sampleId"
    publishDir "${params.outputDir}/mapping", mode: 'copy', overwrite: true

    input:
    tuple val(sampleId), val(sampleName), path(fastqPath), val(cellrangerPath), val(transcriptome)

    output:
    tuple val(sampleId), path("${sampleId}_mapped"), emit: mapped_dirs

    script:
    """
    set -euo pipefail

    # Resolve the cellranger binary: accept either a directory or a direct path to the binary.
    if [ -x "${cellrangerPath}/cellranger" ]; then
      CELLRANGER_BIN="${cellrangerPath}/cellranger"
    elif [ -x "${cellrangerPath}" ]; then
      CELLRANGER_BIN="${cellrangerPath}"
    else
      echo "ERROR: Cannot find cellranger binary at ${cellrangerPath}"
      exit 1
    fi

    if [ ! -d "${transcriptome}" ]; then
      echo "ERROR: Transcriptome directory does not exist: ${transcriptome}"
      exit 1
    fi

    echo "Mapping sample ${sampleId}"
    echo "FASTQ path: ${fastqPath}"
    echo "Cell Ranger binary: \$CELLRANGER_BIN"
    echo "Transcriptome: ${transcriptome}"

    "\$CELLRANGER_BIN" count \
      --id="${sampleId}_mapped" \
      --create-bam true \
      --fastqs="${fastqPath}" \
      --sample="${sampleName}" \
      --transcriptome="\$TRANSCRIPTOME_DIR"

    # Keep only outs/ to reduce footprint in published mapping outputs.
    if [ -d "${sampleId}_mapped" ]; then
      find "${sampleId}_mapped" -mindepth 1 -maxdepth 1 ! -name outs -exec rm -rf {} + || true
    fi
    """
}
