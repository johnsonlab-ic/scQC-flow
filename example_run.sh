#!/usr/bin/env bash
set -euo pipefail

# Example run using bundled mini data.
# Edit the paths below before running.

# ---------------------------------------------------------------------------
# Paths — edit to match your environment
# ---------------------------------------------------------------------------
CELLRANGER_PATH="/mnt/data/projects/scQC-flow/resources/downloads/cellranger-10.0.0/"
GENOME_FASTA="/mnt/data/projects/scQC-flow/resources/downloads/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
GENOME_GTF="/mnt/data/projects/scQC-flow/resources/downloads/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"

# ---------------------------------------------------------------------------
# Output locations
# ---------------------------------------------------------------------------
PROJECT="data/test_run_mapping"
WORK_DIR="${PROJECT}/work"
OUTDIR="${PROJECT}/outputs"

rm -rf "${WORK_DIR}" "${OUTDIR}"
mkdir -p "${WORK_DIR}" "${OUTDIR}/pipeline_info"

export NXF_LOG_FILE="${OUTDIR}/pipeline_info/nextflow.log"

echo "Starting scQC-flow MAPPING workflow"

nextflow run main.nf \
  -profile offline \
  -w "${WORK_DIR}" \
  --raw_data_dir    data/mini_landmark_data/raw \
  --cellrangerPath  "${CELLRANGER_PATH}" \
  --genome_fasta    "${GENOME_FASTA}" \
  --genome_gtf      "${GENOME_GTF}" \
  --chemistry       3v4 \
  --outputDir       "${OUTDIR}"

echo ""
echo "Run complete. Outputs:"
find "${OUTDIR}" -maxdepth 3 -type f | sort
