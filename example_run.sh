#!/usr/bin/env bash
set -euo pipefail

# Example smoke-test run using bundled mini data.
# Usage:
#   bash example_run.sh
#
# Notes:
# - This runs mapping + QC against bundled mini data.
# - Execution artifacts are collected in outputs/pipeline_info/
# - Work directory can be safely deleted after run completion.

# Set up clean project structure
MY_PROJECT="data/test_run_cellranger"
WORK_DIR="${MY_PROJECT}/work"
OUTDIR="${MY_PROJECT}/outputs"
export NXF_LOG_FILE="${OUTDIR}/pipeline_info/nextflow.log"

#clean up from previous runs
rm -rf "${WORK_DIR}" "${OUTDIR}"
mkdir -p "${WORK_DIR}" "${OUTDIR}" "${OUTDIR}/pipeline_info"

echo "Starting scQC-flow test run"
echo "Project directory: ${MY_PROJECT}"
echo "  Work directory: ${WORK_DIR}"
echo "  Output directory: ${OUTDIR}"
echo ""

nextflow run main.nf \
  -profile offline \
  -w "${WORK_DIR}" \
  --run_mode mapping \
  --raw_data_dir data/mini_landmark_data/raw \
  --cellrangerPath /mnt/data/projects/scQC-flow/resources/downloads/cellranger-10.0.0 \
  --transcriptome /mnt/data/projects/scQC-flow/resources/downloads/refdata-gex-GRCh38-2024-A \
  --outputDir "${OUTDIR}"

echo ""
echo "Run complete. Top-level output structure:"
find "${OUTDIR}" -maxdepth 2 -type d | sort

echo ""
echo "Execution artifacts (timeline, report, trace, DAG):"
ls -lh "${OUTDIR}/pipeline_info/" 2>/dev/null || echo "  (pipeline_info not yet generated)"

# =============================================================================
# Example 2: alevin-fry mapping with SIMPLEAF_INDEX (build index on-the-fly)
# =============================================================================
# This uses the pipeline's built-in SIMPLEAF_INDEX process by passing FASTA+GTF.
# Do NOT pass --alevinfry_index or --alevinfry_t2g for this mode.
#
# Required inputs:
#   --alevinfry_fasta
#   --alevinfry_gtf
#   --alevinfry_chemistry
#   --cellrangerPath (used to extract matching whitelist)
#
# Where to get these (10x):
#   1) Download a 10x reference bundle (contains FASTA + GTF):
#      https://www.10xgenomics.com/support/software/cell-ranger/downloads
#   2) Extract and use:
#      FASTA: <ref>/fasta/genome.fa
#      GTF:   <ref>/genes/genes.gtf.gz
#   3) Whitelist is extracted automatically from:
#      <cellranger>/lib/python/cellranger/barcodes/
#      based on --alevinfry_chemistry and written as chemistry_whitelist.txt
#
# You can also use local resources already in this repo under resources/ or dev/sc_data/.

AF_PROJECT="data/test_run_alevinfry"
AF_WORK="${AF_PROJECT}/work"
AF_OUTDIR="${AF_PROJECT}/outputs"
AF_REF="/mnt/data/projects/scQC-flow/resources/downloads/refdata-gex-GRCh38-2024-A"
AF_FASTA="${AF_REF}/fasta/genome.fa"
AF_GTF="${AF_REF}/genes/genes.gtf.gz"
AF_CELLRANGER="/mnt/data/projects/scQC-flow/resources/downloads/cellranger-10.0.0"
AF_CHEMISTRY="3v4"

#clean outs
rm -rf "${AF_WORK}" "${AF_OUTDIR}"
mkdir -p "${AF_WORK}" "${AF_OUTDIR}" "${AF_OUTDIR}/pipeline_info"

echo ""
echo "Starting scQC-flow alevin-fry test run (index build mode)"
echo "Project directory: ${AF_PROJECT}"
echo "  Work directory: ${AF_WORK}"
echo "  Output directory: ${AF_OUTDIR}"
echo "  FASTA: ${AF_FASTA}"
echo "  GTF: ${AF_GTF}"
echo "  Cell Ranger path: ${AF_CELLRANGER}"
echo "  Chemistry: ${AF_CHEMISTRY}"
echo ""

nextflow run main.nf \
  -profile offline \
  -w "${AF_WORK}" \
  --run_mode mapping \
  --mapper alevinfry \
  --raw_data_dir data/mini_landmark_data/raw \
  --cellrangerPath "${AF_CELLRANGER}" \
  --alevinfry_fasta "${AF_FASTA}" \
  --alevinfry_gtf "${AF_GTF}" \
  --alevinfry_chemistry "${AF_CHEMISTRY}" \
  --outputDir "${AF_OUTDIR}"

echo ""
echo "alevin-fry run complete. Output structure:"
find "${AF_OUTDIR}" -maxdepth 2 -type d | sort


##with qc
AF_PROJECT="data/test_run_alevinfry"
AF_WORK="${AF_PROJECT}/work"
AF_OUTDIR="${AF_PROJECT}/outputs"
AF_REF="/mnt/data/projects/scQC-flow/resources/downloads/refdata-gex-GRCh38-2024-A"
AF_FASTA="${AF_REF}/fasta/genome.fa"
AF_GTF="${AF_REF}/genes/genes.gtf.gz"
AF_CELLRANGER="/mnt/data/projects/scQC-flow/resources/downloads/cellranger-10.0.0"
AF_CHEMISTRY="3v4"

#clean outs
rm -rf "${AF_WORK}" "${AF_OUTDIR}"
mkdir -p "${AF_WORK}" "${AF_OUTDIR}" "${AF_OUTDIR}/pipeline_info"

echo ""
echo "Starting scQC-flow alevin-fry test run (index build mode)"
echo "Project directory: ${AF_PROJECT}"
echo "  Work directory: ${AF_WORK}"
echo "  Output directory: ${AF_OUTDIR}"
echo "  FASTA: ${AF_FASTA}"
echo "  GTF: ${AF_GTF}"
echo "  Cell Ranger path: ${AF_CELLRANGER}"
echo "  Chemistry: ${AF_CHEMISTRY}"
echo ""

nextflow run main.nf \
  -profile offline \
  -w "${AF_WORK}" \
  --run_mode both \
  --mapper alevinfry \
  --raw_data_dir data/mini_landmark_data/raw \
  --cellrangerPath "${AF_CELLRANGER}" \
  --alevinfry_fasta "${AF_FASTA}" \
  --alevinfry_gtf "${AF_GTF}" \
  --alevinfry_chemistry "${AF_CHEMISTRY}" \
  --outputDir "${AF_OUTDIR}"

echo ""
echo "alevin-fry run complete. Output structure:"
find "${AF_OUTDIR}" -maxdepth 2 -type d | sort
