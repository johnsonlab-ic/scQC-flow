#!/bin/bash

# Test script for scQC-flow with Quarto reporting
# Usage: ./test_report.sh [sample_name]

SAMPLE=${1:-"test_sample"}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Testing scQC-flow with Quarto reporting for sample: $SAMPLE"
echo "Working directory: $SCRIPT_DIR"

# Create a simple test CSV file for one sample
cat > test_mapping.csv << EOF
samplename,path
$SAMPLE,/rds/general/user/ah3918/ephemeral/Pilot1/$SAMPLE/mapping
EOF

echo "Created test mapping file:"
cat test_mapping.csv

echo ""
echo "Running Nextflow pipeline with report generation..."
echo "Command: nextflow run main.nf --mapping_dirs test_mapping.csv --report -profile imperial"

# Run the pipeline
nextflow run main.nf \
  --mapping_dirs test_mapping.csv \
  --report \
  -profile imperial \
  -resume

echo ""
echo "Pipeline execution completed."
echo "Check the work/ directory for outputs and any HTML reports generated."
