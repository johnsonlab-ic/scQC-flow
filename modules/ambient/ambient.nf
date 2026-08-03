// ambient.nf — ambient RNA removal

process DECONTX {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
  publishDir "${params.outputDir}/ambient/decontx/${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5_file), path(knee_csv)
    path  script

    output:
    tuple val(sampleId), path("filt_counts_${sampleId}.h5"),                emit: h5_files
    tuple val(sampleId), path("cell_barcodes_${sampleId}.csv"),             emit: barcodes
    tuple val(sampleId), path("dcx_params_${sampleId}.csv"),                emit: dcx_params
    tuple val(sampleId), path("barcodes_qc_metrics_${sampleId}.csv.gz"),    emit: qc_metrics
    tuple val(sampleId), path("dcx_summary_${sampleId}.csv"),               emit: summaries
    tuple val(sampleId),
          env('DCX_EXPECTED_CELLS'),
          env('DCX_TOTAL_INCLUDED'),
          env('DCX_MEAN_CONTAMINATION'),
          emit: ambient_params

    script:
    """
    set -euo pipefail

    Rscript ${script} "${sampleId}" "${h5_file}" "${knee_csv}"

    source "ambient_params_${sampleId}.env"
    export DCX_EXPECTED_CELLS
    export DCX_TOTAL_INCLUDED
    export DCX_MEAN_CONTAMINATION
    """
}

process CELLBENDER {
    label     "process_cellbender"
    tag       "$sampleId"
  publishDir "${params.outputDir}/ambient/cellbender/${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(h5_file), path(knee_csv)
    path  postprocess_script

    output:
    tuple val(sampleId), path("filt_counts_${sampleId}.h5"),             emit: h5_files
    tuple val(sampleId), path("cell_barcodes_${sampleId}.csv"),          emit: barcodes
    tuple val(sampleId), path("cb_summary_${sampleId}.csv"),             emit: summaries
    tuple val(sampleId), path("cb_barcode_labels_${sampleId}.csv.gz"),   emit: labels

    script:
    def migDevice = params.cellbender_mig_device?.trim() ?: ''
    def cbEnvName = params.cellbender_env_name?.trim() ?: ''
    """
    set -euo pipefail

    ${migDevice ? "export CUDA_VISIBLE_DEVICES=\"${migDevice}\"" : '# CUDA_VISIBLE_DEVICES not set; PyTorch will use default GPU'}

    if [[ -n "${cbEnvName}" ]]; then
      CB_PREFIX=(conda run --no-capture-output -n "${cbEnvName}")
    else
      CB_PREFIX=()
    fi

    python3 - <<'PY' "${knee_csv}" > cb_params.env
import csv
import sys

knee_csv = sys.argv[1]
with open(knee_csv, newline='') as fh:
    rows = list(csv.DictReader(fh))
if not rows:
    raise SystemExit("knee CSV is empty")

row0 = rows[0]
def as_int(v, default=0):
    if v is None or v == "":
        return default
    return int(float(v))

print(f"CB_EXPECTED_CELLS={as_int(row0.get('expected_cells', '0'))}")
print(f"CB_TOTAL_DROPLETS_INCLUDED={as_int(row0.get('total_droplets_included', '0'))}")
print(f"CB_LOW_COUNT_THRESHOLD={as_int(row0.get('low_count_threshold', '5'), 5)}")
PY

    source cb_params.env

    if [[ "${params.cellbender_expected_cells}" != "0" ]]; then
        CB_EXPECTED_CELLS=${params.cellbender_expected_cells}
    fi
    if [[ "${params.cellbender_total_droplets_included}" != "0" ]]; then
        CB_TOTAL_DROPLETS_INCLUDED=${params.cellbender_total_droplets_included}
    fi
    if [[ "${params.cellbender_low_count_threshold}" != "0" ]]; then
        CB_LOW_COUNT_THRESHOLD=${params.cellbender_low_count_threshold}
    fi

    CB_OUT="cellbender_${sampleId}.h5"

    CMD=(
      \${CB_PREFIX[@]} cellbender remove-background
      --input "${h5_file}"
      --output "\$CB_OUT"
      --total-droplets-included "\$CB_TOTAL_DROPLETS_INCLUDED"
      --low-count-threshold "\$CB_LOW_COUNT_THRESHOLD"
      --learning-rate "${params.cellbender_learning_rate}"
      --empty-drop-training-fraction "${params.cellbender_empty_drop_training_fraction}"
      --checkpoint-mins 30.0
    )

    if [[ "${params.cellbender_posterior_batch_size}" != "0" ]]; then
      CMD+=(--posterior-batch-size "${params.cellbender_posterior_batch_size}")
    fi

    if [[ "\$CB_EXPECTED_CELLS" -gt 0 ]]; then
      CMD+=(--expected-cells "\$CB_EXPECTED_CELLS")
    fi

    CMD+=(--cuda)

    "\${CMD[@]}"

    if [[ ! -f "cellbender_${sampleId}_filtered.h5" ]]; then
      echo "Expected filtered output not found: cellbender_${sampleId}_filtered.h5" >&2
      exit 1
    fi

    cp "cellbender_${sampleId}_filtered.h5" "filt_counts_${sampleId}.h5"

    python3 "${postprocess_script}" \
      --sample_id "${sampleId}" \
      --knee_csv "${knee_csv}" \
      --filtered_h5 "filt_counts_${sampleId}.h5" \
      --barcodes_out "cell_barcodes_${sampleId}.csv" \
      --summary_out "cb_summary_${sampleId}.csv" \
      --labels_out "cb_barcode_labels_${sampleId}.csv.gz"
    """
}
