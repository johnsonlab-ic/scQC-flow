// ambient_de.nf — ambient RNA differential expression analysis
//
// Runs EdgeR comparison of empty droplets vs. cells to identify ambient genes.
// Pools all samples into a multi-sample comparison, building pseudobulk matrices
// from raw H5 (empty barcodes) and decontX-filtered H5 (cell barcodes).

// ---------------------------------------------------------------------------
// Lightweight staging step: symlink af_counts_mat.h5 → barcode_matrix_<id>.h5
// so that collected files have unique names for AMBIENT_DE.
// ---------------------------------------------------------------------------

process STAGE_RAW_H5 {
    label     "process_low"
    tag       "$sampleId"

    input:
    tuple val(sampleId), path(h5_file)

    output:
    path "barcode_matrix_${sampleId}.h5", emit: h5

    script:
    """
    ln -s ${h5_file} barcode_matrix_${sampleId}.h5
    """
}

process AMBIENT_DE {
    label     "process_high"
    tag       "ambient_de"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/ambient_de", mode: 'copy', overwrite: true

    input:
    path raw_h5_files      // collected barcode_matrix_*.h5 (uniquely named via STAGE_RAW_H5)
    path knee_csvs         // collected knee_plot_data_*.csv from BARCODE_ESTIMATION
    path filt_h5_files     // collected filt_counts_*.h5 from DECONTX
    path genome_gtf        // reference GTF
    path script            // ambient_de.R

    output:
    path "edger_dt.csv.gz",   emit: de_table
    path "pb_empties.rds",    emit: pb_empties

    script:
    """
    set -euo pipefail

    Rscript ${script} \
        --raw_h5_pattern  'barcode_matrix_*.h5' \
        --knee_pattern    'knee_plot_data_*.csv' \
        --filt_h5_pattern 'filt_counts_*.h5' \
        --gtf             ${genome_gtf} \
        --out_de          edger_dt.csv.gz \
        --out_pb_empties  pb_empties.rds
    """
}
