// cellsweep.nf — CellSweep ambient/bulk decontamination (post-annotation refinement)
//
// Per sample: builds the CellSweep input from pipeline outputs and runs the EM
// denoiser. Inputs are the post-QC singlet nuclei (with their integration Leiden
// cluster as `celltype`) + a subsample of is_empty barcodes (raw counts) for the
// ambient profile. Outputs the decontaminated AnnData + per-cell alpha_hat.

process CELLSWEEP {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/cellsweep", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(filt_h5), path(raw_h5), path(empty_csv)
    path  integration_dt
    path  script

    output:
    tuple val(sampleId), path("${sampleId}_alpha_hat.csv.gz"), emit: alpha
    tuple val(sampleId), path("${sampleId}_cellsweep.h5ad"),   emit: h5ad

    script:
    def n_emp = params.cellsweep_n_empties     ?: 30000
    def ctcol = params.cellsweep_celltype_col  ?: 'RNA_snn_res.0.5'
    """
    # writable HOME + numba cache (task runs as a UID with no home dir; scanpy/numba
    # imported by cellsweep otherwise fail with 'no locator available' cache errors)
    export HOME="\$PWD"
    export NUMBA_CACHE_DIR="\$PWD/.numba_cache"
    export MPLCONFIGDIR="\$PWD/.mpl"

    python3 ${script} \\
      --sample_id ${sampleId} \\
      --filt_h5 ${filt_h5} \\
      --raw_h5 ${raw_h5} \\
      --empty_barcodes ${empty_csv} \\
      --integration_dt ${integration_dt} \\
      --out_dir . \\
      --n_empties ${n_emp} \\
      --celltype_col ${ctcol} \\
      --keep_h5ad
    """
}
