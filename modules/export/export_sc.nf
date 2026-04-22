process EXPORT_SCANPY {
    label     "process_high"
    tag       "export_anndata"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/export", mode: params.publish_mode_nonreport, overwrite: true

    input:
    path h5_files
    path qc_csvs
    path integration_csv
    path annotation_labels_csv
    path genome_gtf
    path script

    output:
    path "scqcflow_all_cells.h5ad", emit: all_cells
    path "scqcflow_clean_cells.h5ad", emit: clean_cells

    script:
    def annotationArg = annotation_labels_csv.name == 'NO_FILE' ? 'NO_FILE' : annotation_labels_csv.name
    """
    set -euo pipefail
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export NUMBA_CACHE_DIR="\$PWD/.numba"
    export PYTHONUNBUFFERED=1

    python3 -u ${script} \
        --h5_pattern 'filt_counts_*.h5' \
        --qc_pattern 'qc_metrics_*.csv.gz' \
        --integration_csv ${integration_csv} \
        --annotation_csv ${annotationArg} \
        --genome_gtf ${genome_gtf} \
        --out_all scqcflow_all_cells.h5ad \
        --out_clean scqcflow_clean_cells.h5ad
    """
}

process EXPORT_SEURAT {
    label     "process_high"
    tag       "export_seurat"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/export", mode: params.publish_mode_nonreport, overwrite: true

    input:
    path h5_files
    path qc_csvs
    path integration_csv
    path annotation_labels_csv
    path genome_gtf
    path script
    path utils_r

    output:
    path "scqcflow_all_cells.rds", emit: all_cells
    path "scqcflow_clean_cells.rds", emit: clean_cells

    script:
    def annotationArg = annotation_labels_csv.name == 'NO_FILE' ? 'NO_FILE' : annotation_labels_csv.name
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --h5_pattern 'filt_counts_*.h5' \
        --qc_pattern 'qc_metrics_*.csv.gz' \
        --integration_csv ${integration_csv} \
        --annotation_csv ${annotationArg} \
        --genome_gtf ${genome_gtf} \
        --utils_r ${utils_r} \
        --out_all scqcflow_all_cells.rds \
        --out_clean scqcflow_clean_cells.rds
    """
}