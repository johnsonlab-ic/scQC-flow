process EXPORT_SCANPY {
    label     "process_high"
    tag       "export_anndata"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/export", mode: params.publish_mode_nonreport, overwrite: true

    input:
    path h5_files
    path qc_csvs
    path integration_csv
    path annotation_metadata_csvs
    path genome_gtf
    path script

    output:
    path "anndata/sample_*_all.h5ad", emit: per_sample_all
    path "anndata/sample_*_clean.h5ad", emit: per_sample_clean
    path "anndata/combined_all.h5ad", optional: true, emit: combined_all
    path "anndata/combined_clean.h5ad", optional: true, emit: combined_clean

    script:
    def writeCombined = params.export_write_combined ?: true
    """
    set -euo pipefail
    export MPLCONFIGDIR="/tmp/.mplconfig"
    export NUMBA_CACHE_DIR="/tmp/.numba"
    export PYTHONUNBUFFERED=1

    python3 -u ${script} \
        --h5_pattern 'filt_counts_*.h5' \
        --qc_pattern 'qc_metrics_*.csv.gz' \
        --integration_csv ${integration_csv} \
        --annotation_pattern 'annotation_export_*.csv.gz' \
        --genome_gtf ${genome_gtf} \
        --out_dir . \
        --write_combined ${writeCombined}
    """
}
    
process EXPORT_CELL_METADATA {
    label     "process_low"
    tag       "export_cell_metadata"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/export", mode: params.publish_mode_nonreport, overwrite: true

    input:
    path qc_csvs
    path integration_csv
    path annotation_metadata_csvs
    path script
    path utils_r

    output:
    path "metadata/cell_metadata_all.csv.gz", emit: all
    path "metadata/cell_metadata_clean.csv.gz", optional: true, emit: clean

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --qc_pattern 'qc_metrics_*.csv.gz' \
        --integration_csv ${integration_csv} \
        --annotation_pattern 'annotation_export_*.csv.gz' \
        --utils_r ${utils_r} \
        --output_all_csv metadata/cell_metadata_all.csv.gz \
        --output_clean_csv metadata/cell_metadata_clean.csv.gz
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
    path annotation_metadata_csvs
    path genome_gtf
    path script
    path utils_r

    output:
    path "seurat/sample_*_all.rds", emit: per_sample_all
    path "seurat/sample_*_clean.rds", emit: per_sample_clean
    path "seurat/combined_all.rds", optional: true, emit: combined_all
    path "seurat/combined_clean.rds", optional: true, emit: combined_clean

    script:
    def writeCombined = params.export_write_combined ?: true
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --h5_pattern 'filt_counts_*.h5' \
        --qc_pattern 'qc_metrics_*.csv.gz' \
        --integration_csv ${integration_csv} \
        --annotation_pattern 'annotation_export_*.csv.gz' \
        --genome_gtf ${genome_gtf} \
        --utils_r ${utils_r} \
        --out_dir . \
        --write_combined ${writeCombined}
    """
}

process EXPORT_SCANPY_ZOOM {
    label     "process_high"
    tag       "export_anndata_zoom_${zoom_name}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/zoom", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${zoom_name}/exports/${filename}" }

    input:
    tuple val(zoom_name), val(spec_b64), path(zoom_qc_metrics), path(zoom_selection), path(zoom_integration_csv)
    path h5_files
    path annotation_metadata_csvs
    path genome_gtf
    path script

    output:
    path "anndata/sample_*_all.h5ad", emit: per_sample_all
    path "anndata/sample_*_clean.h5ad", optional: true, emit: per_sample_clean
    path "anndata/combined_all.h5ad", optional: true, emit: combined_all
    path "anndata/combined_clean.h5ad", optional: true, emit: combined_clean

    script:
    def writeCombined = params.export_write_combined ?: true
    """
    set -euo pipefail
    export MPLCONFIGDIR="/tmp/.mplconfig"
    export NUMBA_CACHE_DIR="/tmp/.numba"
    export PYTHONUNBUFFERED=1

    python3 -u ${script} \
        --h5_pattern 'filt_counts_*.h5' \
        --qc_pattern 'zoom_qc_metrics.csv.gz' \
        --integration_csv ${zoom_integration_csv} \
        --annotation_pattern 'annotation_export_*.csv.gz' \
        --genome_gtf ${genome_gtf} \
        --out_dir . \
        --write_combined ${writeCombined}
    """
}

process EXPORT_SEURAT_ZOOM {
    label     "process_high"
    tag       "export_seurat_zoom_${zoom_name}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/zoom", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${zoom_name}/exports/${filename}" }

    input:
    tuple val(zoom_name), val(spec_b64), path(zoom_qc_metrics), path(zoom_selection), path(zoom_integration_csv)
    path h5_files
    path annotation_metadata_csvs
    path genome_gtf
    path script
    path utils_r

    output:
    path "seurat/sample_*_all.rds", emit: per_sample_all
    path "seurat/sample_*_clean.rds", optional: true, emit: per_sample_clean
    path "seurat/combined_all.rds", optional: true, emit: combined_all
    path "seurat/combined_clean.rds", optional: true, emit: combined_clean

    script:
    def writeCombined = params.export_write_combined ?: true
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --h5_pattern 'filt_counts_*.h5' \
        --qc_pattern 'zoom_qc_metrics.csv.gz' \
        --integration_csv ${zoom_integration_csv} \
        --annotation_pattern 'annotation_export_*.csv.gz' \
        --genome_gtf ${genome_gtf} \
        --utils_r ${utils_r} \
        --out_dir . \
        --write_combined ${writeCombined}
    """
}

process EXPORT_CELL_METADATA_ZOOM {
    label     "process_low"
    tag       "export_cell_metadata_zoom_${zoom_name}"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/zoom", mode: params.publish_mode_nonreport, overwrite: true, saveAs: { filename -> "${zoom_name}/exports/${filename}" }

    input:
    tuple val(zoom_name), val(spec_b64), path(zoom_qc_metrics), path(zoom_selection), path(zoom_integration_csv)
    path annotation_metadata_csvs
    path script
    path utils_r

    output:
    path "metadata/cell_metadata_all.csv.gz", emit: all
    path "metadata/cell_metadata_clean.csv.gz", optional: true, emit: clean

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"

    Rscript ${script} \
        --qc_pattern 'zoom_qc_metrics.csv.gz' \
        --integration_csv ${zoom_integration_csv} \
        --annotation_pattern 'annotation_export_*.csv.gz' \
        --utils_r ${utils_r} \
        --output_all_csv metadata/cell_metadata_all.csv.gz \
        --output_clean_csv metadata/cell_metadata_clean.csv.gz
    """
}
