// mapping.nf
//
// All processes for the alevin-fry MAPPING workflow:
//   SIMPLEAF_INDEX     — build index from FASTA + GTF
//   SIMPLEAF_QUANT     — run simpleaf quant (resolves whitelist from chemistry)
//   BARCODE_ESTIMATION — knee detection, H5 creation, ambient params
//   BARCODE_ESTIMATION_V2 — standalone barcode caller + audit outputs

// ---------------------------------------------------------------------------
// Index building (always runs — takes genome FASTA + GTF)
// ---------------------------------------------------------------------------

process SIMPLEAF_INDEX {
    label  "process_simpleaf"
    tag    "simpleaf_index"
    publishDir "${params.outputDir}/simpleaf_index", mode: params.publish_mode_nonreport, overwrite: true
    conda  "bioconda::simpleaf"

    input:
    path fasta
    path gtf

    output:
    path "simpleaf_index/index", emit: index_dir
    path "simpleaf_index/index/t2g_3col.tsv", emit: t2g

    script:
    """
    set -euo pipefail

    export ALEVIN_FRY_HOME=\$(pwd)/af_home
    mkdir -p \$ALEVIN_FRY_HOME
    simpleaf set-paths

    GTF_FOR_INDEX="${gtf}"
    if [[ "${gtf}" == *.gz ]]; then
        gzip -cd "${gtf}" > reference.gtf
        GTF_FOR_INDEX="reference.gtf"
    fi

    simpleaf index \
        --fasta   "${fasta}" \
        --gtf     "\$GTF_FOR_INDEX" \
        --output  simpleaf_index \
        --threads ${task.cpus}
    """
}

// ---------------------------------------------------------------------------
// simpleaf quantification
//
// Resolves the whitelist from cellrangerPath + chemistry internally.
// No separate detection step needed.
// ---------------------------------------------------------------------------

process SIMPLEAF_QUANT {
    label  "process_simpleaf"
    tag    "$sampleId"
    conda  "bioconda::simpleaf"
    publishDir "${params.outputDir}/mapping/af_${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), val(sampleName), path(fastqPath)
    val   af_chemistry
    val   orientation
    val   whitelist_filename
    val   cellrangerPath
    path  index_dir
    path  t2g_map

    output:
    tuple val(sampleId), path("af_quant"), emit: quant_dirs
    path  "simpleaf_quant_log.json", emit: quant_log

    script:
    """
    set -euo pipefail

    export ALEVIN_FRY_HOME=\$(pwd)/af_home
    mkdir -p \$ALEVIN_FRY_HOME
    simpleaf set-paths

    # Resolve and decompress whitelist from Cell Ranger barcodes dir
    WL_SRC="${cellrangerPath}/lib/python/cellranger/barcodes/${whitelist_filename}"
    if [ ! -f "\$WL_SRC" ]; then
        echo "ERROR: Whitelist not found: \$WL_SRC" >&2
        exit 1
    fi
    if [[ "\$WL_SRC" == *.gz ]]; then
        gzip -cd "\$WL_SRC" > whitelist.txt
    else
        cp "\$WL_SRC" whitelist.txt
    fi

    # Discover R1 / R2 files
    R1_FILES=\$(find -L ${fastqPath} -name "*_R1_*" -type f | sort | paste -sd, -)
    R2_FILES=\$(find -L ${fastqPath} -name "*_R2_*" -type f | sort | paste -sd, -)

    if [ -z "\$R1_FILES" ] || [ -z "\$R2_FILES" ]; then
        echo "ERROR: Could not find R1/R2 FASTQ files in ${fastqPath}" >&2
        exit 1
    fi

    INDEX_PATH="${index_dir}"
    if [ ! -f "\$INDEX_PATH/piscem_idx.ssi" ] && [ -f "\$INDEX_PATH/piscem_idx.sshash" ]; then
        mkdir -p compat_index
        cp -a \$INDEX_PATH/. compat_index/
        ln -sf piscem_idx.sshash compat_index/piscem_idx.ssi
        INDEX_PATH="compat_index"
    fi

    simpleaf quant \
        --reads1        \$R1_FILES \
        --reads2        \$R2_FILES \
        --threads       ${task.cpus} \
        --index         \$INDEX_PATH \
        --chemistry     ${af_chemistry} \
        --resolution    cr-like \
        --expected-ori  ${orientation} \
        --t2g-map       ${t2g_map} \
        --unfiltered-pl whitelist.txt \
        --min-reads     1 \
        --output        af_quant

    cp af_quant/af_quant/quant.json simpleaf_quant_log.json || true
    """
}

// ---------------------------------------------------------------------------
// Barcode estimation: knee detection, H5, ambient params, knee data CSV
// ---------------------------------------------------------------------------

process BARCODE_ESTIMATION {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/mapping/af_${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(quant_dir)
    path  script

    output:
    tuple val(sampleId), path("af_counts_mat.h5"), emit: h5_files
    tuple val(sampleId), path("knee_plot_data_${sampleId}.csv"), emit: knee_data
    tuple val(sampleId), env('CB_TOTAL_DROPLETS_INCLUDED'), env('CB_EXPECTED_CELLS'), env('CB_LOW_COUNT_THRESHOLD'), env('KNEE1'), env('SHIN1'), env('KNEE2'), env('SHIN2'), emit: ambient_params

    script:
    """
    set -euo pipefail

    Rscript ${script} \
        "${sampleId}" \
        "${quant_dir}/af_quant" \
        "af_counts_mat.h5" \
        "knee_plot_data_${sampleId}.csv" \
        "ambient_params_${sampleId}.env"

    source "ambient_params_${sampleId}.env"
    export CB_TOTAL_DROPLETS_INCLUDED
    export CB_EXPECTED_CELLS
    export CB_LOW_COUNT_THRESHOLD
    export KNEE1
    export SHIN1
    export KNEE2
    export SHIN2
    """
}

// ---------------------------------------------------------------------------
// Barcode estimation v2: standalone caller with explicit decision zones.
// ---------------------------------------------------------------------------

process BARCODE_ESTIMATION_V2 {
    label     "process_high"
    tag       "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/mapping_v2/af_${sampleId}", mode: params.publish_mode_nonreport, overwrite: true

    input:
    tuple val(sampleId), path(quant_dir)
    path  script

    output:
    tuple val(sampleId), path("af_counts_mat_v2.h5"), emit: h5_files
    tuple val(sampleId), path("barcode_audit_v2_${sampleId}.csv.gz"), emit: audit
    tuple val(sampleId), path("barcode_summary_v2_${sampleId}.csv"), emit: summaries
    tuple val(sampleId), path("cell_barcodes_v2_${sampleId}.csv"), emit: barcodes
    tuple val(sampleId), env('CB_EXPECTED_CELLS_V2'), env('CB_LOW_COUNT_THRESHOLD_V2'), env('KNEE1_V2'), env('SHIN1_V2'), env('KNEE2_V2'), env('SHIN2_V2'), env('ED_FDR_V2'), emit: ambient_params

    script:
    """
    set -euo pipefail

    Rscript ${script} \
        --sample_id          "${sampleId}" \
        --quant_dir          "${quant_dir}/af_quant" \
        --h5_out             "af_counts_mat_v2.h5" \
        --audit_out          "barcode_audit_v2_${sampleId}.csv.gz" \
        --summary_out        "barcode_summary_v2_${sampleId}.csv" \
        --barcodes_out       "cell_barcodes_v2_${sampleId}.csv" \
        --ambient_env_out    "ambient_params_v2_${sampleId}.env" \
        --min_umis_empty     ${params.barcode_v2_min_umis_empty} \
        --niters             ${params.barcode_v2_ed_niters} \
        --ed_fdr             ${params.barcode_v2_ed_fdr} \
        --splice_context     '${params.barcode_v2_splice_context}' \
        --low_count_strategy '${params.barcode_v2_low_count_strategy}'

    source "ambient_params_v2_${sampleId}.env"
    export CB_EXPECTED_CELLS_V2
    export CB_LOW_COUNT_THRESHOLD_V2
    export KNEE1_V2
    export SHIN1_V2
    export KNEE2_V2
    export SHIN2_V2
    export ED_FDR_V2
    """
}

