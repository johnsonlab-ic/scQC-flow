// alevin-fry mapping module
// Uses simpleaf for FASTQ -> count matrix quantification (alternative to Cell Ranger).
// NOTE: alevin-fry is a pseudo-aligner and does not produce BAM files.
//       Downstream QC modules (DropletQC, scDblFinder, Seurat) are NOT yet
//       compatible with alevin-fry output. Use --run_mode mapping only.

process SIMPLEAF_INDEX {
    label "process_simpleaf"
    tag "simpleaf_index"
    publishDir "${params.outputDir}/alevinfry_index", mode: 'copy', overwrite: true
    conda "bioconda::simpleaf"

    input:
    val(fasta)
    val(gtf)

    output:
    path("simpleaf_index"), emit: index_dir

    script:
    """
    set -euo pipefail

    export ALEVIN_FRY_HOME=\$(pwd)/af_home
    mkdir -p \$ALEVIN_FRY_HOME
    simpleaf set-paths

    if [ ! -f "${fasta}" ]; then
        echo "ERROR: FASTA file does not exist: ${fasta}"
        exit 1
    fi
    if [ ! -f "${gtf}" ]; then
        echo "ERROR: GTF file does not exist: ${gtf}"
        exit 1
    fi

    GTF_FOR_INDEX="${gtf}"
    if [[ "${gtf}" == *.gz ]]; then
        echo "Detected gzipped GTF; decompressing for simpleaf index build"
        gzip -cd "${gtf}" > reference.gtf
        GTF_FOR_INDEX="reference.gtf"
    fi

    echo "Building simpleaf index"
    echo "FASTA: ${fasta}"
    echo "GTF:   \$GTF_FOR_INDEX"

    simpleaf index \
        --fasta "${fasta}" \
        --gtf "\$GTF_FOR_INDEX" \
        --output simpleaf_index \
        --threads ${task.cpus}
    """
}

process SIMPLEAF_PREP_WHITELIST {
    label "process_low"
    tag "${chemistry}_whitelist"
    publishDir "${params.outputDir}/alevinfry_index", mode: 'copy', overwrite: true

    input:
    val(cellrangerPath)
    val(chemistry)

    output:
    path("chemistry_whitelist.txt"), emit: whitelist

    script:
    """
    set -euo pipefail

    # Resolve Cell Ranger home from either install dir or binary path.
    if [ -x "${cellrangerPath}/cellranger" ]; then
        CR_HOME="${cellrangerPath}"
    elif [ -x "${cellrangerPath}" ]; then
        CR_HOME=\$(cd "\$(dirname \"${cellrangerPath}\")/.." && pwd)
    else
        echo "ERROR: Cannot find Cell Ranger at ${cellrangerPath}"
        exit 1
    fi

    BARCODE_DIR="\$CR_HOME/lib/python/cellranger/barcodes"
    if [ ! -d "\$BARCODE_DIR" ]; then
        echo "ERROR: Cell Ranger barcode directory not found: \$BARCODE_DIR"
        exit 1
    fi

    # Map explicit chemistry name to Cell Ranger whitelist file.
    case "${chemistry}" in
        3v2|5v1|5v2)
            WL_FILE="737K-august-2016.txt"
            ;;
        3v3)
            WL_FILE="3M-february-2018_TRU.txt.gz"
            ;;
        5v3)
            WL_FILE="3M-5pgex-jan-2023.txt.gz"
            ;;
        3v4)
            WL_FILE="3M-3pgex-may-2023_TRU.txt.gz"
            ;;
        3LT)
            WL_FILE="9K-LT-march-2021.txt.gz"
            ;;
        multiome)
            WL_FILE="737K-arc-v1.txt.gz"
            ;;
        *)
            echo "ERROR: Unsupported --alevinfry_chemistry '${chemistry}' for whitelist extraction"
            echo "Supported: 3v2, 5v1, 5v2, 3v3, 5v3, 3v4, 3LT, multiome"
            exit 1
            ;;
    esac

    SRC="\$BARCODE_DIR/\$WL_FILE"
    if [ ! -f "\$SRC" ]; then
        echo "ERROR: Expected whitelist file not found: \$SRC"
        exit 1
    fi

    if [[ "\$SRC" == *.gz ]]; then
        gzip -cd "\$SRC" > chemistry_whitelist.txt
    else
        cp "\$SRC" chemistry_whitelist.txt
    fi

    echo "Prepared whitelist for ${chemistry}: \$SRC -> chemistry_whitelist.txt"
    """
}

process SIMPLEAF_QUANT {
    label "process_simpleaf"
    tag "$sampleId"
    publishDir "${params.outputDir}/alevinfry", mode: 'copy', overwrite: true
    conda "bioconda::simpleaf"

    input:
    tuple val(sampleId), val(sampleName), path(fastqPath), path(indexDir), path(t2gMap), path(whitelist), val(chemistry), val(orientation)

    output:
    tuple val(sampleId), path("${sampleId}_af_quant"), emit: quant_dirs

    script:
    """
    set -euo pipefail

    export ALEVIN_FRY_HOME=\$(pwd)/af_home
    mkdir -p \$ALEVIN_FRY_HOME
    simpleaf set-paths

        # Discover R1 and R2 FASTQ files.
        # In Nextflow work dirs, input directories are often staged as symlinks.
        # Use -L so find traverses symlinked directories reliably.
        R1_FILES=\$(find -L ${fastqPath} -name "*_R1_*" -type f | sort | tr '\\n' ',' | sed 's/,\$//')
        R2_FILES=\$(find -L ${fastqPath} -name "*_R2_*" -type f | sort | tr '\\n' ',' | sed 's/,\$//')

    if [ -z "\$R1_FILES" ] || [ -z "\$R2_FILES" ]; then
        echo "ERROR: Could not find R1/R2 FASTQ files in ${fastqPath}"
        echo "Contents:"
        ls -la ${fastqPath}/
        exit 1
    fi

    echo "Mapping sample ${sampleId} (${sampleName})"
    echo "FASTQ path: ${fastqPath}"
    echo "R1 files: \$R1_FILES"
    echo "R2 files: \$R2_FILES"
    echo "Index: ${indexDir}"
    echo "T2G map: ${t2gMap}"
    echo "Whitelist: ${whitelist}"
    echo "Chemistry: ${chemistry}"
    echo "Orientation: ${orientation}"

    INDEX_PATH="${indexDir}"

    # Compatibility shim for older pre-built simpleaf indices.
    # Newer piscem expects piscem_idx.ssi, while some older indices contain
    # piscem_idx.sshash instead. Create a local compatibility copy if needed.
    if [ ! -f "${indexDir}/piscem_idx.ssi" ] && [ -f "${indexDir}/piscem_idx.sshash" ]; then
        echo "Detected legacy index format (piscem_idx.sshash without piscem_idx.ssi)."
        echo "Creating compatibility index view in work directory."
        mkdir -p compat_index
        cp -a ${indexDir}/. compat_index/
        ln -sf piscem_idx.sshash compat_index/piscem_idx.ssi
        INDEX_PATH="compat_index"
    fi

    if [ ! -f "\$INDEX_PATH/piscem_idx.ssi" ]; then
        echo "ERROR: Could not find \$INDEX_PATH/piscem_idx.ssi"
        echo "Index directory contents:"
        ls -la \$INDEX_PATH/ || true
        echo "Hint: rebuild the index with the same simpleaf version used for quant."
        exit 1
    fi

    simpleaf quant \
        --reads1 \$R1_FILES \
        --reads2 \$R2_FILES \
        --threads ${task.cpus} \
        --index \$INDEX_PATH \
        --chemistry ${chemistry} \
        --resolution cr-like \
        --expected-ori ${orientation} \
        --t2g-map ${t2gMap} \
        --unfiltered-pl ${whitelist} \
        --min-reads 1 \
        --output ${sampleId}_af_quant
    """
}

process ALEVINFRY_TO_H5 {
    label "process_medium"
    tag "$sampleId"
    publishDir "${params.outputDir}/alevinfry", mode: 'copy', overwrite: true
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"

    input:
    tuple val(sampleId), path(quantDir)
    path(script)

    output:
    tuple val(sampleId), path("${sampleId}_counts.h5"), emit: h5_files

    script:
    """
    Rscript "${script}" \\
        "${sampleId}" \
        "${quantDir}" \
        "${sampleId}_counts.h5"
    """
}
