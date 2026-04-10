// BARCODE_ESTIMATION
//
// Loads an alevin-fry quant directory and produces four outputs that feed
// the rest of the alevin-fry QC workflow:
//
//   h5_files           — unfiltered total (S+U+A) count matrix (.h5, 10x v3)
//                        used as input for CellBender or decontX
//   nuclear_fraction   — per-barcode nuclear_fraction CSV (U/(S+U));
//                        equivalent to DropletQC output for Cell Ranger samples
//   knee_data          — per-barcode knee estimation data (rank, total, knee1/shin1/knee2/shin2,
//                        expected_cells, total_droplets_included, low_count_threshold,
//                        in_empty_plateau, spliced, unspliced); input for BARCODE_REPORT
//   knee_params        — single-line CSV of knee-derived CellBender/decontX
//                        parameters; unpacked into channel vals by the workflow
//                        (no YAML hand-off — values travel in the channel tuple)
//
// BARCODE_REPORT
//
// Collects knee_data CSVs from all samples, renders barcode_estimation.qmd,
// and publishes a single HTML report. Runs immediately after all BARCODE_ESTIMATION
// tasks complete (blocking — pipeline fails if report fails).

process BARCODE_ESTIMATION {
    label "process_high"
    tag  "$sampleId"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/barcode_estimation", mode: 'copy', overwrite: true,
        saveAs: { filename ->
            // Publish H5, nuclear_fraction, and knee_data; knee_params is
            // an internal artefact used only to bootstrap channel values.
            filename.endsWith("_knee_params.txt") ? null : filename
        }

    input:
    tuple val(sampleId), path(quantDir)
    path(script)

    output:
    tuple val(sampleId), path("${sampleId}_counts.h5"),             emit: h5_files
    tuple val(sampleId), path("${sampleId}_nuclear_fraction.csv"),  emit: nuclear_fraction
    tuple val(sampleId), path("${sampleId}_knee_data.csv"),         emit: knee_data
    tuple val(sampleId), path("${sampleId}_knee_params.txt"),       emit: knee_params

    script:
    """
    Rscript "${script}" \\
        "${sampleId}" \\
        "${quantDir}" \\
        "${sampleId}_counts.h5" \\
        "${sampleId}_nuclear_fraction.csv" \\
        "${sampleId}_knee_data.csv" \\
        "${sampleId}_knee_params.txt"
    """
}

process BARCODE_REPORT {
    label "process_reports"
    tag  "barcode_estimation_report"
    container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/reports", mode: 'copy', overwrite: true

    input:
    path(knee_data_csvs)   // collected list of *_knee_data.csv files (all samples)
    path(nf_csvs)          // collected list of *_nuclear_fraction.csv files (all samples)
    path(report_qmd)       // barcode_estimation.qmd (from projectDir)

    output:
    path("barcode_estimation_report.html"), emit: html_report

    script:
    """
    export HOME="\$PWD"
    cp "${report_qmd}" barcode_estimation_report.qmd
    quarto render barcode_estimation_report.qmd \\
        --output barcode_estimation_report.html
    """
}
