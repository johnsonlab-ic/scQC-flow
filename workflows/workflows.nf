/*
 * scQC-flow Workflow Definitions
 * 
 * This file contains the main workflow logic for both standard single-cell
 * and multiome data processing.
 */

// Import standard modules
include { CELLBENDER; CELLBENDER_GPU; CELLBENDER_H5_CONVERT; CELLBENDER_FROM_H5 } from '../modules/cellbender/cellbender'
include { GENERATE_ATAC_REPORT } from '../modules/reports/reports'
include { CREATE_SEURAT } from '../modules/seurat/seurat'
include { DROPLETQC } from '../modules/dropletqc/dropletqc'
include { SCDBL } from '../modules/scdbl/scdbl'
include { CELLBENDER_COMPARISON; CELLBENDER_COMPARISON_STATS_ONLY } from '../modules/cellbender_comparison/cellbender_comparison'
include { MAP_CELLRANGER } from '../modules/mapping/mapping'
include { SIMPLEAF_INDEX; SIMPLEAF_PREP_WHITELIST; SIMPLEAF_QUANT } from '../modules/alevinfry/alevinfry'
include { BARCODE_ESTIMATION; BARCODE_REPORT } from '../modules/barcode_estimation/barcode_estimation'
include { DECONTX; DECONTX_REPORT } from '../modules/decontx/decontx'

// =============================================================================
// RNA MAPPING WORKFLOW
// =============================================================================
workflow RNA_MAPPING_WORKFLOW {
    take:
        sampleChannelFastq      // tuple(sampleId, sampleName, fastqPath)
        mapper                  // string
        transcriptome           // string: required path to transcriptome reference
        cellrangerPath          // string: required path to cellranger installation

    main:
        if (mapper != 'cellranger') {
            error "Unsupported mapper: ${mapper}. Only 'cellranger' is currently supported."
        }

        if (!transcriptome) {
            error "--transcriptome is required for mapping with Cell Ranger."
        }

        if (!cellrangerPath) {
            error "--cellrangerPath is required for mapping with Cell Ranger."
        }

        map_input_ch = sampleChannelFastq.map { sampleId, sampleName, fastqPath ->
            tuple(sampleId, sampleName, file(fastqPath), cellrangerPath, transcriptome)
        }

        mapping_results = MAP_CELLRANGER(map_input_ch)
        sample_channel = mapping_results.mapped_dirs

    emit:
        sample_channel = sample_channel
}

// =============================================================================
// ALEVIN-FRY MAPPING WORKFLOW
// =============================================================================
workflow ALEVINFRY_MAPPING_WORKFLOW {
    take:
        sampleChannelFastq        // tuple(sampleId, sampleName, fastqPath)
        cellrangerPath            // string: path to Cell Ranger install/bin for whitelist extraction
        alevinfry_index           // string or null: path to pre-built index dir
        alevinfry_t2g             // string or null: path to t2g file
        alevinfry_whitelist       // string or null: optional whitelist override
        alevinfry_chemistry       // string: chemistry (e.g. 10xv3)
        alevinfry_orientation     // string: orientation (fw, rc, both)
        alevinfry_fasta           // string or null: FASTA for building index
        alevinfry_gtf             // string or null: GTF for building index

    main:
        // Map explicit chemistry names to simpleaf chemistry argument.
        def chemistryMap = [
            '3v2': '10xv2',
            '5v1': '10xv2',
            '5v2': '10xv2',
            '3v3': '10xv3',
            '5v3': '10xv3',
            '3v4': '10xv4-3p',
            '3LT': '10xv3LT',
            'multiome': '10xmultiome'
        ]
        af_chemistry = chemistryMap[alevinfry_chemistry]

        // Derive default orientation from explicit chemistry if requested.
        // 5' chemistries are reverse-complement; others default to forward.
        af_orientation = alevinfry_orientation == 'auto'
            ? (['5v1', '5v2', '5v3'].contains(alevinfry_chemistry) ? 'rc' : 'fw')
            : alevinfry_orientation

        // Resolve index: build if not provided, otherwise use pre-built
        if (alevinfry_index) {
            index_ch = channel.value(file("${alevinfry_index}/index"))
            t2g_ch   = channel.value(file(alevinfry_t2g))
        } else {
            SIMPLEAF_INDEX(alevinfry_fasta, alevinfry_gtf)
            index_ch = SIMPLEAF_INDEX.out.index_dir.map { it -> file("${it}/index") }
            t2g_ch   = SIMPLEAF_INDEX.out.index_dir.map { it -> file("${it}/ref/t2g_3col.tsv") }
        }

        if (alevinfry_whitelist) {
            whitelist_ch = channel.value(file(alevinfry_whitelist))
        } else {
            SIMPLEAF_PREP_WHITELIST(cellrangerPath, alevinfry_chemistry)
            whitelist_ch = SIMPLEAF_PREP_WHITELIST.out.whitelist
        }

        // Build per-sample quant input channel
        quant_input_ch = sampleChannelFastq
            .combine(index_ch)
            .combine(t2g_ch)
            .combine(whitelist_ch)
            .map { sampleId, sampleName, fastqPath, idx, t2g, wl ->
                tuple(sampleId, sampleName, file(fastqPath), idx, t2g, wl, af_chemistry, af_orientation)
            }

        quant_results = SIMPLEAF_QUANT(quant_input_ch)

        // Barcode estimation: H5 conversion + nuclear_fraction + knee params.
        // Mirrors scprocess save_alevin_to_h5 — this is the last step of mapping,
        // not the first step of QC.
        barcode_script     = channel.value(file("${projectDir}/modules/barcode_estimation/run_barcode_estimation.R"))
        barcode_report_qmd = channel.value(file("${projectDir}/modules/barcode_estimation/barcode_estimation.qmd"))
        BARCODE_ESTIMATION(quant_results.quant_dirs, barcode_script)
        BARCODE_REPORT(
            BARCODE_ESTIMATION.out.knee_data.map { _id, csv -> csv }.collect(),
            BARCODE_ESTIMATION.out.nuclear_fraction.map { _id, csv -> csv }.collect(),
            barcode_report_qmd
        )

    emit:
        quant_dirs       = quant_results.quant_dirs
        h5_files         = BARCODE_ESTIMATION.out.h5_files
        knee_data        = BARCODE_ESTIMATION.out.knee_data
        knee_params      = BARCODE_ESTIMATION.out.knee_params
        nuclear_fraction = BARCODE_ESTIMATION.out.nuclear_fraction
}

// Import multiome-specific modules
include { SCDBL_MULTIOME } from '../modules/multiome/scdbl_multiome'
include { CREATE_SEURAT_MULTIOME } from '../modules/multiome/seurat_multiome'
include { EXTRACT_MODALITIES } from '../modules/multiome/extract_modalities'
include { CELLBENDER_MULTIOME; CELLBENDER_MULTIOME_GPU } from '../modules/multiome/cellbender_multiome'

// =============================================================================
// STANDARD SINGLE-CELL WORKFLOW
// =============================================================================
workflow STANDARD_WORKFLOW {
    take:
        sampleChannelBase       // tuple(sampleName, mappingDir)
        dropletqc_script_path
        scdbl_script_path
        seurat_script_path
        cellbender              // boolean
        gpu                     // boolean
        max_mito                // double
        min_nuclear             // double
        metadata                // string or null

    main:
        if (cellbender) {
            log.info "Running CellBender workflow for all samples"
            
            // Run CellBender first
            if (gpu) {
                log.info "GPU acceleration enabled for CellBender"
                cellbender_results = CELLBENDER_GPU(sampleChannelBase)
            } else {
                cellbender_results = CELLBENDER(sampleChannelBase)
            }
            
            // Run H5 conversion after CellBender
            cellbender_h5_results = CELLBENDER_H5_CONVERT(cellbender_results.cellbender_output)

            // Prepare DropletQC inputs: BAM file + BAM index + CellBender barcodes
            dropletqc_input_ch = sampleChannelBase
                .join(cellbender_results.cellbender_output)
                .map { sampleName, mappingDir, cellbenderOutput -> 
                    def bamFile = file("${mappingDir}/outs/possorted_genome_bam.bam")
                    def bamIndex = file("${mappingDir}/outs/possorted_genome_bam.bam.bai")
                    def barcodesFile = file("${cellbenderOutput}/cellbender_out_cell_barcodes.csv")
                    tuple(sampleName, bamFile, bamIndex, barcodesFile, dropletqc_script_path)
                }

            // Prepare scDbl inputs: CellBender H5 file
            scdbl_input_ch = cellbender_h5_results.seurat_h5
                .map { sampleName, h5File -> tuple(sampleName, h5File, scdbl_script_path) }

            // Run DropletQC and scDbl with CellBender outputs
            dropletqc_results = DROPLETQC(dropletqc_input_ch)
            scdbl_results = SCDBL(scdbl_input_ch)

            // Run CellBender Comparison to analyze droplet calling differences
            cellbender_comparison_input_ch = sampleChannelBase
                .join(cellbender_h5_results.seurat_h5)
                .map { sampleName, mappingDir, h5File -> 
                    tuple(sampleName, mappingDir, h5File)
                }
            cellbender_comparison_results = CELLBENDER_COMPARISON(cellbender_comparison_input_ch).comparison_results

            // Create Seurat objects using CellBender H5 and updated QC metrics
            seurat_input_ch = sampleChannelBase
                .join(dropletqc_results.metrics)
                .join(scdbl_results.metrics)
                .join(cellbender_h5_results.seurat_h5)
                .map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) }

            seurat_input_with_script = seurat_input_ch.map { sampleName, mappingDir, dropletqc, scdbl, h5_path -> 
                tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, max_mito, min_nuclear, metadata, h5_path) 
            }
            seurat_results = CREATE_SEURAT(seurat_input_with_script)
            
        } else {
            log.info "Running standard workflow without CellBender"
            
            // Prepare DropletQC inputs: BAM file + BAM index + Cell Ranger barcodes
            dropletqc_input_ch = sampleChannelBase.map { sampleName, mappingDir -> 
                def bamFile = file("${mappingDir}/outs/possorted_genome_bam.bam")
                def bamIndex = file("${mappingDir}/outs/possorted_genome_bam.bam.bai")
                def barcodesFile = file("${mappingDir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
                tuple(sampleName, bamFile, bamIndex, barcodesFile, dropletqc_script_path)
            }

            // Prepare scDbl inputs: Cell Ranger H5 file
            scdbl_input_ch = sampleChannelBase.map { sampleName, mappingDir -> 
                def h5File = file("${mappingDir}/outs/filtered_feature_bc_matrix.h5")
                tuple(sampleName, h5File, scdbl_script_path)
            }

            // Run DropletQC and scDbl with Cell Ranger outputs
            dropletqc_results = DROPLETQC(dropletqc_input_ch)
            scdbl_results = SCDBL(scdbl_input_ch)

            // Run CellBender Comparison stats (Cell Ranger only) to analyze droplet calling
            cellbender_comparison_results = CELLBENDER_COMPARISON_STATS_ONLY(sampleChannelBase).comparison_results

            // Use default 10X H5 if CellBender is not run
            h5_path_ch = sampleChannelBase.map { sampleName, mappingDir ->
                def default_h5 = file("${mappingDir}/outs/filtered_feature_bc_matrix.h5")
                tuple(sampleName, default_h5)
            }
            
            seurat_input_ch = sampleChannelBase
                .join(dropletqc_results.metrics)
                .join(scdbl_results.metrics)
                .join(h5_path_ch)
                .map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) }

            seurat_input_with_script = seurat_input_ch.map { sampleName, mappingDir, dropletqc, scdbl, h5_path -> 
                tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, max_mito, min_nuclear, metadata, h5_path) 
            }
            seurat_results = CREATE_SEURAT(seurat_input_with_script)
        }

    emit:
        seurat_results = seurat_results
        cellbender_comparison_results = cellbender_comparison_results
}

// =============================================================================
// MULTIOME WORKFLOW
// =============================================================================
workflow MULTIOME_WORKFLOW {
    take:
        sampleChannelBase       // tuple(sampleName, mappingDir)
        dropletqc_script_path
        scdbl_script_path       // multiome version: run_scdbl_multiome.R
        seurat_script_path      // multiome version: make_seurat_multiome.R
        extract_modalities_script_path // R script to extract Gene Expression and ATAC from multiome H5
        cellbender              // boolean
        gpu                     // boolean
        max_mito                // double
        min_nuclear             // double
        metadata                // string or null

    main:
        if (cellbender) {
            log.info "Running CellBender workflow for multiome samples"
            
            // First, extract Gene Expression and ATAC modalities from multiome raw H5
            // GEX H5 is used for CellBender, ATAC H5 is passed to Seurat for ChromatinAssay
            extract_input_ch = sampleChannelBase.map { sampleName, mappingDir ->
                tuple(sampleName, mappingDir, extract_modalities_script_path)
            }
            extract_results = EXTRACT_MODALITIES(extract_input_ch)
            
            // Run CellBender on the extracted Gene Expression H5
            if (gpu) {
                log.info "GPU acceleration enabled for CellBender"
                cellbender_results = CELLBENDER_MULTIOME_GPU(extract_results.gex_h5)
            } else {
                cellbender_results = CELLBENDER_MULTIOME(extract_results.gex_h5)
            }
            
            // Run H5 conversion after CellBender
            cellbender_h5_results = CELLBENDER_H5_CONVERT(cellbender_results.cellbender_output)

            // Prepare DropletQC inputs: GEX BAM file + BAM index + CellBender barcodes
            // Note: Multiome uses gex_possorted_bam.bam instead of possorted_genome_bam.bam
            dropletqc_input_ch = sampleChannelBase
                .join(cellbender_results.cellbender_output)
                .map { sampleName, mappingDir, cellbenderOutput -> 
                    def bamFile = file("${mappingDir}/outs/gex_possorted_bam.bam")
                    def bamIndex = file("${mappingDir}/outs/gex_possorted_bam.bam.bai")
                    def barcodesFile = file("${cellbenderOutput}/cellbender_out_cell_barcodes.csv")
                    tuple(sampleName, bamFile, bamIndex, barcodesFile, dropletqc_script_path)
                }

            // Prepare scDbl inputs: CellBender H5 file (uses multiome R script)
            scdbl_input_ch = cellbender_h5_results.seurat_h5
                .map { sampleName, h5File -> tuple(sampleName, h5File, scdbl_script_path) }

            // Run DropletQC (same process, different BAM) and scDbl with multiome module
            dropletqc_results = DROPLETQC(dropletqc_input_ch)
            scdbl_results = SCDBL_MULTIOME(scdbl_input_ch)

            // Run CellBender Comparison to analyze droplet calling differences (multiome + cellbender)
            cellbender_comparison_input_ch = sampleChannelBase
                .join(cellbender_h5_results.seurat_h5)
                .map { sampleName, mappingDir, h5File -> 
                    tuple(sampleName, mappingDir, h5File)
                }
            cellbender_comparison_results = CELLBENDER_COMPARISON(cellbender_comparison_input_ch).comparison_results

            // Create Seurat objects using CellBender H5, ATAC H5, and multiome module
            // Join: sampleChannelBase + dropletqc + scdbl + cellbender_h5 + atac_h5
            seurat_input_ch = sampleChannelBase
                .join(dropletqc_results.metrics)
                .join(scdbl_results.metrics)
                .join(cellbender_h5_results.seurat_h5)
                .join(extract_results.atac_h5)
                // Structure: [sampleName, mappingDir, dropletqc_metrics, scdbl_metrics, gex_h5, atac_h5]
                .map { it -> tuple(it[0], it[1], it[2], it[3], it[4], it[5]) }

            seurat_input_with_script = seurat_input_ch.map { sampleName, mappingDir, dropletqc, scdbl, h5_path, atac_h5_path -> 
                tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, max_mito, min_nuclear, metadata, h5_path, atac_h5_path) 
            }
            seurat_multiome_output = CREATE_SEURAT_MULTIOME(seurat_input_with_script)
            seurat_results = seurat_multiome_output.seurat_rds
            atac_files = seurat_multiome_output.atac_files
            
        } else {
            log.info "Running multiome workflow without CellBender"
            
            // For non-CellBender multiome, we still need to extract the ATAC modality
            // from the raw H5 for Seurat ChromatinAssay creation
            extract_input_ch = sampleChannelBase.map { sampleName, mappingDir ->
                tuple(sampleName, mappingDir, extract_modalities_script_path)
            }
            extract_results = EXTRACT_MODALITIES(extract_input_ch)
            
            // Prepare DropletQC inputs: GEX BAM file + BAM index + Cell Ranger barcodes
            // Note: Multiome uses gex_possorted_bam.bam instead of possorted_genome_bam.bam
            dropletqc_input_ch = sampleChannelBase.map { sampleName, mappingDir -> 
                def bamFile = file("${mappingDir}/outs/gex_possorted_bam.bam")
                def bamIndex = file("${mappingDir}/outs/gex_possorted_bam.bam.bai")
                def barcodesFile = file("${mappingDir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
                tuple(sampleName, bamFile, bamIndex, barcodesFile, dropletqc_script_path)
            }

            // Prepare scDbl inputs: Cell Ranger H5 file (multiome version extracts Gene Expression)
            scdbl_input_ch = sampleChannelBase.map { sampleName, mappingDir -> 
                def h5File = file("${mappingDir}/outs/filtered_feature_bc_matrix.h5")
                tuple(sampleName, h5File, scdbl_script_path)
            }

            // Run DropletQC (same process) and scDbl with multiome module
            dropletqc_results = DROPLETQC(dropletqc_input_ch)
            scdbl_results = SCDBL_MULTIOME(scdbl_input_ch)

            // Run CellBender Comparison stats (Cell Ranger only) for multiome
            cellbender_comparison_results = CELLBENDER_COMPARISON_STATS_ONLY(sampleChannelBase).comparison_results

            // Use default 10X H5 if CellBender is not run
            h5_path_ch = sampleChannelBase.map { sampleName, mappingDir ->
                def default_h5 = file("${mappingDir}/outs/filtered_feature_bc_matrix.h5")
                tuple(sampleName, default_h5)
            }
            
            // Join: sampleChannelBase + dropletqc + scdbl + gex_h5 + atac_h5
            seurat_input_ch = sampleChannelBase
                .join(dropletqc_results.metrics)
                .join(scdbl_results.metrics)
                .join(h5_path_ch)
                .join(extract_results.atac_h5)
                // Structure: [sampleName, mappingDir, dropletqc_metrics, scdbl_metrics, gex_h5, atac_h5]
                .map { it -> tuple(it[0], it[1], it[2], it[3], it[4], it[5]) }

            seurat_input_with_script = seurat_input_ch.map { sampleName, mappingDir, dropletqc, scdbl, h5_path, atac_h5_path -> 
                tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, max_mito, min_nuclear, metadata, h5_path, atac_h5_path) 
            }
            seurat_multiome_output = CREATE_SEURAT_MULTIOME(seurat_input_with_script)
            seurat_results = seurat_multiome_output.seurat_rds
            atac_files = seurat_multiome_output.atac_files
        }

    emit:
        seurat_results = seurat_results
        atac_files = atac_files
        sample_channel = sampleChannelBase
        cellbender_comparison_results = cellbender_comparison_results
}

// =============================================================================
// ALEVIN-FRY QC WORKFLOW
// =============================================================================
// Runs QC without BAM files. Barcode estimation (H5 + nuclear_fraction + knee
// params) is performed as the last step of ALEVINFRY_MAPPING_WORKFLOW (mirroring
// scprocess save_alevin_to_h5), so this workflow receives those outputs directly.
// Ambient RNA removal is switchable: params.ambient_method = 'decontx' | 'cellbender'
workflow ALEVINFRY_QC_WORKFLOW {
    take:
        h5_files            // tuple(sampleId, h5)  — raw count matrix from BARCODE_ESTIMATION
        knee_data           // tuple(sampleId, csv) — knee plot data from BARCODE_ESTIMATION
        knee_params         // tuple(sampleId, txt) — CellBender/decontX params from BARCODE_ESTIMATION
        nuclear_fraction    // tuple(sampleId, csv) — nuclear fraction from BARCODE_ESTIMATION
        quant_dirs          // tuple(sampleId, dir) — quant dir; passed as mappingDir to Seurat
        scdbl_script_path
        seurat_script_path
        ambient_method      // string: 'decontx' or 'cellbender'
        gpu                 // boolean: enable CellBender GPU (only used if ambient_method == 'cellbender')
        max_mito            // double
        min_nuclear         // double
        max_nuclear         // double
        metadata            // string or null

    main:
        // ---------------------------------------------------------------
        // 1. Ambient RNA correction (decontX or CellBender)
        //    h5_files and knee outputs come from BARCODE_ESTIMATION,
        //    which ran at the end of ALEVINFRY_MAPPING_WORKFLOW.
        // ---------------------------------------------------------------
        if (ambient_method == 'decontx') {
            log.info "Ambient correction: decontX (barcodeRanks cell calling, CPU)"
            decontx_script     = channel.value(file("${projectDir}/modules/decontx/run_decontx.R"))
            decontx_report_qmd = channel.value(file("${projectDir}/modules/decontx/decontx_report.qmd"))
            decontx_input_ch = h5_files
                .join(knee_data)
                .map { sampleId, h5File, kneeCsv -> tuple(sampleId, h5File, kneeCsv) }
            DECONTX(decontx_input_ch, decontx_script)
            filtered_h5_ch = DECONTX.out.filtered_h5
            DECONTX_REPORT(
                DECONTX.out.usa_metrics.map { _id, csv -> csv }.collect(),
                DECONTX.out.barcodes.map    { _id, csv -> csv }.collect(),
                knee_data.map               { _id, csv -> csv }.collect(),
                decontx_report_qmd
            )

        } else if (ambient_method == 'cellbender') {
            log.info "Ambient correction: CellBender (${gpu ? 'GPU' : 'CPU'})"
            knee_vals_ch = knee_params
                .map { sampleId, paramsFile ->
                    def parts = paramsFile.text.trim().split(',')
                    tuple(sampleId, parts[0] as Integer, parts[1] as Integer, parts[2] as Integer)
                }
            cellbender_input_ch = h5_files
                .join(knee_vals_ch)
                .map { sampleId, h5File, ec, td, lct -> tuple(sampleId, h5File, ec, td, lct) }
            CELLBENDER_FROM_H5(cellbender_input_ch)
            filtered_h5_ch = CELLBENDER_FROM_H5.out.filtered_h5

        } else {
            error "Unknown ambient_method: ${ambient_method}. Choose 'decontx' or 'cellbender'."
        }

        // ---------------------------------------------------------------
        // 2. Doublet detection on the ambient-corrected H5
        // ---------------------------------------------------------------
        scdbl_input_ch = filtered_h5_ch
            .map { sampleId, h5File ->
                tuple(sampleId, h5File, scdbl_script_path)
            }
        scdbl_results = SCDBL(scdbl_input_ch)

        // ---------------------------------------------------------------
        // 3. Seurat object creation
        //    nuclear_fraction substitutes for dropletqc_metrics.
        //    quant_dirs is passed as mappingDir (unused when h5_path is provided).
        // ---------------------------------------------------------------
        seurat_input_ch = filtered_h5_ch
            .join(nuclear_fraction)
            .join(scdbl_results.metrics)
            .join(quant_dirs)
            .map { sampleId, h5File, nf_csv, scdbl_csv, quantDir ->
                tuple(sampleId, quantDir, nf_csv, scdbl_csv,
                      seurat_script_path, max_mito, min_nuclear, metadata, h5File)
            }
        seurat_results = CREATE_SEURAT(seurat_input_ch)

    emit:
        seurat_results   = seurat_results
        nuclear_fraction = nuclear_fraction
        h5_files         = h5_files
        filtered_h5      = filtered_h5_ch
}