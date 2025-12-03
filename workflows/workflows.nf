/*
 * scQC-flow Workflow Definitions
 * 
 * This file contains the main workflow logic for both standard single-cell
 * and multiome data processing.
 */

// Import modules
include { CELLBENDER; CELLBENDER_GPU; CELLBENDER_H5_CONVERT } from './modules/cellbender/cellbender'
include { GENERATE_REPORTS; COMBINE_REPORTS; GENERATE_COMBINED_REPORT } from './modules/reports/reports'
include { CREATE_SEURAT } from './modules/seurat/seurat'
include { DROPLETQC } from './modules/dropletqc/dropletqc'
include { SCDBL } from './modules/scdbl/scdbl'

// =============================================================================
// STANDARD SINGLE-CELL WORKFLOW
// =============================================================================
workflow STANDARD_WORKFLOW {
    take:
        sampleChannelBase       // tuple(sampleName, mappingDir)
        dropletqc_script_path
        scdbl_script_path
        seurat_script_path

    main:
        if (params.cellbender) {
            log.info "Running CellBender workflow for all samples"
            
            // Run CellBender first
            if (params.gpu) {
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

            // Create Seurat objects using CellBender H5 and updated QC metrics
            seurat_input_ch = sampleChannelBase
                .join(dropletqc_results.metrics)
                .join(scdbl_results.metrics)
                .join(cellbender_h5_results.seurat_h5)
                .map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) }

            seurat_input_with_script = seurat_input_ch.map { sampleName, mappingDir, dropletqc, scdbl, h5_path -> 
                tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, params.max_mito, params.min_nuclear, params.metadata, h5_path) 
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
                tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, params.max_mito, params.min_nuclear, params.metadata, h5_path) 
            }
            seurat_results = CREATE_SEURAT(seurat_input_with_script)
        }

    emit:
        seurat_results = seurat_results
}

// =============================================================================
// MULTIOME WORKFLOW
// =============================================================================
workflow MULTIOME_WORKFLOW {
    take:
        sampleChannelBase       // tuple(sampleName, mappingDir)
        dropletqc_script_path
        scdbl_script_path
        seurat_script_path

    main:
        if (params.cellbender) {
            log.info "Running CellBender workflow for all samples"
            
            // Run CellBender first
            if (params.gpu) {
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

            // Create Seurat objects using CellBender H5 and updated QC metrics
            seurat_input_ch = sampleChannelBase
                .join(dropletqc_results.metrics)
                .join(scdbl_results.metrics)
                .join(cellbender_h5_results.seurat_h5)
                .map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) }

            seurat_input_with_script = seurat_input_ch.map { sampleName, mappingDir, dropletqc, scdbl, h5_path -> 
                tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, params.max_mito, params.min_nuclear, params.metadata, h5_path) 
            }
            seurat_results = CREATE_SEURAT(seurat_input_with_script)
            
        } else {
            log.info "Running standard multiome workflow without CellBender"
            
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
                tuple(sampleName, mappingDir, dropletqc, scdbl, seurat_script_path, params.max_mito, params.min_nuclear, params.metadata, h5_path) 
            }
            seurat_results = CREATE_SEURAT(seurat_input_with_script)
        }

    emit:
        seurat_results = seurat_results
}

// =============================================================================
// REPORTING WORKFLOW
// =============================================================================
workflow REPORTING {
    take:
        sampleChannelBase       // tuple(sampleName, mappingDir)
        seurat_results          // tuple(sampleName, pre_rds, post_rds)
        report_template_path
        combined_template_path
        book_template_path

    main:
        // Prepare GENERATE_REPORTS input channel
        report_input_ch = sampleChannelBase
            .join(seurat_results)
            .map { sampleName, mappingDir, pre_rds, post_rds -> 
                tuple(sampleName, mappingDir, pre_rds, post_rds, report_template_path, params.max_mito, params.min_nuclear) 
            }

        if (params.report) {
            reports_output = GENERATE_REPORTS(report_input_ch)

            // Optionally combine all reports into a single book
            if (params.book) {
                log.info "Combining reports into a Quarto book"

                all_html_reports = reports_output.html_report.map { sampleName, htmlFile -> htmlFile }.collect()
                all_qmd_sources = reports_output.qmd_source.map { sampleName, qmdFile -> qmdFile }.collect()

                combined_book = COMBINE_REPORTS(all_html_reports, all_qmd_sources, book_template_path)
            }
        }
}
