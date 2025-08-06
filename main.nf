nextflow.enable.dsl=2

// Default parameters
params.mapping_dirs = "${projectDir}/personal/mapping_dirs.csv"
params.outputDir = "results"
params.gpu = false
params.help = false

// Help message
def helpMessage() {
    log.info """
    ===================================
    scQC-flow
    ===================================
    
    Usage:
    nextflow run main.nf --mapping_dirs <mapping_dirs.csv> --outputDir <output_directory>
    
    Required arguments:
      --mapping_dirs        Path to mapping directories CSV file
      --outputDir           Output directory for results
      
    Optional arguments:
      --gpu                 Use GPU acceleration for cellbender
      --help                Show this help message
      
    Example:
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results --gpu
    """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Process for cellbender
process cellbender {
    label "process_high_memory"
    tag { sampleName }
    container "us.gcr.io/broad-dsde-methods/cellbender:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir)

    output:
    path "${sampleName}_cellbender_output"

    script:
    """
    echo "Running cellbender for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"

    mkdir -p ${sampleName}_cellbender_output

    cellbender remove-background \
                 --input ${mappingDir}/outs/raw_feature_bc_matrix.h5 \
                 --output ${sampleName}_cellbender_output/cellbender_out.h5 

    echo "Cellbender processing completed" > ${sampleName}_cellbender_output/summary.txt
    echo "Cellbender completed for ${sampleName}"
    """
}

process cellbender_gpu {
    label "process_gpu"
    tag { sampleName }
    container "us.gcr.io/broad-dsde-methods/cellbender:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir)

    output:
    path "${sampleName}_cellbender_output"

    script:
    """
    echo "Running cellbender for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"

    mkdir -p ${sampleName}_cellbender_output

    cellbender remove-background \
                 --input ${mappingDir}/outs/raw_feature_bc_matrix.h5 \
                 --output ${sampleName}_cellbender_output/cellbender_out.h5 \
                 --cuda

    echo "Cellbender processing completed" > ${sampleName}_cellbender_output/summary.txt
    echo "Cellbender completed for ${sampleName}"
    """
}



// Main workflow
workflow {
    // Log the pipeline parameters
    log.info """
    ===================================
    scQC-flow pipeline
    ===================================
    Mapping directories: ${params.mapping_dirs}
    Output directory: ${params.outputDir}
    GPU acceleration: ${params.gpu}
    ===================================
    """
    
    // Validate required parameters
    if (!params.mapping_dirs) {
        error "Mapping directories file not provided. Please specify --mapping_dirs"
    }
    
    // Read mapping directory information from the CSV file
    Channel
        .fromPath(params.mapping_dirs)
        .splitCsv(header: true)
        .map { row -> tuple(row.samplename, file(row.path)) }
        .set { sampleChannel }

    // Run cellbender on each sample - choose GPU or CPU version based on params.gpu
    if (params.gpu) {
        cellbender_gpu(sampleChannel)
    } else {
        cellbender(sampleChannel)
    }
}
