nextflow.enable.dsl=2

// Default parameters
params.mapping_dirs = "${projectDir}/personal/mapping_dirs.csv"
params.outputDir = "results"
params.gpu = false
params.report = false
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
      --report              Generate Quarto QC reports for each sample
      --help                Show this help message
      
    Example:
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results --gpu
    nextflow run main.nf --mapping_dirs personal/mapping_dirs.csv --outputDir results --report
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

process quarto_sc_report {
    label "process_medium"
    tag { sampleName }
    container "ah3918/pilot-analyses:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir)

    output:
    path "${sampleName}_qc_report.html"

    script:
    """
    echo "Generating QC report for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"

    # Create QMD file
    cat > ${sampleName}_qc_report.qmd << 'EOF'
---
title: "Seurat QC Report - ${sampleName}"
author: "scQC-flow Pipeline"
date: today
format: 
  html:
    toc: true
    theme: cosmo
    embed-resources: true
execute:
  echo: false
  warning: false
  message: false
---

# Sample Overview: ${sampleName}

**Data path:** `${mappingDir}/outs/filtered_feature_bc_matrix/`

```{r}
#| label: setup

library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(knitr)

# Load the data
data_path <- "${mappingDir}/outs/filtered_feature_bc_matrix/"
count_data <- Read10X(data.dir = data_path)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = count_data, 
                                 project = "${sampleName}", 
                                 min.cells = 3, 
                                 min.features = 200)

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
```

## Basic Statistics

```{r}
#| label: basic-stats

stats <- data.frame(
  Metric = c("Total Cells", "Total Features", "Median UMI per cell", "Median Features per cell"),
  Value = c(
    ncol(seurat_obj),
    nrow(seurat_obj),
    median(seurat_obj\\$nCount_RNA),
    median(seurat_obj\\$nFeature_RNA)
  )
)

kable(stats, caption = "Dataset Statistics")
```

## Quality Control Metrics

```{r}
#| label: qc-metrics

qc_stats <- data.frame(
  Metric = c("Mean UMI per cell", "Mean Features per cell", "Mean % Mitochondrial", "Mean % Ribosomal"),
  Value = c(
    round(mean(seurat_obj\\$nCount_RNA), 2),
    round(mean(seurat_obj\\$nFeature_RNA), 2),
    round(mean(seurat_obj\\$percent.mt), 2),
    round(mean(seurat_obj\\$percent.rb), 2)
  )
)

kable(qc_stats, caption = "Quality Control Metrics")
```

## Visualizations

```{r}
#| label: violin-plots
#| fig-width: 12
#| fig-height: 6

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1)
```

```{r}
#| label: scatter-plots
#| fig-width: 12
#| fig-height: 5

plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_smooth(method = "lm")

# Print plots side by side
library(cowplot)
plot_grid(plot1, plot2, ncol = 2)
```

## Top Expressed Genes

```{r}
#| label: top-genes

# Get top expressed genes
count_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
gene_expression <- Matrix::rowSums(count_matrix)
top_genes <- sort(gene_expression, decreasing = TRUE)[1:10]

top_genes_df <- data.frame(
  Gene = names(top_genes),
  Total_UMI = as.numeric(top_genes),
  Percentage = round((as.numeric(top_genes) / sum(gene_expression)) * 100, 2)
)

kable(top_genes_df, caption = "Top 10 Most Expressed Genes")
```

```{r}
#| label: top-genes-plot
#| fig-width: 10
#| fig-height: 6

ggplot(top_genes_df, aes(x = reorder(Gene, Total_UMI), y = Total_UMI)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Most Expressed Genes", 
       x = "Gene", y = "Total UMI Count") +
  theme_minimal()
```

## Filtering Recommendations

```{r}
#| label: recommendations

min_features <- quantile(seurat_obj\\$nFeature_RNA, 0.01)
max_features <- quantile(seurat_obj\\$nFeature_RNA, 0.99)
max_mt <- quantile(seurat_obj\\$percent.mt, 0.95)

recommendations <- data.frame(
  Filter = c("Minimum features per cell", "Maximum features per cell", "Maximum mitochondrial %"),
  Recommended_Threshold = c(
    round(min_features),
    round(max_features),
    round(max_mt, 1)
  ),
  Rationale = c(
    "Remove low-quality cells (1st percentile)",
    "Remove potential doublets (99th percentile)", 
    "Remove dying cells (95th percentile)"
  )
)

kable(recommendations, caption = "Suggested Filtering Thresholds")
```

## Session Information

```{r}
#| label: session-info

sessionInfo()
```
EOF

    # Render the report
    quarto render ${sampleName}_qc_report.qmd

    echo "QC report completed for ${sampleName}"
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

    // Conditionally run Quarto report generation
    if (params.report) {
        quarto_sc_report(sampleChannel)
    }
}
