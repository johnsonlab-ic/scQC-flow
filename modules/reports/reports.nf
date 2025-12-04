// Reports module for generating QC reports from Cell Ranger outputs
// This module provides Quarto-based report generation

process GENERATE_REPORTS {
    label "process_reports"
    tag { sampleName }
        container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true, pattern: "*.html"

    input:
    tuple val(sampleName), path(mappingDir), path(seurat_preqc_rds), path(seurat_postqc_rds), path(template_qmd), val(max_mito), val(min_nuclear)

    output:
    tuple val(sampleName), path("${sampleName}_qc_report.html"), emit: html_report
    tuple val(sampleName), path("${sampleName}_qc_report.qmd"), emit: qmd_source

    script:
    """
    echo "Generating QC report for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"
    echo "Pre-QC Seurat RDS: ${seurat_preqc_rds}"
    echo "Post-QC Seurat RDS: ${seurat_postqc_rds}"
    echo "QC thresholds: max_mito=${max_mito}, min_nuclear=${min_nuclear}"

    # Copy the template from input path to a new file and replace placeholders
    cp ${template_qmd} ${sampleName}_qc_report.qmd

    # Replace placeholders with actual values (use absolute path)
    ABSOLUTE_PATH=\$(realpath ${mappingDir})
    sed -i "s|SAMPLE_NAME_PLACEHOLDER|${sampleName}|g" ${sampleName}_qc_report.qmd
    sed -i "s|DATA_PATH_PLACEHOLDER|\$ABSOLUTE_PATH|g" ${sampleName}_qc_report.qmd

    # Use the Seurat object paths that were passed as input
    SEURAT_PRE_PATH=\$(realpath ${seurat_preqc_rds})
    SEURAT_POST_PATH=\$(realpath ${seurat_postqc_rds})
    sed -i "s|SEURAT_PRE_PATH_PLACEHOLDER|\$SEURAT_PRE_PATH|g" ${sampleName}_qc_report.qmd
    sed -i "s|SEURAT_POST_PATH_PLACEHOLDER|\$SEURAT_POST_PATH|g" ${sampleName}_qc_report.qmd
    
    # Replace QC threshold placeholders
    sed -i "s|MAX_MITO_PLACEHOLDER|${max_mito}|g" ${sampleName}_qc_report.qmd
    sed -i "s|MIN_NUCLEAR_PLACEHOLDER|${min_nuclear}|g" ${sampleName}_qc_report.qmd

    # Render the report
    echo "Rendering Quarto report..."
    quarto render ${sampleName}_qc_report.qmd

    echo "QC report completed for ${sampleName}"
    """
}

process COMBINE_REPORTS {
        label "process_reports"
        container "ghcr.io/johnsonlab-ic/landmark-sc_image"
        publishDir "${params.outputDir}", mode: 'copy', overwrite: true

        input:
        path html_reports
        path qmd_sources
        path(book_template_dir)

        output:
        path "combined_qc_book/", emit: book_directory
        path "combined_qc_book/_book/", emit: rendered_book

        script:
        """
        echo "Combining all QC reports into a single Quarto book..."

        mkdir -p combined_qc_book/chapters
        cd combined_qc_book

        # Copy the book template and configuration from input
        cp -r ${book_template_dir}/* .

        # Copy all individual QMD reports to chapters directory
        cp ../*.qmd chapters/ || echo "No QMD files to copy"

        # Update the _quarto.yml with all chapters
        echo "Creating book configuration..."

        # Generate chapters list (only if chapters exist)
        if ls chapters/*.qmd >/dev/null 2>&1; then
                CHAPTERS=\$(ls chapters/*.qmd | sed 's|chapters/||g' | sort)
        else
                CHAPTERS=""
        fi

        # Create the chapters section in _quarto.yml
        cat > _quarto.yml << 'EOF'
project:
    type: book
    output-dir: _book

book:
    title: "scQC-flow Quality Control Report"
    author: "scQC-flow Pipeline"
    date: today
  
    chapters:
        - index.qmd
EOF

        # Add each chapter to the YAML
        for chapter in \$CHAPTERS; do
                echo "    - chapters/\$chapter" >> _quarto.yml
        done

        cat >> _quarto.yml << 'EOF'

format:
    html:
        theme: cosmo
        toc: true
        code-fold: true
        code-tools: true
        embed-resources: true

    
execute:
    warning: false
    message: false
EOF

        # Render the book
        echo "Rendering combined QC book..."
        quarto render

        echo "Combined QC book completed"
        """
}

process GENERATE_COMBINED_REPORT {
    label "process_reports"
        container "ghcr.io/johnsonlab-ic/landmark-sc_image"
    publishDir "${params.outputDir}", mode: 'copy', overwrite: true

    input:
    val sample_names
    path mapping_dirs
    path dropletqc_files
    path scdbl_files
    path(combined_template_qmd)

    output:
    path "combined_qc_report.html", emit: combined_report
    path "combined_qc_report.qmd", emit: combined_qmd

    script:
    """
    echo "Generating combined QC report for all samples..."

    # Debug: List all input files
    echo "All files in work directory:"
    ls -la

    # Use the combined template directly from input path
    cp ${combined_template_qmd} combined_qc_report.qmd

    # Create sample information file that the template can read
    echo "sample_name,mapping_dir,dropletqc_file,scdbl_file" > sample_info.csv

    # Convert to arrays (handle the join properly)
    SAMPLE_NAMES="${sample_names.join(' ')}"

    # Get all mapping directories and files in arrays to maintain order
    declare -a sample_array=(\$SAMPLE_NAMES)
    declare -a mapping_dirs=(\$(find . -maxdepth 1 -type l -name "*mapped" | sort))
    declare -a dropletqc_files=(\$(ls *_dropletqc_metrics.csv 2>/dev/null | sort))
    declare -a scdbl_files=(\$(ls *_scdbl_metrics.csv 2>/dev/null | sort))

    # Process files by index to maintain order
    for i in "\${!sample_array[@]}"; do
        sample="\${sample_array[\$i]}"
        echo "Processing sample: \$sample (index \$i)"

        # Get mapping directory by index (since Nextflow preserves order)
        if [ \$i -lt \${#mapping_dirs[@]} ]; then
            mapping_dir=\$(realpath "\${mapping_dirs[\$i]}")
        else
            mapping_dir=""
        fi

        # Get dropletqc file by index
        if [ \$i -lt \${#dropletqc_files[@]} ]; then
            dropletqc_file=\$(realpath "\${dropletqc_files[\$i]}")
        else
            dropletqc_file=""
        fi

        # Get scdbl file by index
        if [ \$i -lt \${#scdbl_files[@]} ]; then
            scdbl_file=\$(realpath "\${scdbl_files[\$i]}")
        else
            scdbl_file=""
        fi

        echo "  Mapping dir: \$mapping_dir"
        echo "  DropletQC file: \$dropletqc_file"
        echo "  scDbl file: \$scdbl_file"

        # Add to CSV
        echo "\$sample,\$mapping_dir,\$dropletqc_file,\$scdbl_file" >> sample_info.csv
    done

    echo "Sample info file contents:"
    cat sample_info.csv

    # Replace placeholder with the sample info file path
    SAMPLE_INFO_PATH=\$(realpath sample_info.csv)
    sed -i "s|SAMPLE_INFO_PLACEHOLDER|\$SAMPLE_INFO_PATH|g" combined_qc_report.qmd

    # Render the combined report
    echo "Rendering combined QC report..."
    quarto render combined_qc_report.qmd

    echo "Combined QC report completed"
    """
}
