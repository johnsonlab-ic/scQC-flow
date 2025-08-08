// Reports module for generating QC reports from Cell Ranger outputs
// This module provides Quarto-based report generation

process GENERATE_REPORTS {
    label "process_reports"
    tag { sampleName }
    container "ah3918/pilot-analyses:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(mappingDir), path(dropletqc_metrics), path(scdbl_metrics)

    output:
    tuple val(sampleName), path("${sampleName}_qc_report.html"), emit: html_report
    tuple val(sampleName), path("${sampleName}_qc_report.qmd"), emit: qmd_source

    script:
    """
    echo "Generating QC report for sample: ${sampleName}"
    echo "Mapping directory: ${mappingDir}"
    echo "DropletQC metrics: ${dropletqc_metrics}"
    echo "scDblFinder metrics: ${scdbl_metrics}"

    # Copy the template and replace placeholders
    cp ${projectDir}/modules/reports/template.qmd ${sampleName}_qc_report.qmd
    
    # Replace placeholders with actual values (use absolute path)
    ABSOLUTE_PATH=\$(realpath ${mappingDir})
    sed -i "s|SAMPLE_NAME_PLACEHOLDER|${sampleName}|g" ${sampleName}_qc_report.qmd
    sed -i "s|DATA_PATH_PLACEHOLDER|\$ABSOLUTE_PATH|g" ${sampleName}_qc_report.qmd
    
    # Use the QC metrics files that were passed as input
    DROPLETQC_PATH=\$(realpath ${dropletqc_metrics})
    SCDBL_PATH=\$(realpath ${scdbl_metrics})
    sed -i "s|DROPLETQC_PATH_PLACEHOLDER|\$DROPLETQC_PATH|g" ${sampleName}_qc_report.qmd
    sed -i "s|SCDBL_PATH_PLACEHOLDER|\$SCDBL_PATH|g" ${sampleName}_qc_report.qmd

    # Render the report
    echo "Rendering Quarto report..."
    quarto render ${sampleName}_qc_report.qmd

    echo "QC report completed for ${sampleName}"
    """
}

process COMBINE_REPORTS {
    label "process_reports"
    container "ah3918/pilot-analyses:latest"
    publishDir "${params.outputDir}", mode: 'copy', overwrite: true
    
    input:
    path html_reports
    path qmd_sources

    output:
    path "combined_qc_book/", emit: book_directory
    path "combined_qc_book/_book/", emit: rendered_book

    script:
    """
    echo "Combining all QC reports into a single Quarto book..."
    
    mkdir -p combined_qc_book/chapters
    cd combined_qc_book

    # Copy the book template and configuration
    cp -r ${projectDir}/modules/reports/book_template/* .
    
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
    container "ah3918/pilot-analyses:latest"
    publishDir "${params.outputDir}", mode: 'copy', overwrite: true
    
    input:
    val sample_names
    path mapping_dirs
    path dropletqc_files
    path scdbl_files

    output:
    path "combined_qc_report.html", emit: combined_report
    path "combined_qc_report.qmd", emit: combined_qmd

    script:
    """
    echo "Generating combined QC report for all samples..."
    
    # Debug: List all input files
    echo "All files in work directory:"
    ls -la
    
    # Copy the combined template
    cp ${projectDir}/modules/reports/combined_template.qmd combined_qc_report.qmd
    
    # Create sample information file that the template can read
    echo "sample_name,mapping_dir,dropletqc_file,scdbl_file" > sample_info.csv
    
    # Convert to arrays (handle the join properly)
    SAMPLE_NAMES="${sample_names.join(' ')}"
    
    # Process files in order - they should be named with sample prefixes now
    declare -a sample_array=(\$SAMPLE_NAMES)
    
    for sample in \${sample_array[@]}; do
        echo "Processing sample: \$sample"
        
        # Find mapping directory for this sample
        mapping_dir=""
        for dir in */; do
            if [[ "\$dir" == *"\$sample"* ]]; then
                mapping_dir=\$(realpath "\$dir")
                break
            fi
        done
        
        # Find dropletqc file for this sample
        dropletqc_file=""
        for file in \${sample}_dropletqc_metrics.csv; do
            if [[ -f "\$file" ]]; then
                dropletqc_file=\$(realpath "\$file")
                break
            fi
        done
        
        # Find scdbl file for this sample
        scdbl_file=""
        for file in \${sample}_scdbl_metrics.csv; do
            if [[ -f "\$file" ]]; then
                scdbl_file=\$(realpath "\$file")
                break
            fi
        done
        
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
