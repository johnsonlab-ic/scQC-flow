

Aesthetic improvements (Reports)
Mapping diagnostics report; 

[]Called cells vs ddroplet limits table is ugly. make this an interactive DT table, max 10 rows. Same type of table as in Annotation Diagnostics (Top 50 markers)

Ambient RNA report; 

[] first table should also be an interactive table, as per above.

Integration Diagnostics report; 

[] "Barcharts" section; currently, the legend is below all barcharts. this squishes the barcharts when there are too many samples. Place legend on the side, and make legend smaller.
[] "Metadata across UMAP space" section; Never show sample_id as a tab. Too many samples. secondly, don't use hexbins; use pixelbins, as done in scprocess "plot_clusters_annotated_by_densities()". Example function below.


for (meta_var in metadata_vars) {
  cat('####', meta_var, '\n')
  print(plot_clusters_annotated_by_densities(sel_cl_int_dt, meta_var))
  cat('\n\n')
}


Annotation Diagnostics report;


[] Cells used for annotation is wrong, given that doublets are removed.
[] "Heatmaps of marker genes" section; in the Top markers by cluster tab, thhere are too many genes and/or not enough space. Reduce genes to top 5 markers. Legends on the side please. 
[] "Top markers genes by cluster" ; currently each cluster is a tab, and then theres another tabset. let's simplify; do one tab per cluster, and show both the expression dotplots + UMAP feature plots by a patchwork or something similar. p1 + p2 . Reduce to top 5 markers.
[] Move the "Top 50 marker tables by cluster at resolution 0.2" section below the top marker genes by cluster section. 
[] "Cell-type composition by sample" section; the "_harmony_batch" metadata variable is redundant, please remove. It is uninformative. Again, legends for all barcharts on the side, instead of below.
[] "Metadata across UMAP space" is already in integration diagnostics. Remove from this report.
[] Remove the memory usage prints from the rendered report. Write memory usage to a text file for diagnostic usage.


Zoom diagnostics;

In general, the zoom report should be much more similar to the scprocess zoom report. It should show HVG selection, Ambient gene estimation, clustering diagnostics, cluster distribution across samples, qc metrics by cluster, cluster split by metadata variable, highly variable genes, top marker genes. please implement if possible and if the data is being generated through the pipeline. If not, flag whatever can't be implemented in the report as something that needs to be changed in the pipeline.

[] "Selected cells per sample" ; instead of a table, create a single patchworked figure. p1 will be a UMAP of the full dataset, with all cells greyed out except the selected "zoom" value, in colour. p2 will be a barchart of the number of cells taken per individual (labelled) and as usual a black border on barchart. Follow style used in previous barcharts across other reports.
[] "Clusters over UMAP section"; resolution tabset is shown as 1.0 , 0.5, and 0.2 . Please reverse.
[] "Metadata across UMAP space" section; Never show sample_id as a tab. Too many samples. secondly, don't use hexbins; use pixelbins, as done in scprocess "plot_clusters_annotated_by_densities()". Example function below.


for (meta_var in metadata_vars) {
  cat('####', meta_var, '\n')
  print(plot_clusters_annotated_by_densities(sel_cl_int_dt, meta_var))
  cat('\n\n')
}

