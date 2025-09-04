rm(list = ls()[!grepl("seurat|x.genes", ls())])

source("~/functions/functions_dge_pathways.R", local = knitr::knit_global())

# ==============================================================================
## Load Seurat Object ##
seurat.merged <- readRDS("~/SeuratObjects/Merged_SeuratObjects/merged_samples_all_clustered.rds")
seurat.merged$seurat_clusters <- seurat.merged$RNA_snn_res.0.1
cluster_res <- "RNA_snn_res.0.1"

## Set and Create Directories for Results and Plots ##
title_prefix <- paste("Velsera", " - ", sep = "")
analysis_directory <- "velsera_object/"

results_directory <- paste0("~/Results/", analysis_directory)
results_plots_directory <- paste0("~/Results_plots/", analysis_directory)
markers_results_directory <- paste0(results_directory, "Cluster_analysis/")

if(!dir.exists(results_directory)) dir.create(results_directory)
if(!dir.exists(results_plots_directory)) dir.create(results_plots_directory)
if(!dir.exists(markers_results_directory)) dir.create(markers_results_directory)

## Databases in Molecular Signature DB to be used ##
msig_databases <- c("C2", "H")

## Create Empty Dataframes ##
all_results_dataframe <- data.frame()
main_results_dataframe <- data.frame()
parent_results_dataframe <- data.frame()
# ==============================================================================


# # ## DGE between clusters ##
# Idents(seurat.merged) <- "seurat_clusters"
# 
# seurat.merged.markers <- FindAllMarkers(seurat.merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25) %>%
#   dplyr::select(-gene) %>% rownames_to_column(var="gene")
# 
# write.table(seurat.merged.markers, file = paste0(markers_results_directory, "all_clusters_markers.tsv"))



## Pathway Analysis ##
x.genes <- wilcoxauc(seurat.merged, "seurat_clusters")

colnames(x.genes)[colnames(x.genes) == "feature"] <- "gene_name"
colnames(x.genes)[colnames(x.genes) == "pval"] <- "PValue"

for(cl in sort(unique(x.genes$group))) {
  cat(paste0("\n*** Cluster ", cl, "\n\n"), sep = "")
  
  clstr.genes <- x.genes %>% dplyr::filter(group == cl)
  
  reactome_list <- getReactome(de.results.df = clstr.genes, rank_metric = "auc")
  
  all_path_df <- as.data.frame(reactome_list$all_react_list[["Reactome"]])
  # all_path_df$leadingEdge <- NULL
  all_path_df$leadingEdge <- sapply(all_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
  all_path_df <- complete_DE_results(current_de_results = all_path_df, group1 = "other_cl", group2 = cl)
  all_path_df$database <- "Reactome"
  all_results_dataframe <- update_DE_results_df(all_results_dataframe, all_path_df)
  
  main_path_df <- as.data.frame(reactome_list$main_react_list[["Reactome"]])
  if(nrow(main_path_df) > 0){
    # main_path_df$leadingEdge <- NULL
    main_path_df$leadingEdge <- sapply(main_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
    main_path_df <- complete_DE_results(current_de_results = main_path_df, group1 = "other_cl", group2 = cl)
    main_path_df$database <- "Reactome"
    main_results_dataframe <- update_DE_results_df(main_results_dataframe, main_path_df)
  }
  
  parent_path_df <- as.data.frame(reactome_list$parent_react_list[["Reactome"]])
  if(nrow(parent_path_df) > 0) {
    # parent_path_df$leadingEdge <- NULL
    parent_path_df$leadingEdge <- sapply(parent_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
    parent_path_df <- complete_DE_results(current_de_results = parent_path_df, group1 = "other_cl", group2 = cl)
    parent_path_df$database <- "Reactome"
    parent_results_dataframe <- update_DE_results_df(parent_results_dataframe, parent_path_df)
  }
  
  for(db in msig_databases){
    cat(paste0("==> ", db, "\n"), sep = "")
    
    msigdb_list <- getMSigDbPathways(de.results.df = clstr.genes, db, rank_metric = "auc")
    
    all_path_df <- as.data.frame(msigdb_list$all_msigdb_list[[db]])
    if(nrow(all_path_df) > 0){
      all_path_df <- complete_DE_results(current_de_results = all_path_df, group1 = "other_cl", group2 = cl)
      all_path_df$leadingEdge <- NA
      # all_path_df$leadingEdge <- sapply(all_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
      all_path_df$pathway_level <- NA
      all_path_df$database <- db
      all_results_dataframe <- update_DE_results_df(all_results_dataframe, all_path_df)
    }
    
    main_path_df <- as.data.frame(msigdb_list$main_msigdb_list[[db]])
    if(nrow(main_path_df) > 0){
      main_path_df <- complete_DE_results(current_de_results = main_path_df, group1 = "other_cl", group2 = cl)
      main_path_df$leadingEdge <- NA
      # main_path_df$leadingEdge <- sapply(main_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
      main_path_df$pathway_level <- NA
      main_path_df$database <- db
      main_results_dataframe <- update_DE_results_df(main_results_dataframe, main_path_df)
    }
    
    parent_path_df <- as.data.frame(msigdb_list$parent_msigdb_list[[db]])
    if(nrow(parent_path_df) > 0){
      parent_path_df <- complete_DE_results(current_de_results = parent_path_df, group1 = "other_cl", group2 = cl)
      parent_path_df$leadingEdge <- NA
      # parent_path_df$leadingEdge <- sapply(parent_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
      parent_path_df$pathway_level <- NA
      parent_path_df$database <- db
      parent_results_dataframe <- update_DE_results_df(parent_results_dataframe, parent_path_df)
    }
  }
  
}


### ===============
### ALL PATHWAYS
### ===============

### Save Pseudobulk Level Pathways
params_for_post_process <- list(
  colnames_for_empty_df = c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "pathway_level", "group1", "group2", "group", "database"),
  ncols_for_empty_df = 11,
  results_file_name = paste0(markers_results_directory, "All_path_results_total_pseudobulk.tsv")
)

all_results_dataframe <- post_process_and_save_results(all_results_dataframe, params_for_post_process)




### ===============
### MAIN PATHWAYS
### ===============

### Save Pseudobulk Level Pathways
params_for_post_process <- list(
  colnames_for_empty_df = c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "pathway_level", "group1", "group2", "group", "database"),
  ncols_for_empty_df = 11,
  results_file_name = paste0(markers_results_directory, "Main_path_results_total_pseudobulk.tsv")
)

main_results_dataframe <- post_process_and_save_results(main_results_dataframe, params_for_post_process)




### ===============
### PARENT PATHWAYS
### ===============

### Save Pseudobulk Level Pathways
params_for_post_process <- list(
  colnames_for_empty_df = c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "pathway_level", "group1", "group2", "group", "database"),
  ncols_for_empty_df = 11,
  results_file_name = paste0(markers_results_directory, "Parent_path_results_total_pseudobulk.tsv")
)

parent_results_dataframe <- post_process_and_save_results(parent_results_dataframe, params_for_post_process)