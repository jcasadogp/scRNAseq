rm(list = ls()[!grepl("seurat", ls())])

source("~/functions/functions_dge_pathways.R", local = knitr::knit_global())

# ==============================================================================
## Load Seurat Object ##
cluster_res <- "RNA_snn_res.0.3"

## Set and Create Directories for Results and Plots ##
title_prefix <- paste("Velsera", " - ", sep = "")
analysis_directory <- "velsera_object/"

results_directory <- paste0("~/Results/", analysis_directory)
results_plots_directory <- paste0("~/Results_plots/", analysis_directory)
pathway_results_directory <- paste0(results_directory, "pathway_results/")

if(!dir.exists(results_directory)) dir.create(results_directory)
if(!dir.exists(results_plots_directory)) dir.create(results_plots_directory)
if(!dir.exists(pathway_results_directory)) dir.create(pathway_results_directory)

## Databases in Molecular Signature DB to be used ##
msig_databases <- c("C2", "H")

## Create Empty Dataframes ##
all_results_dataframe <- data.frame()
all_results_dataframe_cluster <- data.frame()

main_results_dataframe <- data.frame()
main_results_dataframe_cluster <- data.frame()

parent_results_dataframe <- data.frame()
parent_results_dataframe_cluster <- data.frame()
# ==============================================================================

de_results_dataframe <- read.table(paste0("~/Results/", analysis_directory, "DGE_results/Weighted_DGE_results_total_pseudobulk_CONTRAST.tsv"), header = TRUE)
de_results_dataframe_cluster <- read.table(paste0("~/Results/", analysis_directory, "DGE_results/Weighted_DGE_results_per_cluster_CONTRAST.tsv"), header = TRUE)

for(gr in unique(de_results_dataframe$group)){
  
  de.results.df <- de_results_dataframe %>% dplyr::filter(group == gr)
  de.results.df.cluster <- de_results_dataframe_cluster %>% dplyr::filter(group == gr)
  
  # Determine groups and parameters
  g1 <- unique(de.results.df$group1)
  g2 <- unique(de.results.df$group2)
  
  cat("\n********************************************************************************************", sep="\n")
  cat(paste(g2, "vs", g1, "\n"), sep = "")
  cat("********************************************************************************************\n", sep="\n")
  
  # ===================================================================================================================================================
  # TOTAL PSEUDOBULK
  # ===================================================================================================================================================
  cat("\n=> Total Pseudobulk\n", sep = "\n")
  
  cat("==> Reactome", sep = "\n")
  reactome_list <- getReactome(de.results.df)
  
  all_path_df <- as.data.frame(reactome_list$all_react_list[["Reactome"]])
  # all_path_df$leadingEdge <- NULL
  all_path_df$leadingEdge <- sapply(all_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
  all_path_df <- complete_DE_results(current_de_results = all_path_df, group1 = g1, group2 = g2)
  all_path_df$database <- "Reactome"
  all_results_dataframe <- update_DE_results_df(all_results_dataframe, all_path_df)
  
  main_path_df <- as.data.frame(reactome_list$main_react_list[["Reactome"]])
  if(nrow(main_path_df) > 0){
    # main_path_df$leadingEdge <- NULL
    main_path_df$leadingEdge <- sapply(main_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
    main_path_df <- complete_DE_results(current_de_results = main_path_df, group1 = g1, group2 = g2)
    main_path_df$database <- "Reactome"
    main_results_dataframe <- update_DE_results_df(main_results_dataframe, main_path_df)
  }
  
  parent_path_df <- as.data.frame(reactome_list$parent_react_list[["Reactome"]])
  if(nrow(parent_path_df) > 0) {
    # parent_path_df$leadingEdge <- NULL
    parent_path_df$leadingEdge <- sapply(parent_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
    parent_path_df <- complete_DE_results(current_de_results = parent_path_df, group1 = g1, group2 = g2)
    parent_path_df$database <- "Reactome"
    parent_results_dataframe <- update_DE_results_df(parent_results_dataframe, parent_path_df)
  }

  for(db in msig_databases){
    cat(paste0("==> ", db, "\n"), sep = "")
    
    msigdb_list <- getMSigDbPathways(de.results.df, db)
    
    all_path_df <- as.data.frame(msigdb_list$all_msigdb_list[[db]])
    if(nrow(all_path_df) > 0){
      all_path_df <- complete_DE_results(current_de_results = all_path_df, group1 = g1, group2 = g2)
      all_path_df$leadingEdge <- NA
      # all_path_df$leadingEdge <- sapply(all_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
      all_path_df$pathway_level <- NA
      all_path_df$database <- db
      all_results_dataframe <- update_DE_results_df(all_results_dataframe, all_path_df)
    }

    main_path_df <- as.data.frame(msigdb_list$main_msigdb_list[[db]])
    if(nrow(main_path_df) > 0){
      main_path_df <- complete_DE_results(current_de_results = main_path_df, group1 = g1, group2 = g2)
      main_path_df$leadingEdge <- NA
      # main_path_df$leadingEdge <- sapply(main_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
      main_path_df$pathway_level <- NA
      main_path_df$database <- db
      main_results_dataframe <- update_DE_results_df(main_results_dataframe, main_path_df)
    }

    parent_path_df <- as.data.frame(msigdb_list$parent_msigdb_list[[db]])
    if(nrow(parent_path_df) > 0){
      parent_path_df <- complete_DE_results(current_de_results = parent_path_df, group1 = g1, group2 = g2)
      parent_path_df$leadingEdge <- NA
      # parent_path_df$leadingEdge <- sapply(parent_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
      parent_path_df$pathway_level <- NA
      parent_path_df$database <- db
      parent_results_dataframe <- update_DE_results_df(parent_results_dataframe, parent_path_df)
    }
  }
  
  # ===================================================================================================================================================
  # PER CLUSTER ANALYSIS
  # ===================================================================================================================================================
  cat("\n=> Per Cluster\n", sep = "\n")

  for (cl in sort(unique(de.results.df.cluster$cluster))) {
    cat(paste0("\n*** Cluster ", cl, "\n\n"), sep = "")

    de.results.df <- de.results.df.cluster %>% dplyr::filter(cluster == cl)

    if(nrow(de.results.df) == 0) next

    cluster_res <- unique(de.results.df$clust_res)

    cat("==> Reactome", sep = "\n")
    reactome_list <- getReactome(de.results.df)

    all_path_df <- as.data.frame(reactome_list$all_react_list[["Reactome"]])
    if(nrow(all_path_df) != 0){
      # all_path_df$leadingEdge <- NULL
      all_path_df$leadingEdge <- sapply(all_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
      all_path_df <- complete_DE_results(current_de_results = all_path_df, group1 = g1, group2 = g2, cluster = cl, clust_res = cluster_res)
      all_path_df$database <- "Reactome"
      all_results_dataframe_cluster <- update_DE_results_df(all_results_dataframe_cluster, all_path_df)
    }

    main_path_df <- as.data.frame(reactome_list$main_react_list[["Reactome"]])
    if(nrow(main_path_df) != 0){
      # main_path_df$leadingEdge <- NULL
      main_path_df$leadingEdge <- sapply(main_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
      main_path_df <- complete_DE_results(current_de_results = main_path_df, group1 = g1, group2 = g2, cluster = cl, clust_res = cluster_res)
      main_path_df$database <- "Reactome"
      main_results_dataframe_cluster <- update_DE_results_df(main_results_dataframe_cluster, main_path_df)
    }

    parent_path_df <- as.data.frame(reactome_list$parent_react_list[["Reactome"]])
    if(nrow(parent_path_df) != 0){
      # parent_path_df$leadingEdge <- NULL
      parent_path_df$leadingEdge <- sapply(parent_path_df$leadingEdge, function(x) {paste(as.character(x), collapse = ",")})
      parent_path_df <- complete_DE_results(current_de_results = parent_path_df, group1 = g1, group2 = g2, cluster = cl, clust_res = cluster_res)
      parent_path_df$database <- "Reactome"
      parent_results_dataframe_cluster <- update_DE_results_df(parent_results_dataframe_cluster, parent_path_df)
    }

    for(db in msig_databases){
      cat(paste0("==> ", db, "\n"), sep = "")
      msigdb_list <- getMSigDbPathways(de.results.df, db)
      
      all_path_df <- as.data.frame(msigdb_list$all_react_list[[db]])
      if(nrow(all_path_df) != 0){
        all_path_df <- complete_DE_results(current_de_results = all_path_df, group1 = g1, group2 = g2, cluster = cl, clust_res = cluster_res)
        all_path_df$leadingEdge <- NA
        all_path_df$pathway_level <- NA
        all_path_df$database <- db
        all_results_dataframe_cluster <- update_DE_results_df(all_results_dataframe_cluster, all_path_df)
      }

      main_path_df <- as.data.frame(msigdb_list$main_react_list[[db]])
      if(nrow(main_path_df) != 0){
        main_path_df <- complete_DE_results(current_de_results = main_path_df, group1 = g1, group2 = g2, cluster = cl, clust_res = cluster_res)
        main_path_df$leadingEdge <- NA
        main_path_df$pathway_level <- NA
        main_path_df$database <- db
        main_results_dataframe_cluster <- update_DE_results_df(main_results_dataframe_cluster, main_path_df)
      }

      parent_path_df <- as.data.frame(msigdb_list$parent_react_list[[db]])
      if(nrow(parent_path_df) != 0){
        parent_path_df <- complete_DE_results(current_de_results = parent_path_df, group1 = g1, group2 = g2, cluster = cl, clust_res = cluster_res)
        parent_path_df$leadingEdge <- NA
        parent_path_df$pathway_level <- NA
        parent_path_df$database <- db
        parent_results_dataframe_cluster <- update_DE_results_df(parent_results_dataframe_cluster, parent_path_df)
      }
    } #Close db

  } #Close cluster
  
} #Close main loop



### ===============
### ALL PATHWAYS
### ===============

### Save Pseudobulk Level Pathways
params_for_post_process <- list(
  colnames_for_empty_df = c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "pathway_level", "group1", "group2", "group", "database"),
  ncols_for_empty_df = 11,
  results_file_name = paste0(pathway_results_directory, "Weighted_All_path_results_total_pseudobulk_CONTRAST.tsv")
)

all_results_dataframe <- post_process_and_save_results(all_results_dataframe, params_for_post_process)


### Save Single Cell Level Pathways
params_for_post_process <- list(
  colnames_for_empty_df = c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "pathway_level", "group1", "group2", "group", "database", "cluster", "cluster_res"),
  ncols_for_empty_df = 13,
  results_file_name = paste0(pathway_results_directory, "Weighted_All_path_results_per_cluster_CONTRAST.tsv")
)

all_results_dataframe_cluster <- post_process_and_save_results(all_results_dataframe_cluster, params_for_post_process)





### ===============
### MAIN PATHWAYS
### ===============

### Save Pseudobulk Level Pathways
params_for_post_process <- list(
  colnames_for_empty_df = c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "pathway_level", "group1", "group2", "group", "database"),
  ncols_for_empty_df = 11,
  results_file_name = paste0(pathway_results_directory, "Weighted_Main_path_results_total_pseudobulk_CONTRAST.tsv")
)

main_results_dataframe <- post_process_and_save_results(main_results_dataframe, params_for_post_process)


### Save Single Cell Level Pathways
params_for_post_process <- list(
  colnames_for_empty_df = c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "pathway_level", "group1", "group2", "group", "database", "cluster", "cluster_res"),
  ncols_for_empty_df = 13,
  results_file_name = paste0(pathway_results_directory, "Weighted_Main_path_results_per_cluster_CONTRAST.tsv")
)

main_results_dataframe_cluster <- post_process_and_save_results(main_results_dataframe_cluster, params_for_post_process)





### ===============
### PARENT PATHWAYS
### ===============

### Save Pseudobulk Level Pathways
params_for_post_process <- list(
  colnames_for_empty_df = c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "pathway_level", "group1", "group2", "group", "database"),
  ncols_for_empty_df = 11,
  results_file_name = paste0(pathway_results_directory, "Weighted_Parent_path_results_total_pseudobulk_CONTRAST.tsv")
)

parent_results_dataframe <- post_process_and_save_results(parent_results_dataframe, params_for_post_process)


### Save Single Cell Level Pathways
params_for_post_process <- list(
  colnames_for_empty_df = c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "pathway_level", "group1", "group2", "group", "database", "cluster", "cluster_res"),
  ncols_for_empty_df = 13,
  results_file_name = paste0(pathway_results_directory, "Weighted_Parent_path_results_per_cluster_CONTRAST.tsv")
)

parent_results_dataframe_cluster <- post_process_and_save_results(parent_results_dataframe_cluster, params_for_post_process)

