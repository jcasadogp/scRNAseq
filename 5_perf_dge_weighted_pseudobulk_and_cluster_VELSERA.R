rm(list = ls()[!grepl("seurat", ls())])

source("~/functions/functions_dge_pathways.R", local = knitr::knit_global())

# ==============================================================================
## Load Seurat Object ##
# seurat.merged <- readRDS("~/SeuratObjects/Merged_SeuratObjects/merged_samples_all_clustered.rds")
seurat.merged$seurat_clusters <- seurat.merged$RNA_snn_res.0.1
cluster_res <- "RNA_snn_res.0.1"

## Set and Create Directories for Results and Plots ##
title_prefix <- paste("Velsera", " - ", sep = "")
analysis_directory <- "velsera_object/"

results_directory <- paste0("~/Results/", analysis_directory)
results_plots_directory <- paste0("~/Results_plots/", analysis_directory)
dge_results_directory <- paste0(results_directory, "DGE_results/")
volcano_plots_directory <- paste0(results_plots_directory, "volcano_plots/")

if(!dir.exists(results_directory)) dir.create(results_directory)
if(!dir.exists(results_plots_directory)) dir.create(results_plots_directory)
if(!dir.exists(dge_results_directory)) dir.create(dge_results_directory)
if(!dir.exists(volcano_plots_directory)) dir.create(volcano_plots_directory)

## Set Parameters fro the Analysis ##
group_condition <- "Batch"
perf_values <- sort(unique(seurat.merged$Batch))
perf_combinations_df <- data.frame(
  group2 = perf_values,
  group1 = sapply(perf_values, function(x) paste(perf_values[perf_values != x], collapse = ""))
)

## Create Empty Dataframes ##
de_results_dataframe <- data.frame()
de_results_dataframe_cluster <- data.frame()
# ==============================================================================

# SCE
sce <- as.SingleCellExperiment(seurat.merged, assay = "RNA")

for(i in seq(nrow(perf_combinations_df))){
  
  # Determine groups and parameters
  group1 <- as.character(perf_combinations_df[i, ]$group1)
  group2 <- as.character(perf_combinations_df[i, ]$group2)
  
  cat("\n********************************************************************************************", sep="\n")
  cat(paste(group2, "vs", group1, "\n"), sep = "")
  cat("********************************************************************************************", sep="\n")
  
  # ===================================================================================================================================================
  # TOTAL PSEUDOBULK
  # ===================================================================================================================================================
  cat("=> Total Pseudobulk", sep = "\n")
  sce_aggregated <- aggregateAcrossCells(sce, id = colData(sce)[, c("FiralisID")], use.assay.type = "counts")
  dim(sce_aggregated)
  colData(sce_aggregated)[, c("FiralisID", "Batch", "Donor", "LibraryPreparationChip")]
  
  # current_sce <- sce_aggregated[,sce_aggregated[[group_condition]] %in% c(group1, group2)]
  # dim(current_sce)
  # colData(current_sce)[, c("FiralisID", "Batch", "Donor", "LibraryPreparationChip")]
  
  current_sce <- sce_aggregated
  current_sce[[group_condition]] <- ifelse(current_sce[[group_condition]] == group2, group2, group1)
  current_sce[[group_condition]] <- relevel(factor(current_sce[[group_condition]]), ref=group1)
  
  de.results <- pseudoBulkDGE(
    current_sce, 
    label = as.factor(rep("all", dim(current_sce)[2])),
    design = as.formula(paste0("~ LibraryPreparationChip + ", group_condition)),
    coef = paste0(group_condition, group2),
    condition = current_sce[[group_condition]]
  )
  
  de.results.df <- de.results@listData[["all"]]
  
  if (is.null(de.results.df)) next
  
  de.results.df <- de.results.df[!is.na(de.results.df$PValue),]
  de.results.df$gene_name <- rownames(de.results.df)
  rownames(de.results.df) <- NULL
  
  de.results.df.filtered <- de.results.df[de.results.df$FDR <= 0.05,]
  de.results.df.filtered.sorted <- de.results.df.filtered[order(de.results.df.filtered$logFC, decreasing=TRUE),]
  
  label_head <- head(de.results.df.filtered.sorted@listData[["gene_name"]], n=10)
  label_tail <- tail(de.results.df.filtered.sorted@listData[["gene_name"]], n=10)
  
  label_genes <- c(label_head, label_tail)
  
  de.results.df <- complete_DE_results(de.results.df, group1, group2)
  de_results_dataframe <- update_DE_results_df(de_results_dataframe, de.results.df)
  
  # Volcano Plot
  filename <- paste0(volcano_plots_directory, group2, "_vs_", group1, "-dge_volcano_plot_total_pseudobulk.png")
  plotDGEVolcano(de.results.df, label_genes, variable_name = "gene_name", x = "logFC", y = "PValue", title = paste0(title_prefix, group2," vs ",group1), filename)
  
  # ===================================================================================================================================================
  # PER CLUSTER ANALYSIS
  # ===================================================================================================================================================
  cat("=> Per Cluster", sep = "\n")
  sce_aggregated <- aggregateAcrossCells(sce, id = colData(sce)[, c("FiralisID", "seurat_clusters")], use.assay.type = "counts")
  
  # current_sce <- sce_aggregated[,sce_aggregated[[group_condition]] %in% c(group1, group2)]
  # dim(current_sce)
  # colData(current_sce)[, c("FiralisID", "Batch", "Donor", "LibraryPreparationChip")]
  
  current_sce <- sce_aggregated
  current_sce[[group_condition]] <- ifelse(current_sce[[group_condition]] == group2, group2, group1)
  current_sce[[group_condition]] <- relevel(factor(current_sce[[group_condition]]), ref=group1)
  
  de.results <- pseudoBulkDGE(
    current_sce, 
    label = current_sce$seurat_clusters,
    design = as.formula(paste0("~ LibraryPreparationChip + ", group_condition)),
    coef = paste0(group_condition, group2),
    condition = current_sce[[group_condition]])
  
  de_results_clusters <- de.results
  
  for (cluster in names(de_results_clusters)) {
    
    print(paste0("Extracting results for cluster ", cluster))
    
    de.results.df <- de_results_clusters@listData[[cluster]]
    
    if (is.null(de.results.df)) next
    
    de.results.df <- de.results.df[!is.na(de.results.df$PValue),]
    de.results.df$gene_name <- rownames(de.results.df)
    rownames(de.results.df) <- NULL
    
    de.results.df.filtered <- de.results.df[de.results.df$FDR <= 0.05,]
    de.results.df.filtered.sorted <- de.results.df.filtered[order(de.results.df.filtered$logFC, decreasing=TRUE),]
    
    label_head <- head(de.results.df.filtered.sorted@listData[["gene_name"]], n=10)
    label_tail <- tail(de.results.df.filtered.sorted@listData[["gene_name"]], n=10)
    
    label_genes <- c(label_head, label_tail)
    
    cluster <- toString(cluster)
    
    de.results.df <- complete_DE_results(current_de_results = de.results.df, group1 = group1, group2 = group2, cluster = cluster, clust_res = cluster_res)
    de_results_dataframe_cluster <- update_DE_results_df(de_results_dataframe_cluster, de.results.df)
    
    filename <- paste0(volcano_plots_directory, group2, "_vs_", group1, "-dge_volcano_plot_cluster_", cluster, ".png")
    
    plotDGEVolcano(de.results.df, label_genes, variable_name = "gene_name", x = "logFC", y = "PValue", 
                   title = paste0(title_prefix, group2," vs ", group1, " - Cluster ", cluster),
                   filename)
    
  } # close cluster
  
} # close main loop - mnf step


## Save Files ##

# Save the whole PSEUDOBULK data in one file
params_for_post_process <- list(
  colnames_for_empty_df = c("logFC", "PValue", "FDR", "group1", "group2", "gene_name"),
  ncols_for_empty_df = 6,
  results_file_name = paste0(dge_results_directory, "Weighted_DGE_results_total_pseudobulk.tsv")
)

de_results_dataframe <- post_process_and_save_results(de_results_dataframe, params_for_post_process)

# Save the whole CLUSTER data in one file
params_for_post_process <- list(
  colnames_for_empty_df = c("logFC", "PValue", "FDR", "group1", "group2", "gene_name", "cluster", "clust_res"),
  ncols_for_empty_df = 8,
  results_file_name = paste0(dge_results_directory, "Weighted_DGE_results_per_cluster.tsv")
)

de_results_dataframe_cluster <- post_process_and_save_results(de_results_dataframe_cluster, params_for_post_process)
