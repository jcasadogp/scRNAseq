rm(list = ls()[!grepl("seurat", ls())])

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
dge_results_directory <- paste0(results_directory, "DGE_results/")
volcano_plots_directory <- paste0(results_plots_directory, "volcano_plots/")

if(!dir.exists(results_directory)) dir.create(results_directory)
if(!dir.exists(results_plots_directory)) dir.create(results_plots_directory)
if(!dir.exists(dge_results_directory)) dir.create(dge_results_directory)
if(!dir.exists(volcano_plots_directory)) dir.create(volcano_plots_directory)

## Set Parameters fro the Analysis ##
group_condition <- "Batch"
group_levels <- sort(unique(seurat.merged$Batch))

## Create Empty Dataframes ##
de_results_dataframe <- data.frame()
de_results_dataframe_cluster <- data.frame()
# ==============================================================================

# SCE
sce <- as.SingleCellExperiment(seurat.merged, assay = "RNA")

for(i in seq(group_levels)){
  
  # Determine groups and parameters
  group_i <- as.character(group_levels[i]) # A
  other_groups <- as.character(group_levels[group_levels != group_i]) # c("B", "C", ...)
  other_groups_prefix <- paste(other_groups, collapse = "") # BC...
  
  cat("\n********************************************************************************************", sep="\n")
  cat(paste(group_i, "vs", paste(other_groups, collapse = " and "), "\n"), sep = "")
  cat("********************************************************************************************\n", sep="\n")
  
  analysis_prefix <- paste(group_i, "vs", other_groups_prefix, sep = "_") # A_vs_BC...
  
  results_plots_directory <- paste0("~/Results_plots/", analysis_directory, analysis_prefix)
  
  if(!dir.exists(results_plots_directory)) dir.create(results_plots_directory)
  
  # ===================================================================================================================================================
  # TOTAL PSEUDOBULK
  # ===================================================================================================================================================
  cat("\n=> Total Pseudobulk\n", sep = "\n")
  
  sce_aggregated <- aggregateAcrossCells(sce, id = colData(sce)[, c("FiralisID")], use.assay.type = "counts")
  dim(sce_aggregated)
  colData(sce_aggregated)[, c("FiralisID", group_condition, "Donor", "LibraryPreparationChip")]
  
  current_sce <- sce_aggregated
  colData(current_sce)[[group_condition]] <- factor(colData(current_sce)[[group_condition]], levels = group_levels)
  current_sce[[group_condition]] <- relevel(factor(current_sce[[group_condition]]), ref=group_i)
  
  design_formula <- as.formula(paste0("~ 0 + LibraryPreparationChip + ", group_condition)) # ~0 + LibraryPreparationChip + Batch
  design_mat <- model.matrix(design_formula, data = colData(current_sce))
  colnames(design_mat) <- c("LibraryPreparationChip", levels(colData(current_sce)[[group_condition]]))
  print(design_mat)
  
  contrast_string <- paste0(group_i, " - (", paste(other_groups, collapse = " + "), ")/", length(other_groups)) # "A - (B + C + ...)/n"
  contrast_formula <- paste0(analysis_prefix, " = ", contrast_string)
  contrast_mat <- makeContrasts(contrasts = contrast_formula, levels = design_mat)
  print(contrast_mat)
  
  de.results <- pseudoBulkDGE(current_sce, 
                              label = analysis_prefix,
                              design = design_formula,
                              coef = analysis_prefix,
                              contrast = contrast_mat)
  
  
  de.results.df <- de.results@listData[[analysis_prefix]]
  
  if (is.null(de.results.df)) next
  
  de.results.df <- de.results.df[!is.na(de.results.df$PValue),]
  de.results.df$gene_name <- rownames(de.results.df)
  rownames(de.results.df) <- NULL
  
  de.results.df.filtered <- de.results.df[de.results.df$FDR <= 0.05,]
  de.results.df.filtered.sorted <- de.results.df.filtered[order(de.results.df.filtered$logFC, decreasing=TRUE),]
  
  label_head <- head(de.results.df.filtered.sorted@listData[["gene_name"]], n=10)
  label_tail <- tail(de.results.df.filtered.sorted@listData[["gene_name"]], n=10)
  
  label_genes <- c(label_head, label_tail)
  
  de.results.df <- complete_DE_results(de.results.df, group1 = other_groups_prefix, group2 = group_i)
  de_results_dataframe <- update_DE_results_df(de_results_dataframe, de.results.df)
  
  # Volcano Plot
  filename <- paste0(volcano_plots_directory, group_i, "_vs_", other_groups_prefix, "-dge_volcano_plot_total_pseudobulk.png")
  plotDGEVolcano(de.results.df, label_genes, variable_name = "gene_name", x = "logFC", y = "PValue", title = paste0(title_prefix, group_i," vs ", other_groups_prefix), filename)
  
  # ===================================================================================================================================================
  # PER CLUSTER ANALYSIS
  # ===================================================================================================================================================
  cat("\n=> Per Cluster\n", sep = "\n")
  
  sce_aggregated <- aggregateAcrossCells(sce, id = colData(sce)[, c("FiralisID", "seurat_clusters")],use.assay.type = "counts")
  dim(sce_aggregated)
  colData(sce_aggregated)[, c("FiralisID", group_condition, "Donor", "LibraryPreparationChip", "seurat_clusters")]
  
  
  current_sce <- sce_aggregated
  colData(current_sce)[[group_condition]] <- factor(colData(current_sce)[[group_condition]], levels = group_levels)
  current_sce[[group_condition]] <- relevel(factor(current_sce[[group_condition]]), ref=group_i)
  
  design_formula <- as.formula(paste0("~ 0 + LibraryPreparationChip + ", group_condition)) # ~0 + LibraryPreparationChip + Batch
  design_mat <- model.matrix(design_formula, data = colData(current_sce))
  colnames(design_mat) <- c("LibraryPreparationChip", levels(colData(current_sce)[[group_condition]]))
  print(design_mat)
  
  contrast_string <- paste0(group_i, " - (", paste(other_groups, collapse = " + "), ")/", length(other_groups)) # "A - (B + C + ...)/n"
  contrast_formula <- paste0(analysis_prefix, " = ", contrast_string)
  contrast_mat <- makeContrasts(contrasts = contrast_formula, levels = design_mat)
  print(contrast_mat)
  
  de.results <- pseudoBulkDGE(current_sce, 
                              label = current_sce$seurat_clusters,
                              design = design_formula,
                              coef = analysis_prefix,
                              contrast = contrast_mat)

  de_results_clusters <- de.results

  for (cluster in names(de_results_clusters)) {

    cat(paste0("\n*** Cluster ", cluster, "\n\n"), sep = "")

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

    de.results.df <- complete_DE_results(current_de_results = de.results.df, group1 = other_groups_prefix, group2 = group_i, cluster = cluster, clust_res = cluster_res)
    de_results_dataframe_cluster <- update_DE_results_df(de_results_dataframe_cluster, de.results.df)

    filename <- paste0(volcano_plots_directory, group_i, "_vs_", other_groups_prefix, "-dge_volcano_plot_cluster_", cluster, ".png")
    plotDGEVolcano(de.results.df, label_genes, variable_name = "gene_name", x = "logFC", y = "PValue",
                   title = paste0(title_prefix, group_i," vs ", other_groups_prefix, " - Cluster ", cluster),
                   filename)

  } # close cluster
  
} # close main loop - comparison

# Save the whole PSEUDOBULK data in one file
params_for_post_process <- list(
  colnames_for_empty_df = c("logFC", "PValue", "FDR", "group1", "group2", "gene_name"),
  ncols_for_empty_df = 6,
  results_file_name = paste0(dge_results_directory, "Weighted_DGE_results_total_pseudobulk_CONTRAST.tsv")
)

de_results_dataframe <- post_process_and_save_results(de_results_dataframe, params_for_post_process)

# Save the whole CLUSTER data in one file
params_for_post_process <- list(
  colnames_for_empty_df = c("logFC", "PValue", "FDR", "group1", "group2", "gene_name", "cluster", "clust_res"),
  ncols_for_empty_df = 8,
  results_file_name = paste0(dge_results_directory, "Weighted_DGE_results_per_cluster_CONTRAST.tsv")
)

de_results_dataframe_cluster <- post_process_and_save_results(de_results_dataframe_cluster, params_for_post_process)