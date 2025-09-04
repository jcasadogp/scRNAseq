generateSeuratObject <- function(sample_id, data_dir, output_dir, sample_metadata) {
  
  ## Generates a Seurat Object from the 10X files stored in workspace
  ## Function inputs:
  #' @param sample_id: Sample identifier (string).
  #' @param sample_dir Directory containing the 10X files (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz).
  #' @param sample_metadata Metadata to be included in the Seurat object. It must contains columns: Sample.ID, Batch, Donor.ID, treatment and timepoint.
  #' @return sample.seurat object.
  
  cat(sample_id, sep = "\n")
  
  # Data Loading:
  sample.data <- Read10X(data.dir = data_dir)
  cat("=> (1/5) 10X files read", sep = "\n")
  sample.seurat <- CreateSeuratObject(counts = sample.data, project = sample_id, min.cells = 1)
  rm(sample.data)
  cat("=> (2/5) Seurat object created", sep = "\n")
  
  # Setting the metadata for all cells in Seurat Object:
  sample.seurat@meta.data[, "Sample"] <- sample_metadata['SEQUENCED.ID']
  sample.seurat@meta.data[, "Project"] <- sample_metadata['PROJECT']
  sample.seurat@meta.data[, "donor_id"] <- sample_metadata['DONOR']
  sample.seurat@meta.data[, "mnf_step"] <- sample_metadata['MNF.STEP']
  sample.seurat@meta.data[, "fds"] <- sample_metadata['FDS']
  sample.seurat@meta.data[, "treatment"] <- sample_metadata['STIMULATION']
  sample.seurat@meta.data[, "timepoint"] <- sample_metadata['STIMULATION.TIME']
  sample.seurat@meta.data[, "region"] <- sample_metadata['REGION']
  sample.seurat@meta.data[, "library_chip"] <- sample_metadata['LIBRARY.CHIP.IN.PROJECT']
  
  cat("=> (3/5) Seurat metadata asigned (including donor_FDS)", sep = "\n")
  
  # Compute the percentage of mitochondrial genes:
  sample.seurat[["percent.mt"]] <- PercentageFeatureSet(sample.seurat, pattern = "^MT-")
  cat("=> (4/5) Seurat MT percentage calculated", sep = "\n")
  
  # Save Seurat Object in RDS file:
  saveRDS(sample.seurat, paste0(output_dir, sample_id,"_before_QC.rds"))
  cat("=> (5/5) Seurat RDS file saved", sep = "\n")
  
  return(sample.seurat)
}

generateSeuratObjectH5 <- function(sample_id, file_route, output_dir, sample_metadata) {
  
  ## Generates a Seurat Object from the 10X files stored in workspace
  ## Function inputs:
  #' @param sample_id: Sample identifier (string).
  #' @param sample_dir Directory containing the 10X files (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz).
  #' @param sample_metadata Metadata to be included in the Seurat object. It must contains columns: Sample.ID, Batch, Donor.ID, treatment and timepoint.
  #' @return sample.seurat object.
  
  cat(sample_id, sep = "\n")
  
  # Data Loading:
  sample.data <- Read10X_h5(filename = file_route)
  cat("=> (1/5) 10X H5 files read", sep = "\n")
  sample.seurat <- CreateSeuratObject(counts = sample.data, project = sample_id, min.cells = 1)
  rm(sample.data)
  cat("=> (2/5) Seurat object created", sep = "\n")
  
  # Setting the metadata for all cells in Seurat Object:
  sample.seurat@meta.data[, "Sample"] <- sample_metadata['SEQUENCED.ID']
  sample.seurat@meta.data[, "Project"] <- sample_metadata['PROJECT']
  sample.seurat@meta.data[, "donor_id"] <- sample_metadata['DONOR']
  sample.seurat@meta.data[, "mnf_step"] <- sample_metadata['MNF.STEP']
  sample.seurat@meta.data[, "fds"] <- sample_metadata['FDS']
  sample.seurat@meta.data[, "treatment"] <- sample_metadata['STIMULATION']
  sample.seurat@meta.data[, "timepoint"] <- sample_metadata['STIMULATION.TIME']
  sample.seurat@meta.data[, "region"] <- sample_metadata['REGION']
  sample.seurat@meta.data[, "library_chip"] <- sample_metadata['LIBRARY.CHIP.IN.PROJECT']
  
  cat("=> (3/5) Seurat metadata asigned (including donor_FDS)", sep = "\n")
  
  # Compute the percentage of mitochondrial genes:
  sample.seurat[["percent.mt"]] <- PercentageFeatureSet(sample.seurat, pattern = "^MT-")
  cat("=> (4/5) Seurat MT percentage calculated", sep = "\n")
  
  # Save Seurat Object in RDS file:
  saveRDS(sample.seurat, paste0(output_dir, sample_id,"_before_QC.rds"))
  cat("=> (5/5) Seurat RDS file saved", sep = "\n")
  
  return(sample.seurat)
}

qcSeuratObject <- function(sample.seurat, sample_id, output_dir, count_max = NA, count_min = NA, feat_max = NA, feat_min = NA, mt_max = NA) {
  
  ## Perform QC steps on Seurat object
  ## Function inputs:
  #' @param sample.seurat: Seurat object before QC.
  #' @param sample_id: Sample identifier of the Seurat object
  #' @param output_dir: Directory to save the Seurat's after QC.
  #' @param count_max: Maximum number of acceptable counts.
  #' @param count_max: Minimum number of acceptable counts.
  #' @param feat_max: Maximum number of acceptable features. More features indicate duplets in the 10X library generation.
  #' @param feat_min: Maximum number of acceptable features. Few features indicate empty droplets in the 10X library generation.
  #' @param mt_max: Maximum percentage of MT material accepted.
  #' @return sample.seurat QCed object.
  
  cat(sample_id, sep = "\n")
  
  if(is.na(count_max)){count_max = max(sample.seurat$nCount_RNA)}
  if(is.na(count_min)){count_min = min(sample.seurat$nCount_RNA)}
  if(is.na(feat_max)){feat_max = max(sample.seurat$nFeature_RNA)}
  if(is.na(feat_min)){feat_min = min(sample.seurat$nFeature_RNA)}
  if(is.na(mt_max)){mt_max = max(sample.seurat$percent.mt)}
  
  sample.seurat <- subset(x = sample.seurat, 
                          subset= (nCount_RNA <= count_max) &
                            (nCount_RNA >= count_min) &
                            (nFeature_RNA <= feat_max) &
                            (nFeature_RNA >= feat_min) &
                            (percent.mt <= mt_max))
  cat("=> (1/2) Seurat object filtered", sep = "\n")
  
  # saveRDS(sample.seurat, paste0("~/SeuratObjects/PerSample_SeuratObjects/",sample_id,"_after_QC.rds"))
  saveRDS(sample.seurat, paste0(output_dir, "/", sample_id,"_after_QC.rds"))
  cat("=> (2/2) Seurat RDS file saved", sep = "\n")
  
  return(sample.seurat)
}

qcViolinPlots <- function(sample.seurat) {
  
  ## Plot violin plots for QC with the number of counts, features and percentage of MT material
  ## Function inputs:
  #' @param sample.seurat: Seurat object before QC.
  #' @return p1, p2, p3: List with the three plots.
  
  common_part <- list(
    geom_violin(alpha = 1, color = "black", fill = 5, 
                scale = "count",
                draw_quantiles = c(0.25, 0.5, 0.75)),
    # geom_jitter(color = "black", width = 0.5, alpha = 0.025, size = 1.2),
    labs(x = NULL),
    theme_bw()
  )
  
  p1 <- sample.seurat@meta.data %>% 
    ggplot(aes(y=nCount_RNA, x=orig.ident)) +
    ggtitle("# UMIs per cell") + common_part
  
  p2 <- sample.seurat@meta.data %>% 
    ggplot(aes(y=nFeature_RNA, x=orig.ident)) +
    ggtitle("# genes per cell") + common_part
  
  p3 <- sample.seurat@meta.data %>% 
    ggplot(aes(y=percent.mt, x=orig.ident)) +
    ggtitle("% mitochondrial reads") + common_part
  
  return(list(count.violin = p1, feature.violin = p2, mt.violin = p3))
}

qcFinalViolinPlot <- function(plots_list, plot_title){
  
  qc_title <- ggdraw() +
    draw_label(plot_title, fontface = 'bold', x = 0.5, hjust = 0.5, size = 16) +
    theme(
      plot.margin = margin(0, 0, 0, 7),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  final_plot <- plot_grid(qc_title,
                          plot_grid(plots_list$count.violin, 
                                    plots_list$feature.violin, 
                                    plots_list$mt.violin, 
                                    ncol = 3),
                          ncol = 1, rel_heights = c(0.1, 1))
  
  return(final_plot)
}

preprocessSeuratObject <- function(sample.seurat, file_prefix, remove_cell_cycle = FALSE) {
  
  ## This function normalizes, finds most variable features and scales, with cell cycle variables, a Seurat Object
  ## It then calculates PCA and UMAP (dimensionality reduction).
  ## Function inputs:
  #' @param sample.seurat
  #' @return sample.seurat
  
  # Join Layers before preprocess
  sample.seurat[["RNA"]] <- JoinLayers(sample.seurat[["RNA"]])
  
  # Normalize Data
  sample.seurat <- NormalizeData(sample.seurat, 
                                 normalization.method = "LogNormalize", 
                                 scale.factor = 10000)
  
  cat("=> (1/8) Data Normalized", sep = "\n")
  
  # Find Variable Features
  sample.seurat <- FindVariableFeatures(sample.seurat, 
                                        selection.method = "vst", 
                                        nfeatures=2000)
  
  cat("=> (2/8) Most Variable Features calculated", sep = "\n")
  
  # Calculate S and G2M Scores
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  sample.seurat <- CellCycleScoring(sample.seurat, 
                                    s.features = s.genes, 
                                    g2m.features = g2m.genes, set.ident = TRUE)
  cat("=> (3/8) Cell Cycle Scores calculated", sep = "\n")
  
  # Scale Data with Cell Cycle
  if(remove_cell_cycle) {
    
    sample.seurat <- ScaleData(sample.seurat, vars.to.regress = c("S.Score", "G2M.Score"))
    cat("=> (4/8) Data Scaled with Cell Cycle Scores", sep = "\n")
    cat("=> (5/8) Saving intermediate Seurat Object...", sep = "\n")
    saveRDS(sample.seurat, paste0("~/SeuratObjects/Merged_SeuratObjects/",file_prefix,"_2_NFS_no_cc.rds"))
    
  } else {
    
    sample.seurat <- ScaleData(sample.seurat)
    cat("=> (4/8) Data Scaled (no Cell Cycle removal)", sep = "\n")
    cat("=> (5/8) Saving intermediate Seurat Object...", sep = "\n")
    saveRDS(sample.seurat, paste0("~/SeuratObjects/Merged_SeuratObjects/",file_prefix,"_2_NFS_cc.rds"))
    
  }
  
  # Dimensionality Reduction
  sample.seurat <- RunPCA(sample.seurat, verbose = FALSE, npcs = 30)
  cat("=> (6/8) PCA calculated", sep = "\n")
  
  sample.seurat <- RunUMAP(sample.seurat, dims = 1:20)
  cat("=> (7/8) UMAP calculated", sep = "\n")
  
  cat("=> (8/8) Saving final Seurat Object...", sep = "\n")
  if(remove_cell_cycle) {
    
    saveRDS(sample.seurat, paste0("~/SeuratObjects/Merged_SeuratObjects/",file_prefix,"_3_PCA_UMAP_no_cc.rds"))
    
  } else {
    
    saveRDS(sample.seurat, paste0("~/SeuratObjects/Merged_SeuratObjects/",file_prefix,"_3_PCA_UMAP_cc.rds"))
    
  }
  
  return(sample.seurat)
}

# Standard workflow
preprocessSeuratObject_standard <- function(sample.seurat, file_prefix) {
  
  ## This function normalizes, finds most variable features and scales, with cell cycle variables, a Seurat Object
  ## It then calculates PCA and UMAP (dimensionality reduction).
  ## Function inputs:
  #' @param sample.seurat
  #' @return sample.seurat
  
  cat("=> (1/8) Joining layers...", sep = "\n")
  sample.seurat[["RNA"]] <- JoinLayers(sample.seurat[["RNA"]])
  
  cat("=> (2/8) Normalizing Data...", sep = "\n")
  sample.seurat <- NormalizeData(sample.seurat, 
                                 normalization.method = "LogNormalize", 
                                 scale.factor = 10000)
  
  
  cat("=> (3/8) Finding Most Variable Features...", sep = "\n")
  sample.seurat <- FindVariableFeatures(sample.seurat, 
                                        selection.method = "vst", 
                                        nfeatures=2000)
  
  cat("=> (4/8) Calculating Cell Cycle Scores...", sep = "\n")
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  sample.seurat <- CellCycleScoring(sample.seurat, 
                                    s.features = s.genes, 
                                    g2m.features = g2m.genes, set.ident = TRUE)
  
  cat("=> (5/8) Scaling Data (no Cell Cycle removal)...", sep = "\n")
  sample.seurat <- ScaleData(sample.seurat)
  
  cat("=> (6/8) Calculating PCA...", sep = "\n")
  sample.seurat <- RunPCA(sample.seurat, verbose = FALSE, npcs = 50)
  
  cat("=> (7/8) Calculating UMAP...", sep = "\n")
  sample.seurat <- RunUMAP(sample.seurat, dims = 1:30)
  
  cat("=> (8/8) Saving final Seurat Object...", sep = "\n")
  saveRDS(sample.seurat, paste0("~/SeuratObjects/Merged_SeuratObjects/", file_prefix,"_2_standard_preprocess.rds"))
  return(sample.seurat)
}

# Remove Batch Effect
preprocessSeuratObject_batch <- function(sample.seurat, file_prefix) {
  
  ## This function normalizes, finds most variable features and scales, with cell cycle variables, a Seurat Object
  ## It then calculates PCA and UMAP (dimensionality reduction).
  ## Function inputs:
  #' @param sample.seurat
  #' @return sample.seurat
  
  
  saveRDS(sample.seurat, paste0("~/SeuratObjects/Merged_SeuratObjects/", file_prefix,"_2_standard_preprocess.rds"))
  return(sample.seurat)
}

# Remove Cell Cycle Effect
preprocessSeuratObject_cc <- function(sample.seurat, file_prefix) {
  
  ## This function removes Cell Cycle
  ## It then calculates PCA and UMAP (dimensionality reduction).
  ## Function inputs:
  #' @param sample.seurat
  #' @return sample.seurat
  
  cat("=> (1/4) Scaling Data to remove Cell Cycle Effect...", sep = "\n")
  sample.seurat <- ScaleData(sample.seurat, vars.to.regress = c("S.Score", "G2M.Score"))
  
  cat("=> (2/4) PCA Calculation...", sep = "\n")
  sample.seurat <- RunPCA(sample.seurat, verbose = FALSE, npcs = 30)
  
  cat("=> (3/4) UMAP Calculation...", sep = "\n")
  sample.seurat <- RunUMAP(sample.seurat, dims = 1:20)
  
  cat("=> (4/4) Saving final Seurat Object...", sep = "\n")
  saveRDS(sample.seurat, paste0("~/SeuratObjects/Merged_SeuratObjects/", file_prefix,"_4_cell_cycle_removed.rds"))
  return(sample.seurat)
}

# Calculate Clusters
preprocessSeuratObject_cluster <- function(sample.seurat, res, file_prefix) {
  
  ## This function removes Cell Cycle
  ## It then calculates PCA and UMAP (dimensionality reduction).
  ## Function inputs:
  #' @param sample.seurat
  #' @return sample.seurat
  
  cat("=> (1/3) Finding neighbors...", sep = "\n")
  sample.seurat <- FindNeighbors(sample.seurat, , reduction = "pca")
  
  cat("=> (2/3) Finding clusters...", sep = "\n")
  sample.seurat <- FindClusters(sample.seurat, resolution = res)
  
  cat("=> (3/3) Saving final Seurat Object...", sep = "\n")
  saveRDS(sample.seurat, paste0("~/SeuratObjects/Merged_SeuratObjects/", file_prefix,"_5_with_clusters.rds"))
  return(sample.seurat)
}

getGeneExpression <- function(object, gene_name, metadata_fields, assay_name = "RNA") {
  
  #' Get Gene Expression Summary from a Seurat Object
  #'
  #' This function retrieves expression data for a specific gene from a Seurat object and computes 
  #' the number of cells expressing (positive) and not expressing (negative) that gene. The results 
  #' are grouped by specified metadata fields (e.g., donor ID, treatment).
  #'
  #' @param object A Seurat object containing gene expression data and metadata.
  #' @param gene_name A character string specifying the gene of interest.
  #' @param metadata_fields A character vector specifying metadata fields to group the summary by.
  #' @param assay_name The name of the assay within the Seurat object to use (default is "RNA").
  #'
  #' @return A data frame summarizing the number of positive and negative cells, total cells, 
  #' and percentage of positive cells, grouped by the specified metadata fields. 
  #' If the gene is not found, an empty data frame is returned with the expected column structure.
  #'
  #' @export
  
  # Check if the gene exists in the assay
  if (!(gene_name %in% rownames(object@assays[[assay_name]]$counts))) {
    cat(paste0(gene_name, " is not in rownames of Seurat\n"))
    
    # Create and return an empty dataframe with expected columns
    empty_columns <- setNames(replicate(length(metadata_fields), character(), simplify = FALSE),
                              metadata_fields)
    empty_columns$positive_cells <- integer()
    empty_columns$negative_cells <- integer()
    empty_df <- as.data.frame(empty_columns, stringsAsFactors = FALSE)
    return(empty_df)
  }
  
  # Get count matrix for the indicated gene
  gene_counts <- GetAssayData(object, assay = assay_name, layer = "counts")[gene_name, ]
  
  # Determine positive and negative cells
  positive_cells <- names(gene_counts[gene_counts > 0])
  negative_cells <- setdiff(colnames(object), positive_cells)
  
  # Extract metadata for both groups
  positive_meta <- FetchData(object, vars = metadata_fields, cells = positive_cells)
  negative_meta <- FetchData(object, vars = metadata_fields, cells = negative_cells)
  
  # Tabulate counts
  positive_table <- as.data.frame(table(positive_meta))
  colnames(positive_table) <- c(metadata_fields, "positive_cells")
  
  negative_table <- as.data.frame(table(negative_meta))
  colnames(negative_table) <- c(metadata_fields, "negative_cells")
  
  # Merge and fill missing values with 0
  gene_expression <- merge(negative_table, positive_table, by = metadata_fields, all = TRUE)
  gene_expression[is.na(gene_expression)] <- 0
  
  # Add total and percentage
  gene_expression$total_cells <- gene_expression$positive_cells + gene_expression$negative_cells
  gene_expression$percentage_positive <- 100 * gene_expression$positive_cells / gene_expression$total_cells
  
  # Keep only groups with at least one cell
  gene_expression <- gene_expression[gene_expression$total_cells > 0, ]
  
  # Add gene_name
  gene_expression$gene_name <- gene_name
  gene_expression <- gene_expression[, c("gene_name", setdiff(names(gene_expression), "gene_name"))]
  
  return(gene_expression)
}

getExpressionDataDT <- function(object, gene_list = NULL, metadata_fields, assay_name = "RNA") {
  
  #' Extract Gene Expression and Metadata as a Data Table
  #'
  #' This function retrieves expression data for a specified list of genes from a Seurat object, 
  #' merges it with selected metadata, and returns the result as a `data.table`. Genes not present 
  #' in the assay are silently skipped.
  #'
  #' @param object A Seurat object containing gene expression and metadata.
  #' @param gene_list A character vector of gene names to extract. If NULL, all genes will be used.
  #' @param metadata_fields A character vector of metadata column names to include.
  #' @param assay_name The name of the assay to use (default is "RNA").
  #'
  #' @return A `data.table` with one row per cell, containing expression values for the specified genes 
  #' and the selected metadata columns. Cell names are included in a `Cell` column.
  #'
  #' @import data.table
  #' @export
  
  assay_data <- GetAssayData(object, assay = assay_name, layer = "counts")
  
  # If gene_list is provided, subset to valid genes only
  if (!is.null(gene_list)) {
    valid_genes <- intersect(gene_list, rownames(assay_data))
    if (length(valid_genes) == 0) {
      warning("None of the provided genes are present in the assay. Returning metadata only.")
      expr_dt <- data.table(Cell = colnames(assay_data))
    } else {
      assay_data <- assay_data[valid_genes, , drop = FALSE]
      expr_dt <- as.data.table(t(assay_data), keep.rownames = "Cell")
    }
  } else {
    expr_dt <- as.data.table(t(assay_data), keep.rownames = "Cell")
  }
  
  meta_dt <- as.data.table(object@meta.data, keep.rownames = "Cell")[, c("Cell", ..metadata_fields)]
  gene_expression <- merge(expr_dt, meta_dt, by = "Cell", all.x = TRUE)
  
  return(gene_expression)
}
