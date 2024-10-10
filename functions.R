generateSeuratObject <- function(sample_id, data_dir, sample_metadata) {
  
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
  sample.seurat@meta.data[, "Sample"] <- sample_metadata['Firalis.ID']
  sample.seurat@meta.data[, "donor_id"] <- sample_metadata['Donor.ID']
  sample.seurat@meta.data[, "FDS"] <- sample_metadata['FDS']
  sample.seurat@meta.data[, "donor_FDS"] <- paste(sample_metadata['Donor.ID'], sample_metadata['FDS'], sep = "_")
  sample.seurat@meta.data[, "region"] <- sample_metadata['Region']
  sample.seurat@meta.data[, "treatment"] <- sample_metadata['Treatment']
  sample.seurat@meta.data[, "library_chip"] <- sample_metadata['LibraryPreparationChip']
  sample.seurat@meta.data[, "sequencing_batch"] <- sample_metadata['SequencingBatch']
  cat("=> (3/5) Seurat metadata asigned (including donor_FDS)", sep = "\n")
  
  # Compute the percentage of mitochondrial genes:
  sample.seurat[["percent.mt"]] <- PercentageFeatureSet(sample.seurat, pattern = "^MT-")
  cat("=> (4/5) Seurat MT percentage calculated", sep = "\n")
  
  # Save Seurat Object in RDS file:
  saveRDS(sample.seurat, paste0("~/SeuratObjects/PerSample_SeuratObjects/",sample_id,"_before_QC.rds"))
  cat("=> (5/5) Seurat RDS file saved", sep = "\n")
  
  return(sample.seurat)
}

qcSeuratObject <- function(sample.seurat, sample_id, count_max, count_min, feat_max, feat_min, mt_max) {
  
  ## Perform QC steps on Seurat object
  ## Function inputs:
  #' @param sample.seurat: Seurat object before QC.
  #' @param sample_id: Sample identifier of the Seurat object
  #' @param count_max: Maximum number of acceptable counts.
  #' @param count_max: Minimum number of acceptable counts.
  #' @param feat_max: Maximum number of acceptable features. More features indicate duplets in the 10X library generation.
  #' @param feat_min: Maximum number of acceptable features. Few features indicate empty droplets in the 10X library generation.
  #' @param mt_max: Maximum percentage of MT material accepted.
  #' @return sample.seurat QCed object.
  
  cat(sample_id, sep = "\n")
  
  if(missing(count_max)){count_max = max(sample.seurat$nCount_RNA)}
  if(missing(count_min)){count_min = min(sample.seurat$nCount_RNA)}
  if(missing(feat_max)){feat_max = max(sample.seurat$nFeature_RNA)}
  if(missing(feat_min)){feat_min = min(sample.seurat$nFeature_RNA)}
  if(missing(mt_max)){mt_max = max(sample.seurat$percent.mt)}
  
  sample.seurat <- subset(x = sample.seurat, 
                          subset= (nCount_RNA <= count_max) &
                            (nCount_RNA >= count_min) &
                            (nFeature_RNA <= feat_max) &
                            (nFeature_RNA >= feat_min) &
                            (percent.mt <= mt_max))
  cat("=> (1/2) Seurat object filtered", sep = "\n")
  
  saveRDS(sample.seurat, paste0("~/SeuratObjects/PerSample_SeuratObjects/",sample_id,"_after_QC.rds"))
  cat("=> (2/2) Seurat RDS file saved", sep = "\n")
  
  return(sample.seurat)
}

qcViolinPlots <- function(sample.seurat) {
  
  ## Plot violin plots for QC with the number of counts, features and percentage of MT material
  ## Function inputs:
  #' @param sample.seurat: Seurat object before QC.
  #' @return p1, p2, p3: List with the three plots.
  
  p1 <- sample.seurat@meta.data %>% 
    ggplot(aes(y=nCount_RNA, x=orig.ident)) +
    geom_violin(alpha = 1, color="red", fill="red")+
    geom_jitter(color="black", width = 0.5, alpha=0.1, size=1.2) + 
    theme_bw() +
    ggtitle("# UMIs per cell") + 
    theme(plot.title = element_text(size = 9))
  
  p2 <- sample.seurat@meta.data %>% 
    ggplot(aes(y=nFeature_RNA, x=orig.ident)) +
    geom_violin(alpha = 1, color="red", fill="red")+
    geom_jitter(color="black", width = 0.5, alpha=0.1, size=1.2) + 
    theme_bw() +
    ggtitle("# genes per cell") + 
    theme(plot.title = element_text(size = 9))
  
  p3 <- sample.seurat@meta.data %>% 
    ggplot(aes(y=percent.mt, x=orig.ident)) +
    geom_violin(alpha = 1, color="red", fill="red")+
    geom_jitter(color="black", width = 0.5, alpha=0.1, size=1.2) + 
    theme_bw() +
    ggtitle("% MT reads") + 
    theme(plot.title = element_text(size = 9))
  
  return(list(count.violin = p1, feature.violin = p2, mt.violin = p3))
}

getGeneExpression <- function(sample.seurat, gene_name){
  
  metadata_fields <- c("donor_id", "FDS", "treatment")
  
  # Check if gene_name exists in the counts matrix
  if (!(gene_name %in% rownames(sample.seurat@assays[[assay_name]]$counts))) {

    cat(paste0(gene_name, " is not in rownames of Seurat"), sep = "\n")

    # If gene_name doesn't exist, return an empty data frame
    gene_expression <- data.frame(Donor = character(),
                                  FDS = character(),
                                  Treatment = character(),
                                  Timepoint = character(),
                                  positive_cells = integer(),
                                  negative_cells = integer())

    return(gene_expression)
  }

  cat(paste0(gene_name, " is in rownames of Seurat"), sep = "\n")

  gene_positive_counts <- GetAssayData(sample.seurat, layer = "counts")[gene_name, ]
  gene_positive_counts <- gene_positive_counts[gene_positive_counts != 0]
  gene_positive_cell_names <- rownames(as.data.frame(gene_positive_counts))
  gene_positive_metadata <- sample.seurat@meta.data[gene_positive_cell_names,]

  gene_negative_cell_names <- setdiff(colnames(sample.seurat), gene_positive_cell_names)
  gene_negative_counts <- sample.seurat@assays[[assay_name]]$counts[gene_name,gene_negative_cell_names]
  gene_negative_metadata <- sample.seurat@meta.data[gene_negative_cell_names,]

  gene_negative_cells <- as.data.frame(table(gene_negative_metadata[,metadata_fields]))
  colnames(gene_negative_cells) <- c("Donor", "FDS", "Treatment", "negative_cells")

  if(length(gene_positive_cell_names) > 1) {
    gene_positive_cells <- as.data.frame(table(gene_positive_metadata[,metadata_fields]))
    colnames(gene_positive_cells) <- c("Donor", "FDS", "Treatment", "positive_cells")
    gene_expression <- merge(gene_negative_cells, gene_positive_cells, by = c("Donor", "FDS", "Treatment"), all = TRUE)
    gene_expression[is.na(gene_expression)] <- 0 # replace NA with 0
  } else {
    gene_expression <- gene_negative_cells
    gene_expression$positive_cells = 0
  }
  
  gene_expression <- gene_expression[gene_expression$negative_cells != 0 | gene_expression$positive_cells != 0, ]
  gene_expression$total_cells <- gene_expression$positive_cells + gene_expression$negative_cells
  gene_expression$Percentage <- gene_expression$positive_cells*100 / gene_expression$total_cells

  return(gene_expression)
}

getExpressionData <- function(object, gene_list, metadata_fields){
  
  expression_data <- GetAssayData(object, layer = "counts")[gene_list, ]
  expression_data <- t(as.data.frame(expression_data))
  expression_data <- cbind(expression_data, rownames(expression_data))
  rownames(expression_data) <- NULL
  colnames(expression_data)[length(colnames(expression_data))] <- "Cell"
  
  metadata <- as.data.frame(sample.seurat@meta.data)
  metadata <- cbind(metadata, rownames(metadata))
  colnames(metadata)[length(colnames(metadata))] <- "Cell"
  
  data_combined <- merge(expression_data, metadata[,c(metadata_fields, "Cell")], by="Cell")
}

plotDGEVolcano <- function(de.results.df, label_genes, variable_name, x, y, title, filename){
  
  ## Plots the DGE Volcano Plot
  ## Function inputs
  #' @param de.results.df: DF with DGE results
  #' @param label_genes: Vector with subset of genes that will be labeled in the Volcano Plot
  #' @param variable_name name that indicates the genes on the de.results.df
  #' @param x variable on the x axis
  #' @param y variable on the y axis
  #' @param title of the final plot
  #' @return volcano plot
  
  volcanoPlot <- EnhancedVolcano(de.results.df, 
                  lab=de.results.df[,variable_name],
                  x = x, 
                  y = y,
                  selectLab = label_genes,
                  pointSize = 2,
                  labSize = 3,
                  title = title,
                  titleLabSize = 12,
                  subtitle = NULL,
                  pCutoff = 0.05,
                  FCcutoff = 2)
  
  ggsave(filename, volcanoPlot, device = "png", width = 8, height = 5.5, units = "in", dpi = 300)
}

getReactome_old <- function(de.results.df){
  
  ## Function with the pipeline to obtain the reactome object
  ## Function inputs:
  #' @param de.results.df: DF with DGE results
  #' @return list with several dfs
  #'  gene_list: sorted list based on gene ranks (based on the PValue and the direction of the DE)
  #'  final.React: list of genes for fgsea with Reactome ordered from up to downregulated genes
  #'  fgsea.set.React: list with pathways and leading edge genes
  #'  fgsea.res.React: result reactome df with pathwyas, p values, NES values, leading edge genes...
  #'  fgsea.res.tidy.React (IMPORTANT): sort by NES and without pathways that contain no gene
  #'  collapsed.Pathways.React
  #'  main.pathways: pathways that are also on the collapsed pathways
  
  tic("Pathway analysis - Reactome DB")
  
  all_react_list <- list()
  main_react_list <- list()
  parent_react_list <- list()
  
  entrez.db <- org.Hs.eg.db
  g.entrezid <- mapIds(entrez.db,
                       keys = de.results.df$gene_name,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
  
  # set entrez db ids as features
  de.results.df$feature <- g.entrezid
  fgsea.set.React <- reactomePathways(unique(g.entrezid))
  
  # create a sorted list based on gene ranks
  gene_list <- de.results.df
  gene_list$PValue <- ifelse(gene_list$PValue == 0, 1e-100, gene_list$PValue)
  gene_list$fcsign <- sign(gene_list$logFC)
  gene_list$logP <- -log10(gene_list$PValue)
  gene_list$metric <- gene_list$logP/gene_list$fcsign
  
  # prepare list of genes for fgsea with Reactome
  final.React <- gene_list[,c("feature", "metric")]
  final.React <- na.omit(final.React[order(final.React$metric),])
  final.React <- deframe(final.React)
  
  # perform pathway analysis
  suppressWarnings(
    fgsea.res.React <- fgsea(pathways=fgsea.set.React,
                             stats=final.React,
                             minSize=1, 
                             maxSize = Inf,
                             nproc=1)
  )
  
  fgsea.res.tidy.React <- fgsea.res.React %>%
    arrange(desc(NES)) %>%
    dplyr::select(-ES, -log2err)
  
  fgsea.res.tidy.React.filtered <- fgsea.res.tidy.React[ sapply(fgsea.res.tidy.React$leadingEdge, length) > 0 ]
  nrows_filtered <- nrow(fgsea.res.tidy.React) - nrow(fgsea.res.tidy.React.filtered)
  
  fgsea.res.tidy.React.filtered$leadingEdge <- mapIdsList(entrez.db,
                                                          keys=fgsea.res.tidy.React.filtered$leadingEdge,
                                                          keytype="ENTREZID",
                                                          column="SYMBOL")
  
  fgsea.res.tidy.React.filtered$leadingEdge <- as.character(lapply(fgsea.res.tidy.React.filtered$leadingEdge, toString))
  
  suppressWarnings(
    collapsed.Pathways.React <- collapsePathways(fgsea.res.React[order(pval)][pval < 0.01],
                                                 fgsea.set.React, final.React)
  )
  
  main.pathways <- fgsea.res.tidy.React.filtered[fgsea.res.tidy.React.filtered$pathway %in% collapsed.Pathways.React$mainPathways,]
  parent.pathways <- fgsea.res.tidy.React.filtered[fgsea.res.tidy.React.filtered$pathway %in% collapsed.Pathways.React$parentPathways,]
  
  all_react_list[["Reactome"]] <- fgsea.res.tidy.React
  main_react_list[["Reactome"]] <- main.pathways
  parent_react_list[["Reactome"]] <- parent.pathways
  
  toc()
  
  return(list(all_react_list = all_react_list, 
              main_react_list = main_react_list, 
              parent_react_list = parent_react_list))
}

getReactome <- function(de.results.df){
  
  ## Function with the pipeline to obtain the reactome object
  ## It does not call the function reactomePathways() directly, but rather it replicates it and adds some lines to incorporate Reactome pathway levels
  ## Function inputs:
  #' @param de.results.df: DF with DGE results
  #' @return list with several dfs
  #'  gene_list: sorted list based on gene ranks (based on the PValue and the direction of the DE)
  #'  final.React: list of genes for fgsea with Reactome ordered from up to downregulated genes
  #'  fgsea.set.React: list with pathways and leading edge genes
  #'  fgsea.res.React: result reactome df with pathwyas, p values, NES values, leading edge genes...
  #'  fgsea.res.tidy.React (IMPORTANT): sort by NES and without pathways that contain no gene
  #'  collapsed.Pathways.React
  #'  main.pathways: pathways that are also on the collapsed pathways
  
  levelname_assignments <- c(
    "R-HSA-168256"  = "Immune System",
    "R-HSA-1640170" = "Cell Cycle",
    "R-HSA-5653656" = "Vesicle-mediated transport",
    "R-HSA-1430728" = "Metabolism",
    "R-HSA-1643685" = "Disease",
    "R-HSA-392499"  = "Metabolism of proteins",
    "R-HSA-1266738" = "Developmental Biology",
    "R-HSA-8953854" = "Metabolism of RNA",
    "R-HSA-9748784" = "Drug ADME",
    "R-HSA-74160"   = "Gene expression (Transcription)",
    "R-HSA-8953897" = "Cellular responses to stimuli",
    "R-HSA-162582"  = "Signal Transduction",
    "R-HSA-1500931" = "Cell-Cell communication",
    "R-HSA-5357801" = "Programmed Cell Death",
    "R-HSA-109582"  = "Hemostasis",
    "R-HSA-112316"  = "Neuronal System",
    "R-HSA-73894"   = "DNA Repair",
    "R-HSA-1852241" = "Organelle biogenesis and maintenance",
    "R-HSA-382551"  = "Transport of small molecules",
    "R-HSA-9709957" = "Sensory Perception",
    "R-HSA-1474244" = "Extracellular matrix organization",
    "R-HSA-397014"  = "Muscle contraction",
    "R-HSA-4839726" = "Chromatin organization",
    "R-HSA-1474165" = "Reproduction",
    "R-HSA-9612973" = "Autophagy",
    "R-HSA-9609507" = "Protein localization",
    "R-HSA-400253"  = "Circadian Clock",
    "R-HSA-69306"   = "DNA Replication",
    "R-HSA-8963743" = "Digestion and absorption"
  )
  
  stopifnot(requireNamespace("reactome.db"))
  stopifnot(requireNamespace("AnnotationDbi"))
  
  tic("Pathway analysis - Reactome DB")
  
  all_react_list <- list()
  main_react_list <- list()
  parent_react_list <- list()
  
  entrez.db <- org.Hs.eg.db
  g.entrezid <- mapIds(entrez.db,
                       keys = de.results.df$gene_name,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
  
  # set entrez db ids as features
  de.results.df$feature <- g.entrezid
  genes <- unique(g.entrezid)
  
  cat("** before reactomePathways()", "\n", sep = "")
  
  ### fgsea.set.React <- reactomePathways(unique(g.entrezid))
  ### ****************************************************************************************************
  pathways_fun <- na.omit(AnnotationDbi::select(reactome.db::reactome.db, keys = genes, c("PATHID"), keytype = "ENTREZID"))
  pathways_fun <- split(pathways_fun$ENTREZID, pathways_fun$PATHID)
  pathway2name <- as.data.table(AnnotationDbi::select(reactome.db::reactome.db, names(pathways_fun), c("PATHNAME"), "PATHID"))
  
  PATHID = NULL
  pathway2name <- pathway2name[!duplicated(PATHID)]
  
  PATHNAME = NULL
  pathway2name[, `:=`(PATHNAME, sub("^[^:]*: ", "", PATHNAME))]
  
  # Add pathway category
  pth_level_ids <- read.table("~/InputData/Reactome_Files/pathway_level_ids_good.tsv", header = TRUE)
  pth_with_level <- unique(merge(pathway2name, pth_level_ids, by = "PATHID", all.x = TRUE))
  
  pth_with_level$LEVELID <- ifelse(
    is.na(pth_with_level$LEVELID) & pth_with_level$PATHID %in% pth_level_ids$LEVELID,
    pth_with_level$PATHID, 
    pth_with_level$LEVELID
  )
  
  pth_with_level$LEVELNAME <- levelname_assignments[pth_with_level$LEVELID]
  
  # sum(!complete.cases(pth_with_level))
  # # There are 770 pathways without top level pathway
  # pth_without_level <- pth_with_level[!complete.cases(pth_with_level),]
  
  name2pathways <- split(pathway2name$PATHID, pathway2name$PATHNAME)
  pathways_fun <- lapply(name2pathways, function(x) unique(do.call(c, pathways_fun[x])))
  fgsea.set.React <- pathways_fun[!is.na(names(pathways_fun))]
  ## ****************************************************************************************************
  cat("** after reactomePathways() with pathway levels", "\n", sep = "")
  
  # create a sorted list based on gene ranks
  gene_list <- de.results.df
  gene_list$PValue <- ifelse(gene_list$PValue == 0, 1e-100, gene_list$PValue)
  gene_list$fcsign <- sign(gene_list$logFC)
  gene_list$logP <- -log10(gene_list$PValue)
  gene_list$metric <- gene_list$logP/gene_list$fcsign
  
  cat("** gene list", "\n", sep = "")
  
  # prepare list of genes for fgsea with Reactome
  final.React <- gene_list[,c("feature", "metric")]
  final.React <- na.omit(final.React[order(final.React$metric),])
  final.React <- deframe(final.React)
  
  cat("** final.React", "\n", sep = "")
  
  # perform pathway analysis
  suppressWarnings(
    fgsea.res.React <- fgsea(pathways=fgsea.set.React,
                             stats=final.React,
                             minSize=1, 
                             maxSize = Inf,
                             nproc=1)
  )
  
  cat("** fgsea.res.React", "\n", sep = "")
  
  fgsea.res.tidy.React <- fgsea.res.React %>%
    arrange(desc(NES)) %>%
    dplyr::select(-ES, -log2err)
  
  cat("** fgsea.res.tidy.React", "\n", sep = "")
  
  fgsea.res.tidy.React.filtered <- fgsea.res.tidy.React[ sapply(fgsea.res.tidy.React$leadingEdge, length) > 0 ]
  nrows_filtered <- nrow(fgsea.res.tidy.React) - nrow(fgsea.res.tidy.React.filtered)
  
  fgsea.res.tidy.React.filtered$leadingEdge <- mapIdsList(entrez.db,
                                                          keys=fgsea.res.tidy.React.filtered$leadingEdge,
                                                          keytype="ENTREZID",
                                                          column="SYMBOL")
  
  fgsea.res.tidy.React.filtered$leadingEdge <- as.character(lapply(fgsea.res.tidy.React.filtered$leadingEdge, toString))
  
  cat("** fgsea.res.tidy.React.filtered", "\n", sep = "")
  
  suppressWarnings(
    collapsed.Pathways.React <- collapsePathways(fgsea.res.React[order(pval)][pval < 0.01],
                                                 fgsea.set.React, final.React)
  )
  
  cat("** collapsed.Pathways.React", "\n", sep = "")
  
  main.pathways <- fgsea.res.tidy.React.filtered[fgsea.res.tidy.React.filtered$pathway %in% collapsed.Pathways.React$mainPathways,]
  parent.pathways <- fgsea.res.tidy.React.filtered[fgsea.res.tidy.React.filtered$pathway %in% collapsed.Pathways.React$parentPathways,]
  
  # Add pathway level
  fgsea.res.tidy.React <- merge(fgsea.res.tidy.React, pth_with_level[, .(PATHNAME, LEVELNAME)], by.x = "pathway", by.y = "PATHNAME", all.x = TRUE)
  names(fgsea.res.tidy.React)[names(fgsea.res.tidy.React) == "LEVELNAME"] <- "pathway_level"
  
  main.pathways <- merge(main.pathways, pth_with_level[, .(PATHNAME, LEVELNAME)], by.x = "pathway", by.y = "PATHNAME", all.x = TRUE)
  names(main.pathways)[names(main.pathways) == "LEVELNAME"] <- "pathway_level"
  
  parent.pathways <- merge(parent.pathways, pth_with_level[, .(PATHNAME, LEVELNAME)], by.x = "pathway", by.y = "PATHNAME", all.x = TRUE)
  names(parent.pathways)[names(parent.pathways) == "LEVELNAME"] <- "pathway_level"
  
  
  all_react_list[["Reactome"]] <- fgsea.res.tidy.React
  main_react_list[["Reactome"]] <- main.pathways
  parent_react_list[["Reactome"]] <- parent.pathways
  
  toc()
  
  return(list(all_react_list = all_react_list, 
              main_react_list = main_react_list, 
              parent_react_list = parent_react_list))
}

getMSigDbPathways <- function(de.results.df){
  
  ## Function with the pipeline to obtain pathways from MSigDB
  ## Function inputs:
  #' @param de.results.df: DF with DGE results
  #' @return list with several dfs
  #'  gene_list: sorted list based on gene ranks (based on the PValue and the direction of the DE)
  #'  final: list of genes for fgsea with MSigDB ordered from up to downregulated genes
  #'  fgsea.set: list with pathways and leading edge genes
  #'  fgsea.res: result df with pathways, p values, NES values, leading edge genes...
  #'  fgsea.res.tidy (IMPORTANT): sort by NES and without pathways that contain no gene
  
  tic("Pathway analysis - MSigDB C2 and Hallmarks")
  
  # Read Entrez IDs from the Ensembl database
  symbol.db <- EnsDb.Hsapiens.v86
  
  msig.df <- list()
  msig.df[["C2"]] <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = 'CP')
  msig.df[["H"]] <- msigdbr(species = "Homo sapiens", category = "H")
  
  all_msigdb_list <- list()
  main_msigdb_list <- list()
  parent_msigdb_list <- list()
  
  for(db in names(msig.df)){
    
    cat (paste("=>", db), sep = "\n")
    fgsea.set <- msig.df[[db]] %>% split(x = .$gene_symbol, f = .$gs_name)
    
    gene_list <- de.results.df
    
    # calculate metrics for ranking genes
    gene_list$PValue <- ifelse(gene_list$PValue == 0, 1e-100, gene_list$PValue)
    gene_list$fcsign <- sign(gene_list$logFC)
    gene_list$logP <- -log10(gene_list$PValue)
    gene_list$metric <- gene_list$logP/gene_list$fcsign
    
    # prepare list of genes for fgsea with MSigDb
    final <- gene_list[,c("gene_name", "metric")]
    final <- na.omit(final[order(final$metric),])
    final <- deframe(final)
    
    suppressWarnings(
      fgsea.res <- fgsea(pathways=fgsea.set,
                         stats=final,
                         minSize=1, 
                         maxSize = Inf,
                         nproc=1)
    )
    
    fgsea.res.tidy <- fgsea.res %>%
      as_tibble() %>%
      arrange(desc(NES)) %>%
      dplyr::select(-leadingEdge, -ES, -log2err)
    
    suppressWarnings(
      collapsed.Pathways <- collapsePathways(fgsea.res[order(pval)][pval < 0.01],
                                             fgsea.set, final)
    )
    
    main.pathways <- fgsea.res.tidy[fgsea.res.tidy$pathway %in% collapsed.Pathways$mainPathways,]
    parent.pathways <- fgsea.res.tidy[fgsea.res.tidy$pathway %in% collapsed.Pathways$parentPathways,]
    
    all_msigdb_list[[db]] <- fgsea.res.tidy
    main_msigdb_list[[db]] <- main.pathways
    parent_msigdb_list[[db]] <- parent.pathways
  }
  
  toc()
  
  return(list(all_msigdb_list = all_msigdb_list, 
              main_msigdb_list = main_msigdb_list, 
              parent_msigdb_list = parent_msigdb_list))
}

getImportantPathways <- function(pathways) {
  
  ## Returns a subset of pathways containing specific strings
  ## Function inputs:
  #' @param pathways df
  #' @return list with topic specific pathways
  
  # APOPTOSIS (and antiapoptosis) PATHWAY
  apoptosis_PA_1 <- pathways[grep("poptosis", pathways$pathway, ignore.case=TRUE),]
  apoptosis_PA_2 <- pathways[grep("death", pathways$pathway, ignore.case=TRUE),]
  apoptosis_PA <- rbind(apoptosis_PA_1, apoptosis_PA_2)
  
  # INFLAMMATORY PATHWAY sign
  inflammatory_PA <- pathways[grep("inflamm", pathways$pathway, ignore.case=TRUE),]
  inflammatory_PA <- inflammatory_PA[!grepl("anti", inflammatory_PA$pathway, ignore.case=TRUE),]
  
  # ANTI-INFLAMMATORY PATHWAY sign
  anti_inflammatory_PA <- pathways[grep("nti-inflammatory", pathways$pathway, ignore.case=TRUE),]
  
  # IMMUNE PATHWAY sign
  immune_PA_1 <- pathways[grep("mmune", pathways$pathway, ignore.case=TRUE),]
  immune_PA_2 <- pathways[grep("lymph", pathways$pathway, ignore.case=TRUE),]
  immune_PA_3 <- pathways[grep("b_cell", pathways$pathway, ignore.case=TRUE),]
  immune_PA_4 <- pathways[grep("t_cell", pathways$pathway, ignore.case=TRUE),]
  immune_PA <- rbind(immune_PA_1, immune_PA_2, immune_PA_3, immune_PA_4)
  
  # CYTOKINE PATHWAY sign
  cytokine_PA <- pathways[grep("ytokin", pathways$pathway, ignore.case=TRUE),]
  
  # DEFENSIS PATHWAY sign
  defensins_PA <- pathways[grep("efensins", pathways$pathway, ignore.case=TRUE),]
  
  # INTERLEUKINS PATHWAY sign
  interleukin_PA_1 <- pathways[grep("nterleukin", pathways$pathway, ignore.case=TRUE),]
  interleukin_PA_2 <- pathways[grep("IL", pathways$pathway),]
  interleukin_PA <- rbind(interleukin_PA_1, interleukin_PA_2)
  
  # INTERFERON PATHWAY sign
  interferon_PA_1 <- pathways[grep("nterferon", pathways$pathway, ignore.case=TRUE),]
  interferon_PA_2 <- pathways[grep("IFN", pathways$pathway, ignore.case=TRUE),]
  interferon_PA <- rbind(interferon_PA_1, interferon_PA_2)
  
  # TNF PATHWAY sign
  tnf_PA_1 <- pathways[grep("tnf", pathways$pathway, ignore.case=TRUE),]
  tnf_PA <- rbind(tnf_PA_1)
  
  # COX PATHWAY sign
  cox_PA_1 <- pathways[grep("COX", pathways$pathway, ignore.case=TRUE),]
  cox_PA_2 <- pathways[grep("cyclooxygenase", pathways$pathway, ignore.case=TRUE),]
  cox_PA_3 <- pathways[grep("prostagland", pathways$pathway, ignore.case=TRUE),]
  cox_PA <- rbind(cox_PA_1, cox_PA_2, cox_PA_3)
  
  # IDO PATHWAY sign
  ido_PA_1 <- pathways[grep("IDO", pathways$pathway, ignore.case=TRUE),]
  ido_PA_2 <- pathways[grep("indoleamine", pathways$pathway, ignore.case=TRUE),]
  ido_PA <- rbind(ido_PA_1, ido_PA_2)
  
  # MHC PATHWAY sign
  mhc_PA <- pathways[grep("MHC", pathways$pathway, ignore.case=TRUE),]
  
  # EXTRACELLULAR PATHWAY sign
  extracellular_PA_1 <- pathways[grep("xtracellular", pathways$pathway, ignore.case=TRUE),]
  extracellular_PA_2 <- pathways[grep("ECM", pathways$pathway, ignore.case=TRUE),]
  extracellular_PA_3 <- pathways[grep("tissue", pathways$pathway, ignore.case=TRUE),]
  extracellular_PA_4 <- pathways[grep("regeneration", pathways$pathway, ignore.case=TRUE),]
  extracellular_PA_5 <- pathways[grep("metallo", pathways$pathway, ignore.case=TRUE),]
  extracellular_PA_6 <- pathways[grep("membrane", pathways$pathway, ignore.case=TRUE),]
  extracellular_PA_7 <- pathways[grep("matrisome", pathways$pathway, ignore.case=TRUE),]
  extracellular_PA <- rbind(extracellular_PA_1, extracellular_PA_2, extracellular_PA_3, extracellular_PA_4, extracellular_PA_5, extracellular_PA_6,
                            extracellular_PA_7)
  
  collagen_PA <- pathways[grep("ollagen", pathways$pathway, ignore.case=TRUE),]
  
  laminin_PA <- pathways[grep("aminin", pathways$pathway, ignore.case=TRUE),]
  
  glycoprot_PA_1 <- pathways[grep("glycoprot", pathways$pathway, ignore.case=TRUE),]
  glycoprot_PA_2 <- pathways[grep("proteoglyc", pathways$pathway, ignore.case=TRUE),]
  glycoprot_PA <- rbind(glycoprot_PA_1, glycoprot_PA_2)
  
  stress_PA_1 <- pathways[grep("tress", pathways$pathway, ignore.case=TRUE),]
  stress_PA_2 <- pathways[grep("UPR", pathways$pathway),]
  stress_PA_3 <- pathways[grep("ER", pathways$pathway),]
  stress_PA <- rbind(stress_PA_1, stress_PA_2, stress_PA_3)
  
  proliferation_PA_1 <- pathways[grep("cycle", pathways$pathway, ignore.case=TRUE),]
  proliferation_PA_2 <- pathways[grep("cyclins", pathways$pathway, ignore.case=TRUE),]
  proliferation_PA_3 <- pathways[grep("CDK", pathways$pathway, ignore.case=TRUE),]
  proliferation_PA <- rbind(proliferation_PA_1, proliferation_PA_2, proliferation_PA_3)
  
  growth_PA_1 <- pathways[grep("RTK", pathways$pathway, ignore.case=TRUE),]
  growth_PA_2 <- pathways[grep("growth", pathways$pathway, ignore.case=TRUE),]
  growth_PA_3 <- pathways[grep("factor", pathways$pathway, ignore.case=TRUE),]
  growth_PA <- rbind(growth_PA_1, growth_PA_2, growth_PA_3)
  
  cellresponses_PA <- pathways[grep("CD", pathways$pathway, ignore.case=TRUE),]
  dna_repain_PA <- pathways[grep("repair", pathways$pathway, ignore.case=TRUE),]
  unfolded_prot_PA <- pathways[grep("unfolded", pathways$pathway, ignore.case=TRUE),]
  
  PATHWAYS <- rbind(apoptosis_PA,
                    inflammatory_PA,
                    anti_inflammatory_PA,
                    immune_PA,
                    cytokine_PA,
                    defensins_PA,
                    interleukin_PA,
                    interferon_PA,
                    tnf_PA,
                    cox_PA,
                    ido_PA,
                    mhc_PA,
                    extracellular_PA,
                    collagen_PA,
                    glycoprot_PA,
                    laminin_PA,
                    stress_PA,
                    proliferation_PA,
                    growth_PA,
                    cellresponses_PA,
                    dna_repain_PA,
                    unfolded_prot_PA)
  
  return(PATHWAYS)
}

addPathwayCategory <- function(pathways){
  
  ## Returns a the pathways DF with an extra column of pathway category
  ## This category is calculated using grep with specific patterns to be found in the pathway name
  ## Function inputs:
  #' @param pathways df
  #' @return pathways df with extra column
  
  manual_categories = list(
    "tRNA processing" = "Metabolism of RNA",
    "Transport of the SLBP independent Mature mRNA" = "Metabolism of RNA",
    "Transport of Mature mRNAs Derived from Intronless Transcripts" = "Metabolism of RNA",
    "Metabolism of non-coding RNA" = "Metabolism of RNA",
    "mRNA decay from 3' to 5' " = "Metabolism of RNA",
    "Separation of Sister Chromatids" = "Cell Cycle",
    "Resolution of Sister Chromatid Cohesion" = "Cell Cycle",
    "DNA strand elongation" = "Cell Cycle",
    "Deposition of new CENPA-containing nucleosomes at the centromere" = "Cell Cycle",
    "Activation of the pre-replicative complex" = "Cell Cycle",
    "PERK regulates gene expression|Response to metal ions" = "Cellular responses to stimuli",
    "Senescence-Associated Secretory Phenotype (SASP)" = "Cellular responses to stimuli",
    "Mitochondrial translation termination" = "Metabolism of proteins",
    "Formation of a pool of free 40S subunits" = "Metabolism of proteins",
    "Eukaryotic Translation Elongation" = "Metabolism of proteins",
    "Peptide ligand-binding receptors" = "Signal Transduction",
    "MAPK1 (ERK2) activation" = "Signal Transduction",
    "GPCR ligand binding" = "Signal Transduction",
    "Downregulation of ERBB2 signaling" = "Signal Transduction",
    "Peptide chain elongation" = "Metabolism of proteins",
    "Gap junction degradation" = "Vesicle-mediated transport",
    "PINK1-PRKN Mediated Mitophagy" = "Autophagy"
  )
  
  manual_df <- tibble(
    pathway = names(manual_categories),
    manual_category = unlist(manual_categories)
  )
  
  pathways <- pathways %>%
    left_join(manual_df, by = "pathway") %>%
    mutate(pathway_category = case_when(
      !is.na(manual_category) ~ manual_category,
      grepl("SARS-CoV|HIV|APOBEC3G|Leishmania|Disease|Viral|HCMV", pathway, ignore.case = TRUE) ~ "Disease",
      grepl("immune|antigen|B Cell|BCR|TCR|ER-Phagosome|NK|TLR3|CLEC7A|antibody|NOD1/2 Signaling Pathway|NF-kappa-B|Growth hormone receptor signaling|PKR", pathway, ignore.case = TRUE) ~ "Immune System",
      grepl("interleukin|IFN|interferon|TNF", pathway, ignore.case = TRUE) ~ "Immune System",
      grepl("IL", pathway, ignore.case = FALSE) ~ "Immuse System",
      grepl("Developmental Biology", pathway, ignore.case = TRUE) ~ "Developmental Biology",
      grepl("Chemokine receptors bind chemokines|Signal Transduction|RHO|RND2|NF-kB is activated|GTPase cycle|TGF-beta|relaxin|ERBB4", pathway, ignore.case = TRUE) ~ "Signal Transduction",
      grepl("glycosilation|deubiquitination|GTP hydrolysis|Translation|codon", pathway, ignore.case = TRUE) ~ "Metabolism of proteins",
      grepl("vesicle", pathway, ignore.case = TRUE) ~ "Vesicle-mediated transport",
      grepl("Nicotina|fructuose|fatty acids|SREBP|chondroitin", pathway, ignore.case = TRUE) ~ "Metabolism",
      grepl("CoA", pathway, ignore.case = FALSE) ~ "Metabolism",
      grepl("cell cycle|mitosis|G1|G2|mitotic|Synthesis of DNA|meio|checkpoint|chromosome|Telomeres", pathway, ignore.case = TRUE) ~ "Cell Cycle",
      grepl("chromatin", pathway, ignore.case = TRUE) ~ "Chromatin organization",
      grepl("stress|death|apoptosis|necrosis|pyroptosis|repair", pathway, ignore.case = TRUE) ~ "DNA Repair",
      grepl("cellular response|Metallothioneins|(NLR) signaling pathways|chaperone|senescence", pathway, ignore.case = TRUE) ~ "Cellular responses to stimuli",
      grepl("rRNA|exon", pathway, ignore.case = TRUE) ~ "Metabolism of RNA",
      grepl("neural", pathway, ignore.case = TRUE) ~ "Neural System",
      grepl("RNA Polymerase|DNA methylation|RUNX|Transcript", pathway, ignore.case = TRUE) ~ "Gene expression (Transcription)",
      grepl("plasma", pathway, ignore.case = TRUE) ~ "Transport of small molecules",
      grepl("ECM|membrane|extracellular|glycoprotein|proteoglycan|collagen|ROBO|kerat|cell surface|cell junction|syndecan|gap junction", pathway, ignore.case = TRUE) ~ "Extracellular matrix organization",
      grepl("autophagy", pathway, ignore.case = TRUE) ~ "Autophagy",
      grepl("platelet|PECAM1", pathway, ignore.case = TRUE) ~ "Hemostasis",
      grepl("Digestion and absorption", pathway, ignore.case = TRUE) ~ "Digestion and absorption",
      TRUE ~ "Other"
    )) %>%
    dplyr::select(-manual_category)
  
  pathways <- pathways %>%
    mutate(pathway_subcategory = case_when(
      grepl("Inflam", pathway, ignore.case = TRUE) ~ "Inflammatory Response",
      grepl("interleukin", pathway, ignore.case = TRUE) ~ "IL",
      grepl("IL", pathway, ignore.case = FALSE) ~ "IL",
      grepl("IFN|interferon", pathway, ignore.case = TRUE) ~ "IFN",
      grepl("TNF", pathway, ignore.case = TRUE) ~ "TNF",
      grepl("Immun", pathway_category, ignore.case = TRUE) ~ "Immune System",
      TRUE ~ "Other"
    ))
}

preparePathwayData <- function(pathways, group, n_top, n_tail, per_cluster){
  
  ## Prepares the data for further plotting
  ## It takes the whole pathways DF and returns a subsetted one according to group and N (and cluster)
  ## It also reorders it by the sum of NES values
  ## Function inputs:
  #' @param pathways df and restriction parameters
  #' @return pathways df filtered
  
  if(per_cluster){
    # Step 1: Filter by Cluster
    group <- as.character(group)
    all_pathways <- pathways %>%
      dplyr::filter(Group %in% group) %>%
      distinct()
    
    # Step 2: Remove NA values
    all_pathways <- all_pathways %>%
      dplyr::filter(!is.na(NES))
    
    # Step 3: Keep TOP and TAIL pathways
    top_pathways <- all_pathways %>%
      group_by(Group) %>%
      arrange(desc(NES)) %>%
      slice_head(n = n_top) %>%
      arrange(Group)
    
    bottom_pathways <- all_pathways %>%
      group_by(Group) %>%
      arrange(desc(NES)) %>%
      slice_tail(n = n_tail) %>%
      arrange(Group)
    
  } else {
    
    # Step 1: Filter by Cluster
    all_pathways <- pathways %>%
      distinct()
    
    # Step 2: Remove NA values
    all_pathways <- all_pathways %>%
      dplyr::filter(!is.na(NES))
    
    # Step 3: Keep TOP and TAIL pathways
    top_pathways <- all_pathways %>%
      arrange(desc(NES)) %>%
      slice_head(n = n_top)
    
    bottom_pathways <- all_pathways %>%
      arrange(desc(NES)) %>%
      slice_tail(n = n_tail)
  }
  
  pathways_filtered <- bind_rows(top_pathways, bottom_pathways) %>% 
    arrange(NES) %>%
    mutate(order = sprintf("%03d", row_number())) %>%
    distinct()
  
  # Step 4: Calculate the sum of NES values per pathway across all treatment groups => (special for heatmaps)
  pathway_sum <- pathways_filtered %>%
    group_by(pathway) %>%
    summarize(sum_nes = sum(NES)) %>%
    ungroup()
  
  # Reorder the pathways based on the sum of NES values in all treatment groups
  pathways_filtered <- pathways_filtered %>%
    mutate(pathway = factor(pathway, levels = pathway_sum$pathway[order(pathway_sum$sum_nes, decreasing = TRUE)]))
  
  return(pathways_filtered)
}

calculateColumnField <- function(pathways_filtered){
  
  ## Returns a the pathways DF with an extra column called Column in order to divide the plot into 2 or 3
  ## This Column is calculated usign the median of the data and splitting pathway categories in order to reduce the dif between both columns
  ## Function inputs:
  #' @param pathways df
  #' @return pathways df with extra column
  
  xi <- unique(pathways_filtered$Group)
  fi <- table(pathways_filtered$Group)
  fi <- subset(table(pathways_filtered$Group), subset = fi > 0)
  Fi <- cumsum(fi)
  median_position <- sum(fi) / 2
  
  fi <- as.data.frame(fi)
  colnames(fi) <- c("xi", "fi")
  
  cut <- xi[which.max(Fi > median_position)]
  Column_1 <- ifelse(as.numeric(xi) < as.numeric(cut), 1, 2)
  Column_2 <- ifelse(as.numeric(xi) <= as.numeric(cut), 1, 2)
  
  DF <- data.frame(fi, Fi, Column_1, Column_2)
  
  dif_col_1 <- abs(sum(DF[DF$Column_1 == 1,]$fi) - sum(DF[DF$Column_1 == 2,]$fi))
  dif_col_2 <- abs(sum(DF[DF$Column_2 == 1,]$fi) - sum(DF[DF$Column_2 == 2,]$fi))
  
  new_DF <- data.frame(xi)
  if(dif_col_1 <= dif_col_2){
    new_DF <- cbind(new_DF, Column_1)
  } else if(dif_col_1 > dif_col_2){
    new_DF <- cbind(new_DF, Column_2)
  }
  colnames(new_DF) <- c("Group", "Column")
  
  pathways_filtered <- merge(pathways_filtered, new_DF, by = "Group", all.x = TRUE)
  
  return(pathways_filtered)
}

segregateIntoTwoColumns <- function(pathways, facetParam){
  
  ## Returns a the pathways DF with an extra column called 'Column' in order to divide the plot into 2 or 3
  ## 'Column' is calculated splitting pathway categories in order to reduce the difference between both columns
  ## Function inputs:
  #' @param pathways df
  #' @param facetParam the parameter that will be used to facet in the plot and whose values we should split into 2
  #' @return pathways df with 'Column'
  
  # 1. Select the unique values for facetParam field
  xi <- unique(pathways %>% dplyr::select(!!sym(facetParam)) %>% arrange(!!sym(facetParam)))
  # 2. Get frequencies for each value
  fi <- table(pathways[,facetParam])
  # 3. Remove values with freq = 0
  fi <- subset(table(pathways[,facetParam]), subset = fi > 0)
  fi <- as.data.frame(fi)
  fi <- fi[order(-fi$Freq), ]
  
  # 4. Initialize the sum for each group
  sum_group1 <- 0
  sum_group2 <- 0
  
  # 5. Greedily assign elements to the groups - add column 'group' to the fi DF
  for (i in 1:nrow(fi)) {
    if (sum_group1 <= sum_group2) {
      sum_group1 <- sum_group1 + fi$Freq[i]
      fi$group[i] <- "1"
    } else {
      sum_group2 <- sum_group2 + fi$Freq[i]
      fi$group[i] <- "2"
    }
  }
  
  # 6. Changes Var1 for facteParam name
  if("Var1" %in% colnames(fi)){
    names(fi)[names(fi) == "Var1"] <- facetParam
  }
  
  new_DF <- fi %>% dplyr::select(!!sym(facetParam), group)
  colnames(new_DF) <- c(facetParam, "Column")
  
  # 7. Adds 'Column' to the original pathways DF
  pathways <- merge(pathways, new_DF, by = facetParam, all.x = TRUE)
  
  return(pathways)
}

segregateIntoThreeColumns <- function(pathways, facetParam){
  
  ## Returns the pathways DF with an extra column called 'Column' in order to divide the plot into 3
  ## 'Column' is calculated by splitting pathway categories in order to reduce the difference between the three columns
  ## Function inputs:
  #' @param pathways df
  #' @param facetParam the parameter that will be used to facet in the plot and whose values we should split into 3
  #' @return pathways df with 'Column'
  
  # 1. Select the unique values for facetParam field
  xi <- unique(pathways %>% dplyr::select(!!sym(facetParam)) %>% arrange(!!sym(facetParam)))
  
  # 2. Get frequencies for each value
  fi <- table(pathways[,facetParam])
  
  # 3. Remove values with freq = 0
  fi <- subset(table(pathways[,facetParam]), subset = fi > 0)
  fi <- as.data.frame(fi)
  fi <- fi[order(-fi$Freq), ]  # Sort by frequency in descending order
  
  # 4. Initialize the sum for each group
  sum_group1 <- 0
  sum_group2 <- 0
  sum_group3 <- 0
  
  # 5. Greedily assign elements to the groups - add column 'group' to the fi DF
  for (i in 1:nrow(fi)) {
    # Assign the item to the group with the least sum
    if (sum_group1 <= sum_group2 && sum_group1 <= sum_group3) {
      sum_group1 <- sum_group1 + fi$Freq[i]
      fi$group[i] <- "1"
    } else if (sum_group2 <= sum_group1 && sum_group2 <= sum_group3) {
      sum_group2 <- sum_group2 + fi$Freq[i]
      fi$group[i] <- "2"
    } else {
      sum_group3 <- sum_group3 + fi$Freq[i]
      fi$group[i] <- "3"
    }
  }
  
  # 6. Changes Var1 for facetParam name
  if("Var1" %in% colnames(fi)){
    names(fi)[names(fi) == "Var1"] <- facetParam
  }
  
  # 7. Create a new DF with the facetParam and the group column
  new_DF <- fi %>% dplyr::select(!!sym(facetParam), group)
  colnames(new_DF) <- c(facetParam, "Column")
  
  # 8. Merge 'Column' back into the original pathways DF
  pathways <- merge(pathways, new_DF, by = facetParam, all.x = TRUE)
  
  return(pathways)
}

segregateIntoFourColumns <- function(pathways, facetParam){
  
  ## Returns the pathways DF with an extra column called 'Column' in order to divide the plot into 4
  ## 'Column' is calculated by splitting pathway categories in order to reduce the difference between the four columns
  ## Function inputs:
  #' @param pathways df
  #' @param facetParam the parameter that will be used to facet in the plot and whose values we should split into 4
  #' @return pathways df with 'Column'
  
  # 1. Select the unique values for facetParam field
  xi <- unique(pathways %>% dplyr::select(!!sym(facetParam)) %>% arrange(!!sym(facetParam)))
  
  # 2. Get frequencies for each value
  fi <- table(pathways[,facetParam])
  
  # 3. Remove values with freq = 0
  fi <- subset(table(pathways[,facetParam]), subset = fi > 0)
  fi <- as.data.frame(fi)
  fi <- fi[order(-fi$Freq), ]  # Sort by frequency in descending order
  
  # 4. Initialize the sum for each group
  sum_group1 <- 0
  sum_group2 <- 0
  sum_group3 <- 0
  sum_group4 <- 0
  
  # 5. Greedily assign elements to the groups - add column 'group' to the fi DF
  for (i in 1:nrow(fi)) {
    # Assign the item to the group with the least sum
    if (sum_group1 <= sum_group2 && sum_group1 <= sum_group3 && sum_group1 <= sum_group4) {
      sum_group1 <- sum_group1 + fi$Freq[i]
      fi$group[i] <- "1"
    } else if (sum_group2 <= sum_group1 && sum_group2 <= sum_group3 && sum_group2 <= sum_group4) {
      sum_group2 <- sum_group2 + fi$Freq[i]
      fi$group[i] <- "2"
    } else if (sum_group3 <= sum_group1 && sum_group3 <= sum_group2 && sum_group3 <= sum_group4) {
      sum_group3 <- sum_group3 + fi$Freq[i]
      fi$group[i] <- "3"
    } else {
      sum_group4 <- sum_group4 + fi$Freq[i]
      fi$group[i] <- "4"
    }
  }
  
  # 6. Changes Var1 for facetParam name
  if("Var1" %in% colnames(fi)){
    names(fi)[names(fi) == "Var1"] <- facetParam
  }
  
  # 7. Create a new DF with the facetParam and the group column
  new_DF <- fi %>% dplyr::select(!!sym(facetParam), group)
  colnames(new_DF) <- c(facetParam, "Column")
  
  # 8. Merge 'Column' back into the original pathways DF
  pathways <- merge(pathways, new_DF, by = facetParam, all.x = TRUE)
  
  return(pathways)
}

segregateIntoFiveColumns <- function(pathways, facetParam){
  
  ## Returns the pathways DF with an extra column called 'Column' in order to divide the plot into 5
  ## 'Column' is calculated by splitting pathway categories in order to reduce the difference between the five columns
  ## Function inputs:
  #' @param pathways df
  #' @param facetParam the parameter that will be used to facet in the plot and whose values we should split into 5
  #' @return pathways df with 'Column'
  
  # 1. Select the unique values for facetParam field
  xi <- unique(pathways %>% dplyr::select(!!sym(facetParam)) %>% arrange(!!sym(facetParam)))
  
  # 2. Get frequencies for each value
  fi <- table(pathways[,facetParam])
  
  # 3. Remove values with freq = 0
  fi <- subset(table(pathways[,facetParam]), subset = fi > 0)
  fi <- as.data.frame(fi)
  fi <- fi[order(-fi$Freq), ]  # Sort by frequency in descending order
  
  # 4. Initialize the sum for each group
  sum_group1 <- 0
  sum_group2 <- 0
  sum_group3 <- 0
  sum_group4 <- 0
  sum_group5 <- 0
  
  # 5. Greedily assign elements to the groups - add column 'group' to the fi DF
  for (i in 1:nrow(fi)) {
    # Assign the item to the group with the least sum
    if (sum_group1 <= sum_group2 && sum_group1 <= sum_group3 && sum_group1 <= sum_group4 && sum_group1 <= sum_group5) {
      sum_group1 <- sum_group1 + fi$Freq[i]
      fi$group[i] <- "1"
    } else if (sum_group2 <= sum_group1 && sum_group2 <= sum_group3 && sum_group2 <= sum_group4 && sum_group2 <= sum_group5) {
      sum_group2 <- sum_group2 + fi$Freq[i]
      fi$group[i] <- "2"
    } else if (sum_group3 <= sum_group1 && sum_group3 <= sum_group2 && sum_group3 <= sum_group4 && sum_group3 <= sum_group5) {
      sum_group3 <- sum_group3 + fi$Freq[i]
      fi$group[i] <- "3"
    } else if (sum_group4 <= sum_group1 && sum_group4 <= sum_group2 && sum_group4 <= sum_group3 && sum_group4 <= sum_group5) {
      sum_group4 <- sum_group4 + fi$Freq[i]
      fi$group[i] <- "4"
    } else {
      sum_group5 <- sum_group5 + fi$Freq[i]
      fi$group[i] <- "5"
    }
  }
  
  # 6. Changes Var1 for facetParam name
  if("Var1" %in% colnames(fi)){
    names(fi)[names(fi) == "Var1"] <- facetParam
  }
  
  # 7. Create a new DF with the facetParam and the group column
  new_DF <- fi %>% dplyr::select(!!sym(facetParam), group)
  colnames(new_DF) <- c(facetParam, "Column")
  
  # 8. Merge 'Column' back into the original pathways DF
  pathways <- merge(pathways, new_DF, by = facetParam, all.x = TRUE)
  
  return(pathways)
}

plotPathwayHeatmap <- function(pathways, group, n_top = 10, n_tail = 10, title_PA, subtitle_PA, include_significance = FALSE){
  
  ## Plots a heatmap of NES values for specific pathwyas
  ## It takes only the n_top and n_tail pathways of each of the group values
  ## It sorts the pathways in descendent order of the sum of the NES values of the subset of pathways
  ## Function inputs:
  #' @param pathways df
  #' @param group groups to be included as x axis
  #' @param n_top number of pathways with higher NES values to be included
  #' @param n_tail number of pathways with lower NES values to be included
  #' @param title_PA title of the final plot
  #' @param subtitle_PA subtitle of the final plot
  #' @return list with topic specific pathways
  
  pathways_filtered <- preparePathwayData(pathways, group, n_top, n_tail)
  
  # Determine max and min NES values
  max_nes <- ifelse(all(pathways_filtered$NES < 0), 0, max(pathways_filtered$NES))
  min_nes <- ifelse(all(pathways_filtered$NES > 0), 0, min(pathways_filtered$NES))
  
  # Define colors
  if (min_nes >= 0) {
    
    colors <- hcl.colors(20, "RdYlGn")[c(11:20)]
    lims <- c(0, max_nes)
    
  } else if (max_nes <= 0) {
    
    colors <- hcl.colors(20, "RdYlGn")[c(1:10)]
    lims <- c(min_nes, 0)
    
  } else {
    
    colors <- hcl.colors(20, "RdYlGn")
    max_abs <- max(abs(c(min_nes, max_nes)))
    lims <- c(-max_abs, max_abs)
    
  }
  
  if(include_significance) {pathways_filtered$significance <- ifelse(pathways_filtered$padj <= 0.05, "*", "")
  } else {pathways_filtered$significance <- ""}
  
  # Plot
  plot <- ggplot(pathways_filtered, aes(y = reorder(str_wrap(pathway, width = 60), NES), x = Group, fill = NES)) + 
    geom_tile() +
    geom_text(aes(label = significance), color = "black", size = 4, hjust = 0.5, vjust = 0.5) +
    scale_fill_gradientn(colors = colors, limits = lims) +
    ggtitle(title_PA, subtitle = subtitle_PA) +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
          axis.title = element_text(size = 9, face = "bold"),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 15, hjust = 0),
          plot.subtitle = element_text(size = 9, hjust = 0))
  
  return(plot)
}

plotPathwayBarplot <- function(pathways, group = NULL, n_top = 10, n_tail = 10, title_PA, subtitle_PA, include_significance = FALSE, save = FALSE, filename = NULL, per_cluster = FALSE){
  
  ## Plots a heatmap of NES values for specific pathwyas
  ## It takes only the n_top and n_tail pathways of each of the group values
  ## It sorts the pathways in descendent order of the sum of the NES values of the subset of pathways
  ## Function inputs:
  #' @param pathways df
  #' @param group groups to be included as x axis
  #' @param n_top number of pathways with higher NES values to be included
  #' @param n_tail number of pathways with lower NES values to be included
  #' @param title_PA title of the final plot
  #' @param subtitle_PA subtitle of the final plot
  
  pathways_filtered <- preparePathwayData(pathways, group, n_top, n_tail, per_cluster)
  
  pathways_filtered$pathway <- gsub("HALLMARK_", "", pathways_filtered$pathway)
  
  # Determine max and min NES values
  max_nes <- ifelse(all(pathways_filtered$NES < 0), 0, max(pathways_filtered$NES))
  min_nes <- ifelse(all(pathways_filtered$NES > 0), 0, min(pathways_filtered$NES))
  
  # Define colors
  if (min_nes >= 0) {
    
    colors <- hcl.colors(20, "RdYlGn")[c(11:20)]
    lims <- c(0, max_nes)
    
  } else if (max_nes <= 0) {
    
    colors <- hcl.colors(20, "RdYlGn")[c(1:10)]
    lims <- c(min_nes, 0)
    
  } else {
    
    colors <- hcl.colors(20, "RdYlGn")
    max_abs <- max(abs(c(min_nes, max_nes)))
    lims <- c(-max_abs, max_abs)
    
  }
  
  if(include_significance) {pathways_filtered$significance <- ifelse(pathways_filtered$padj <= 0.05, "*", "")
  } else {pathways_filtered$significance <- ""}
  
  if(per_cluster){
    pathways_filtered <- calculateColumnField(pathways_filtered)
    
    plot1 <- ggplot(pathways_filtered[pathways_filtered$Column == 1,], aes(x = gsub("HALLMARK_", "", order), y = NES, fill = NES)) +
      facet_grid(Group~., scales = "free", space = "free") +
      geom_bar(stat = "identity") +
      ylim(lims) +
      scale_fill_gradientn(colors = colors, limits = lims) +
      labs(x = "Pathway", y = "NES") +
      geom_text(aes(label = significance), position = position_stack(vjust = 1), vjust = 0.5, hjust = -0.1, size = 7) +
      theme_bw() %+replace%
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 14),
            strip.text = element_text(size = 13),
            legend.position = "none") +
      scale_x_discrete(
        breaks = pathways_filtered[pathways_filtered$Column == 1,]$order,
        labels = pathways_filtered[pathways_filtered$Column == 1,]$pathway,
        expand = c(0,0)) +
      coord_flip()
    
    plot2 <- ggplot(pathways_filtered[pathways_filtered$Column == 2,], aes(x = gsub("HALLMARK_", "", order), y = NES, fill = NES)) +
      facet_grid(Group~., scales = "free", space = "free") +
      geom_bar(stat = "identity") +
      ylim(lims) +
      scale_fill_gradientn(colors = colors, limits = lims) +
      labs(x = "Pathway", y = "NES") +
      geom_text(aes(label = significance), position = position_stack(vjust = 1), vjust = 0.5, hjust = -0.1, size = 7) +
      theme_bw() %+replace%
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 14),
            strip.text = element_text(size = 13),
            legend.position = "none") +
      scale_x_discrete(
        breaks = pathways_filtered[pathways_filtered$Column == 2,]$order,
        labels = pathways_filtered[pathways_filtered$Column == 2,]$pathway,
        expand = c(0,0)) +
      coord_flip()
    
    plot <- plot1 | plot2
    
  } else {
    
    plot <- ggplot(pathways_filtered, aes(x = gsub("HALLMARK_", "", order), y = NES, fill = NES)) +
      geom_bar(stat = "identity") +
      ylim(lims) +
      scale_fill_gradientn(colors = colors, limits = lims) +
      labs(x = "Pathway", y = "NES") +
      geom_text(aes(label = significance), position = position_stack(vjust = 1), vjust = 0.5, hjust = -0.1, size = 7) +
      theme_bw() %+replace%
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 14),
            strip.text = element_text(size = 13),
            legend.position = "none") +
      scale_x_discrete(
        breaks = pathways_filtered$order,
        labels = pathways_filtered$pathway,
        expand = c(0,0)) +
      coord_flip()
  }
  
  plot <- plot + 
    plot_annotation(title = title_PA, subtitle = subtitle_PA, theme = theme(
      plot.title = element_text(size = 25, hjust = 0, margin = margin(l = 10, t = 10, b = 10)),
      plot.subtitle = element_text(size = 20, hjust = 0, margin = margin(l = 10, t = 10, b = 20))))
  
  if(save){
    n <- max(dim(unique(pathways_filtered[pathways_filtered$Column == 1, c("pathway", "Group")]))[1],
             dim(unique(pathways_filtered[pathways_filtered$Column == 2, c("pathway", "Group")]))[1])
    
    png(file = filename, width = 1000, height = 35*n)
    print(plot)
    dev.off()
  }
  
  return(plot)
}

complete_DE_results_treatment <- function(current_de_results, group1, group2, treatment, cluster = NULL, clust_res = NULL) {
  
  ## Completes the results DF with metadata columns
  ## This is for donor comparability so the subset parameter is TREATMENT
  ## Function inputs:
  #' @param current_de_results
  #' @param group1 - metadata columns
  #' @param group2 - metadata columns
  #' @param treatment - metadata columns
  #' @param cluster - metadata columns
  #' @param clust_res - metadata columns
  #' @return completed DF
  
  current_de_results <- as.data.frame(current_de_results)
  current_de_results$logCPM <- NULL
  current_de_results$F <- NULL
  
  current_de_results$group1 <- group1
  current_de_results$group2 <- group2
  current_de_results$treatment <- treatment
  current_de_results$group <- paste0(current_de_results$group1, "-", current_de_results$group2)
  if (!is.null(cluster) && !is.null(clust_res)) {
    current_de_results$cluster <- cluster
    current_de_results$clust_res <- clust_res
  }
  # current_de_results <- current_de_results[complete.cases(current_de_results), ]
  # Many pathways won't have pathway_level for the moment and that is OK
  current_de_results <- current_de_results[complete.cases(
    if("pathway_level" %in% names(current_de_results)) {
      current_de_results[ , !(names(current_de_results) %in% "pathway_level")]
      } else {
        current_de_results
        }
    ), ]
  
  return(current_de_results)
}

complete_DE_results_donor <- function(current_de_results, group1, group2, donor, cluster = NULL, clust_res = NULL) {
  
  ## Completes the results DF with metadata columns
  ## This is for treatment comparability so the subset parameter is DONOR
  ## Function inputs:
  #' @param current_de_results
  #' @param group1 - metadata columns
  #' @param group2 - metadata columns
  #' @param treatment - metadata columns
  #' @param cluster - metadata columns
  #' @param clust_res - metadata columns
  #' @return completed DF
  
  current_de_results <- as.data.frame(current_de_results)
  current_de_results$logCPM <- NULL
  current_de_results$F <- NULL
  
  current_de_results$group1 <- group1
  current_de_results$group2 <- group2
  current_de_results$donor <- donor
  current_de_results$group <- paste0(current_de_results$group1, "-", current_de_results$group2)
  if (!is.null(cluster) && !is.null(clust_res)) {
    current_de_results$cluster <- cluster
    current_de_results$clust_res <- clust_res
  }
  
  # current_de_results <- current_de_results[complete.cases(current_de_results), ]
  # Many pathways won't have pathway_level for the moment and that is OK
  current_de_results <- current_de_results[complete.cases(
    if("pathway_level" %in% names(current_de_results)) {
      current_de_results[ , !(names(current_de_results) %in% "pathway_level")]
    } else {
      current_de_results
    }
  ), ]
  
  return(current_de_results)
}

update_DE_results_df <- function(de_results_dataframe, current_de_results) {
  
  ## Adds the results to the global results DF
  ## Function inputs:
  #' @param de_results_dataframe - global results
  #' @param current_de_results - current results
  #' @return global resuls merged
  
  if (dim(de_results_dataframe)[1] == 0) {
    de_results_dataframe <- current_de_results
  } else {
    de_results_dataframe <- rbind(de_results_dataframe, current_de_results)
  }
  return(de_results_dataframe)
}

post_process_and_save_results <- function(de_results_dataframe, params_for_post_process) {
  
  ## Saves the global results in accordance to indicated parameters
  ## Function inputs:
  #' @param de_results_dataframe - global results to be saved
  #' @param params_for_post_process - parameters to process, like the saving directory
  #' @return the same DF de_results_dataframe
  
  if (!all(c("ncols_for_empty_df", "colnames_for_empty_df", "results_file_name") %in% names(params_for_post_process))) {
    rlang::abort("Please, provide all required parameters within the list: ncol_for_empty_df, colnames_for_empty_df, results_file_name.")
  }
  
  # Remove unnecessary columns
  de_results_dataframe <- as.data.frame(de_results_dataframe)
  de_results_dataframe$logCPM <- NULL
  de_results_dataframe$F <- NULL
  
  # Check if resulting dataframe is empty
  if (dim(de_results_dataframe)[1] == 0) {
    de_results_dataframe <- data.frame(matrix(vector(), ncol = params_for_post_process$ncol_for_empty_df))
    colnames(de_results_dataframe) <- params_for_post_process$colnames_for_empty_df
    de_metrics_list <- list()
  }
  
  # save results to file
  write.table(de_results_dataframe,
              file = params_for_post_process$results_file_name
  )
  return(de_results_dataframe)
}

getPercentageCoincidentPathways <- function(df, donors, treatments){
  
  ## Generates a DF with the percentage of coincident pathways in terms of regulation sign
  ## Function inputs:
  #' @param df: Pathways DF result of pathway analysis. It contains the following columns: pathway, NES, sign_NEs, group2, donor, database
  #' @param donors List of possible donors
  #' @param treatments List of possible treatments
  #' @return result_final
  
  tic("Get Percentage of Coincident Pathways FUNCTION")
  
  result_final <- data.frame(
    Treatment = character(),
    Donor = character(),
    Donor_sign_NES = numeric(),
    Other_donor = character(),
    OtherDonorPos = numeric(),
    OtherDonorNeg = numeric(),
    OtherNoData = numeric()
  )
  
  for(treatment in treatments){
    cat(paste(treatment, "================================="), sep = "\n")
    df_treat <- df %>% dplyr::filter(group2 == treatment)
    
    for(don in donors){
      
      other_donors <- setdiff(donors, don)
      
      pos_pathways_don <- df_treat %>% dplyr::filter(donor == don) %>% dplyr::filter(sign_NES == 1) %>% dplyr::select(pathway)
      neg_pathways_don <- df_treat %>% dplyr::filter(donor == don) %>% dplyr::filter(sign_NES == -1) %>% dplyr::select(pathway)
      
      other_donors_data <- df_treat %>% dplyr::filter(donor %in% other_donors)
      
      cat(paste0(don, ": UP pathways: ", length(pos_pathways_don$pathway), ", DOWN pathways: ", length(neg_pathways_don$pathway)), sep = "\n")
      
      N <- 2*length(other_donors)
      
      result <- data.frame(
        Treatment = rep(treatment, N),
        Donor = rep(don, N),
        Donor_sign_NES = rep(c(1, -1), length(other_donors)),
        Other_donor = rep(other_donors, each = 2),
        OtherDonorPos = numeric(N),
        OtherDonorNeg = numeric(N),
        OtherNoData = numeric(N)
      )
      
      for(other_don in other_donors){
        
        for(pth in pos_pathways_don$pathway){
          n_pos <- other_donors_data %>% dplyr::filter(donor == other_don & pathway == pth & sign_NES == 1) %>% nrow()
          n_neg <- other_donors_data %>% dplyr::filter(donor == other_don & pathway == pth & sign_NES == -1) %>% nrow()
          
          result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == 1, "OtherDonorPos"] <- 
            result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == 1, "OtherDonorPos"] + n_pos
          
          result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == 1, "OtherDonorNeg"] <- 
            result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == 1, "OtherDonorNeg"] + n_neg
          
          if(n_pos == 0 & n_neg == 0){
            result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == 1, "OtherNoData"] <- 
              result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == 1, "OtherNoData"] + 1
          }
        }
        
        for(pth in neg_pathways_don$pathway){
          n_pos <- other_donors_data %>% dplyr::filter(donor == other_don & pathway == pth & sign_NES == 1) %>% nrow()
          n_neg <- other_donors_data %>% dplyr::filter(donor == other_don & pathway == pth & sign_NES == -1) %>% nrow()
          
          result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == -1, "OtherDonorPos"] <- 
            result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == -1, "OtherDonorPos"] + n_pos
          
          result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == -1, "OtherDonorNeg"] <- 
            result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == -1, "OtherDonorNeg"] + n_neg
          
          if(n_pos == 0 & n_neg == 0){
            result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == -1, "OtherNoData"] <- 
              result[result$Donor == don & result$Other_donor == other_don & result$Donor_sign_NES == -1, "OtherNoData"] + 1
          }
        }
      }
      
      result_final <- rbind(result_final, result)
    } 
  }
  cat(paste("======================================"), sep = "\n")
  toc()
  return(result_final)
}

extract_corr_params <- function(data, x_var, y_var) {
  fit <- lm(as.formula(paste(y_var, "~", x_var)), data = data)
  summary_fit <- summary(fit)
  r_squared <- summary_fit$r.squared
  r <- sqrt(r_squared)
  slope <- summary_fit$coefficients[2, 1]
  r <- ifelse(slope < 0, -r, r)
  p_value <- summary_fit$coefficients[2, 4]
  
  # r_pearson <- cor(data[[x_var]], data[[y_var]])
  # print(paste("R from model: ", round(r, 2)))
  # print(paste("-- Pearson R: ", round(r_pearson, 2)))
  
  return(list(R = r, R2 = r_squared, p_value = p_value))
}

safe_extract_corr_params <- function(data, x_var, y_var) {
  tryCatch({
    fit <- lm(as.formula(paste(y_var, "~", x_var)), data = data)
    summary_fit <- summary(fit)
    r_squared <- summary_fit$r.squared
    r <- sqrt(r_squared)
    p_value <- summary_fit$coefficients[2, 4]
    
    return(list(R = r, R2 = r_squared, p_value = p_value))
  }, error = function(e) {
    return(list())
  }, warning = function(w) {
    return(list())
  })
}
