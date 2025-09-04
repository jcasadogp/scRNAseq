plotDGEVolcano <- function(de.results.df, label_genes, variable_name, x, y, title, filename = NULL){
  
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
                                 FCcutoff = 2,
                                 drawConnectors = TRUE,
                                 max.overlaps = Inf)
  
  if(!is.null(filename)){
    ggsave(filename, volcanoPlot, device = "png", width = 7, height = 6, units = "in", dpi = 300) 
  } else {
    return(volcanoPlot) 
  }
}

getReactomePthLevelDF <- function(){
  
  #' Reads and processes Reactome pathway data to extract pathway-level associations.
  #'
  #' This function reads a tab-delimited file containing human pathway data from Reactome,
  #' filters pathways that start with "R-HSA-", and selects relevant columns.
  #'
  #' @return A data frame with two columns:
  #'   \itemize{
  #'     \item PATHID: Reactome pathway identifier (e.g., "R-HSA-XXXXXX").
  #'     \item LEVELID: Corresponding top-level pathway identifier.
  #'   }
  
  cat("---- Remember to check if Reactome has newer versions of the pathway level files!\n", sep = "")
  
  # Check if newer version: https://reactome.org/download-data
  pth_level_ids <- read.table("~/Data/Reactome_Files/Complex_2_Pathway_human.txt", header = TRUE)
  
  pth_level_ids <- pth_level_ids %>% 
    dplyr::select(pathway, top_level_pathway) %>%
    dplyr::filter(str_starts(pathway, "R-HSA-"))
  
  colnames(pth_level_ids) <- c("PATHID", "LEVELID")
  
  return(pth_level_ids)
}

addManualPathwayLevels <- function(df) {
  
  # Check and install openxlsx if needed
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    cat("Installing required package 'openxlsx'...\n")
    install.packages("openxlsx")
  }
  
  # Load the package
  library(openxlsx)
  
  manual <- read.xlsx("~/Data/Reactome_Files/manual_categories_excel.xlsx", sheet = "all", colNames = TRUE) %>%
    dplyr::filter(pathway_level != "NONE") %>%
    dplyr::select(pathway, pathway_level)
  manual$pathway <- gsub("\u00A0", " ", manual$pathway)
  
  result_df <- df %>% left_join(manual, by = "pathway", relationship = "many-to-many") 
  
  result_df$pathway_level <- ifelse(is.na(result_df$pathway_level.x), 
                                    result_df$pathway_level.y, 
                                    result_df$pathway_level.x)
  result_df <- result_df %>% dplyr::select(-pathway_level.x, -pathway_level.y)
  
  return(result_df)
}

getReactome <- function(de.results.df, rank_metric = NULL){
  
  #' Generate Reactome Pathway Data for Differential Gene Expression Analysis
  #'
  #' This function processes differential gene expression (DGE) results to generate 
  #' several data frames related to Reactome pathway analysis. It mimics the `reactomePathways()` 
  #' function and adds additional lines to incorporate pathway levels. It prepares the data 
  #' for downstream pathway analysis using `fgsea` and generates a list of pathways with corresponding 
  #' gene sets and statistical results.
  #'
  #' @param de.results.df A data frame containing the differential gene expression results. 
  #'   It should include columns for gene identifiers, p-values, fold changes, and other statistics 
  #'   relevant to DGE analysis.
  #'
  #' @return A list containing several data frames:
  #' \item{gene_list}{A sorted list of genes based on ranks (derived from p-values and the direction of differential expression).}
  #' \item{final.React}{A list of genes for `fgsea` with Reactome pathways, ordered from upregulated to downregulated genes.}
  #' \item{fgsea.set.React}{A list of pathways and their corresponding leading edge genes.}
  #' \item{fgsea.res.React}{A data frame of Reactome pathway results, including pathway names, p-values, NES values, and leading edge genes.}
  #' \item{fgsea.res.tidy.React}{A tidied version of the Reactome pathway results, sorted by NES, excluding pathways without any genes.}
  #' \item{collapsed.Pathways.React}{A list of collapsed Reactome pathways, aggregating related pathways.}
  #' \item{main.pathways}{A list of main pathways that are included in the collapsed pathways.}
  
  source("~/Data/Reactome_Files/reactome_variables.R")
  
  cat("- 1/11 => Check if Reactome and AnnotationDbi packages are loaded\n", sep = "")
  
  stopifnot(requireNamespace("reactome.db"))
  stopifnot(requireNamespace("AnnotationDbi"))
  
  # tic("Pathway analysis - Reactome DB")
  
  all_react_list <- list()
  main_react_list <- list()
  parent_react_list <- list()
  
  cat("- 2/11 => Load the DB and map out genes to the ENTREZID codes\n", sep = "")
  
  entrez.db <- org.Hs.eg.db
  g.entrezid <- mapIds(entrez.db,
                       keys = de.results.df$gene_name,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
  
  de.results.df$feature <- g.entrezid
  genes <- unique(g.entrezid)
  
  cat("- 3/11 => Calculate pathway levels\n", sep = "")
  
  ### The following code substitutes reactomePathways(unique(g.entrezid))
  ### fgsea.set.React <- reactomePathways(unique(g.entrezid))
  
  ### ****************************************************************************************************
  pathways_fun <- na.omit(AnnotationDbi::select(reactome.db::reactome.db, keys = genes, c("PATHID"), keytype = "ENTREZID"))
  pathways_fun <- split(pathways_fun$ENTREZID, pathways_fun$PATHID)
  pathway2name <- as.data.table(AnnotationDbi::select(reactome.db::reactome.db, names(pathways_fun), c("PATHNAME"), "PATHID"))
  
  PATHID <- NULL
  pathway2name <- pathway2name[!duplicated(PATHID)]
  
  PATHNAME <- NULL
  pathway2name[, `:=`(PATHNAME, sub("^[^:]*: ", "", PATHNAME))]
  ## There are duplicates in pathway2name: same PATHNAME with different PATHID
  
  # Add pathway category
  pth_level_ids <- getReactomePthLevelDF()
  pathway2name_with_level <- unique(merge(pathway2name, pth_level_ids, by = "PATHID", all.x = TRUE))
  ## Consequently, there are duplicates in pathway2name_with_level: same PATHNAME with different PATHID and hence, different LEVELID
  
  pathway2name_with_level$LEVELID <- ifelse(
    is.na(pathway2name_with_level$LEVELID) & pathway2name_with_level$PATHID %in% pth_level_ids$LEVELID,
    pathway2name_with_level$PATHID, 
    pathway2name_with_level$LEVELID
  )
  ## Consequently, there are duplicates in pathway2name_with_level: same PATHNAME with different PATHID and hence, different LEVELID
  
  pathway2name_with_level$LEVELNAME <- levelname_assignments[pathway2name_with_level$LEVELID]
  pathway2name_with_level$PATHNAME <- gsub("\u00A0", " ", pathway2name_with_level$PATHNAME)
  
  ## END ##
  
  name2pathways <- split(pathway2name$PATHID, pathway2name$PATHNAME)
  pathways_fun <- lapply(name2pathways, function(x) unique(do.call(c, pathways_fun[x])))
  fgsea.set.React <- pathways_fun[!is.na(names(pathways_fun))]
  ## ****************************************************************************************************
  
  cat("- 4/11 => Generate the Ranked Gene List\n", sep = "")
  # create a sorted list based on gene ranks
  gene_list <- de.results.df
  if(is.null(rank_metric)) {
    gene_list$PValue <- ifelse(gene_list$PValue == 0, 1e-100, gene_list$PValue)
    gene_list$fcsign <- sign(gene_list$logFC)
    gene_list$logP <- -log10(gene_list$PValue)
    # gene_list$metric <- gene_list$logP/gene_list$fcsign
    # multiplying instead of dividing does not give error in case logFC = 0 => fcsign = 0 (0 in denominator)
    gene_list$metric <- gene_list$logP * gene_list$fcsign
  } else {
    colnames(gene_list)[colnames(gene_list) == rank_metric] <- "metric"
    gene_list <- gene_list %>% arrange(desc(metric))
  }
  
  cat("- 5/11 => Prepare list of genes for fgsea with Reactome\n", sep = "")
  final.React <- gene_list[,c("feature", "metric")]
  final.React <- na.omit(final.React[order(final.React$metric),])
  final.React <- deframe(final.React)
  
  cat("- 6/11 => Perform pathway analysis -> fgsea.res.React\n", sep = "")
  suppressWarnings(
    fgsea.res.React <- fgsea(pathways=fgsea.set.React,
                             stats=final.React,
                             minSize=1, 
                             maxSize = Inf,
                             nproc=1)
  )
  
  cat("- 7/11 => Sort, select desired columns and filter by leadingEdge>0 -> fgsea.res.tidy.React.filtered\n", sep = "")
  fgsea.res.tidy.React <- fgsea.res.React %>% arrange(desc(NES)) %>% dplyr::select(-ES, -log2err)
  fgsea.res.tidy.React.filtered <- fgsea.res.tidy.React[ sapply(fgsea.res.tidy.React$leadingEdge, length) > 0 ]
  
  nrows_filtered <- nrow(fgsea.res.tidy.React) - nrow(fgsea.res.tidy.React.filtered)
  
  cat("- 8/11 => Map leadingEdge genes back to gene names -> fgsea.res.tidy.React.filtered\n", sep = "")
  fgsea.res.tidy.React.filtered$leadingEdge <- mapIdsList(entrez.db,
                                                          keys=fgsea.res.tidy.React.filtered$leadingEdge,
                                                          keytype="ENTREZID",
                                                          column="SYMBOL")
  
  fgsea.res.tidy.React.filtered$leadingEdge <- as.character(lapply(fgsea.res.tidy.React.filtered$leadingEdge, toString))
  
  cat("- 9/11 => Find collapsed pathways: main and parent -> collapsed.Pathways.React\n", sep = "")
  suppressWarnings(
    collapsed.Pathways.React <- collapsePathways(fgsea.res.React[order(pval)][pval < 0.01],
                                                 fgsea.set.React, final.React)
  )
  
  main.pathways <- fgsea.res.tidy.React.filtered[fgsea.res.tidy.React.filtered$pathway %in% collapsed.Pathways.React$mainPathways,]
  parent.pathways <- fgsea.res.tidy.React.filtered[fgsea.res.tidy.React.filtered$pathway %in% collapsed.Pathways.React$parentPathways,]
  
  cat("- 10/11 => Add pathway levels\n", sep = "")
  
  # Before the merge operations at the end, add:
  fgsea.res.tidy.React$pathway <- gsub("\u00A0", " ", fgsea.res.tidy.React$pathway)
  main.pathways$pathway <- gsub("\u00A0", " ", main.pathways$pathway)
  parent.pathways$pathway <- gsub("\u00A0", " ", parent.pathways$pathway)
  
  # Merge
  fgsea.res.tidy.React <- merge(fgsea.res.tidy.React, unique(pathway2name_with_level[, .(PATHNAME, LEVELNAME)]), by.x = "pathway", by.y = "PATHNAME", all.x = TRUE)
  names(fgsea.res.tidy.React)[names(fgsea.res.tidy.React) == "LEVELNAME"] <- "pathway_level"
  
  main.pathways <- merge(main.pathways, unique(pathway2name_with_level[, .(PATHNAME, LEVELNAME)]), by.x = "pathway", by.y = "PATHNAME", all.x = TRUE)
  names(main.pathways)[names(main.pathways) == "LEVELNAME"] <- "pathway_level"
  
  parent.pathways <- merge(parent.pathways, unique(pathway2name_with_level[, .(PATHNAME, LEVELNAME)]), by.x = "pathway", by.y = "PATHNAME", all.x = TRUE)
  names(parent.pathways)[names(parent.pathways) == "LEVELNAME"] <- "pathway_level"
  
  cat("- 11/11 => Add manual pathway levels in the Excel file\n", sep = "")
  fgsea.res.tidy.React <- addManualPathwayLevels(fgsea.res.tidy.React)
  main.pathways <- addManualPathwayLevels(main.pathways)
  parent.pathways <- addManualPathwayLevels(parent.pathways)
  
  # Save into list
  all_react_list[["Reactome"]] <- fgsea.res.tidy.React
  main_react_list[["Reactome"]] <- main.pathways
  parent_react_list[["Reactome"]] <- parent.pathways
  
  # toc()
  
  return(list(all_react_list = all_react_list, 
              main_react_list = main_react_list, 
              parent_react_list = parent_react_list))
}

getMSigDbPathways <- function(de.results.df, database, rank_metric = NULL){
  #' Generate MSigDB Pathway Data for Differential Gene Expression Analysis
  #'
  #' This function processes differential gene expression (DGE) results to generate 
  #' several data frames related to MSigDB pathway analysis. It prepares the data for 
  #' downstream pathway analysis using `fgsea`, generating a list of pathways with 
  #' corresponding gene sets and statistical results.
  #'
  #' @param de.results.df A data frame containing the differential gene expression results. 
  #'   It should include columns for gene identifiers, p-values, fold changes, and other statistics 
  #'   relevant to DGE analysis.
  #'
  #' @return A list containing several data frames:
  #' \item{gene_list}{A sorted list of genes based on ranks (derived from p-values and the direction of differential expression).}
  #' \item{final}{A list of genes for `fgsea` with MSigDB pathways, ordered from upregulated to downregulated genes.}
  #' \item{fgsea.set}{A list of pathways and their corresponding leading edge genes.}
  #' \item{fgsea.res}{A data frame of MSigDB pathway results, including pathway names, p-values, NES values, and leading edge genes.}
  #' \item{fgsea.res.tidy}{A tidied version of the MSigDB pathway results, sorted by NES, excluding pathways without any genes.}
  
  # tic("Pathway analysis - MSigDB C2 and Hallmarks")
  
  all_msigdb_list <- list()
  main_msigdb_list <- list()
  parent_msigdb_list <- list()
  
  cat("- 1/6 => Load the DB and map out genes to the ENTREZID codes\n", sep = "")
  symbol.db <- EnsDb.Hsapiens.v86
  msig.df <- list()
  
  if(database == "C2") {
    msig.df[["C2"]] <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP")
  } else if(database == "H") {
    msig.df[["H"]] <- msigdbr(species = "Homo sapiens", collection = "H")
  }
  
  # fgsea.set <- msig.df[[database]] %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea.set <- split(x = msig.df[[database]]$gene_symbol, f = msig.df[[database]]$gs_name)
  
  cat("- 2/6 => Generate the Ranked Gene List\n", sep = "")
  gene_list <- de.results.df
  if(is.null(rank_metric)) {
    gene_list$PValue <- ifelse(gene_list$PValue == 0, 1e-100, gene_list$PValue)
    gene_list$fcsign <- sign(gene_list$logFC)
    gene_list$logP <- -log10(gene_list$PValue)
    gene_list$metric <- gene_list$logP/gene_list$fcsign
  } else {
    colnames(gene_list)[colnames(gene_list) == rank_metric] <- "metric"
    gene_list <- gene_list %>% arrange(desc(metric))
  }
  
  cat("- 3/6 => Prepare list of genes for fgsea with MSigDB\n", sep = "")
  final <- gene_list[,c("gene_name", "metric")]
  final <- na.omit(final[order(final$metric),])
  final <- deframe(final)
  
  cat("- 4/6 => Perform pathway analysis -> fgsea.res\n", sep = "")
  suppressWarnings(
    fgsea.res <- fgsea(pathways=fgsea.set,
                       stats=final,
                       minSize=1, 
                       maxSize = Inf,
                       nproc=1)
  )
  
  cat("- 5/6 => Sort and select desired columns -> fgsea.res.tidy\n", sep = "")
  fgsea.res.tidy <- fgsea.res %>% as_tibble() %>% arrange(desc(NES)) %>% dplyr::select(-leadingEdge, -ES, -log2err)
  
  cat("- 6/6 => Find collapsed pathways: main and parent -> collapsed.Pathways\n", sep = "")
  suppressWarnings(
    collapsed.Pathways <- collapsePathways(fgsea.res[order(pval)][pval < 0.01], fgsea.set, final)
  )
  
  main.pathways <- fgsea.res.tidy[fgsea.res.tidy$pathway %in% collapsed.Pathways$mainPathways,]
  parent.pathways <- fgsea.res.tidy[fgsea.res.tidy$pathway %in% collapsed.Pathways$parentPathways,]
  
  all_msigdb_list[[database]] <- fgsea.res.tidy
  main_msigdb_list[[database]] <- main.pathways
  parent_msigdb_list[[database]] <- parent.pathways
  
  # toc()
  
  return(list(all_msigdb_list = all_msigdb_list, 
              main_msigdb_list = main_msigdb_list, 
              parent_msigdb_list = parent_msigdb_list))
}

getImportantPathways <- function(pathways) {
  
  #' Get Important Pathways from MSigDB or Reactome
  #'
  #' This function filters and extracts important pathways from a given data frame of pathways. 
  #' It identifies pathways related to specific biological processes such as apoptosis, inflammation, 
  #' immune response, and more. The function applies case-insensitive regular expression matching 
  #' to identify relevant pathways based on specific keywords.
  #'
  #' The function uses `dplyr` functions for efficient filtering and combining pathway data.
  #'
  #' @param pathways A data frame containing pathway information. It must have at least two columns: 
  #'   - `pathway`: A character vector with pathway names or descriptions.
  #'
  #' @return A data frame containing the filtered pathways related to the following categories:
  #'   - Apoptosis and anti-apoptosis
  #'   - Inflammatory pathways (excluding anti-inflammatory)
  #'   - Anti-inflammatory pathways
  #'   - Immune pathways (e.g., T-cell, B-cell)
  #'   - Cytokine pathways
  #'   - Defensin pathways
  #'   - Interleukin pathways
  #'   - Interferon pathways
  #'   - TNF pathways
  #'   - COX pathways
  #'   - IDO pathways
  #'   - MHC pathways
  #'   - Extracellular matrix pathways
  #'   - Collagen pathways
  #'   - Laminin pathways
  #'   - Glycoprotein and proteoglycan pathways
  #'   - Stress response pathways (e.g., UPR, ER stress)
  #'   - Proliferation-related pathways (e.g., cell cycle, cyclins, CDKs)
  #'   - Growth-related pathways (e.g., RTK, growth factors)
  #'   - DNA repair and unfolded protein response pathways
  
  # Helper function to filter pathways using case-insensitive pattern matching
  filter_pathways <- function(pattern, pathways) {
    pathways %>%
      filter(grepl(pattern, pathway, ignore.case = TRUE))
  }
  
  # APOPTOSIS (and antiapoptosis) PATHWAY
  apoptosis_PA <- bind_rows(
    filter_pathways("poptosis", pathways),
    filter_pathways("death", pathways)
  )
  
  # INFLAMMATORY PATHWAY sign (excluding anti-inflammatory)
  inflammatory_PA <- pathways %>%
    filter(grepl("inflamm", pathway, ignore.case = TRUE)) %>%
    filter(!grepl("anti", pathway, ignore.case = TRUE))
  
  # ANTI-INFLAMMATORY PATHWAY sign
  anti_inflammatory_PA <- filter_pathways("nti-inflammatory", pathways)
  
  # IMMUNE PATHWAY sign
  immune_PA <- bind_rows(
    filter_pathways("mmune", pathways),
    filter_pathways("lymph", pathways),
    filter_pathways("b_cell", pathways),
    filter_pathways("t_cell", pathways)
  )
  
  # CYTOKINE PATHWAY sign
  cytokine_PA <- filter_pathways("ytokin", pathways)
  
  # DEFENSINS PATHWAY sign
  defensins_PA <- filter_pathways("efensins", pathways)
  
  # INTERLEUKINS PATHWAY sign
  interleukin_PA <- bind_rows(
    filter_pathways("nterleukin", pathways),
    filter_pathways("IL", pathways)
  )
  
  # INTERFERON PATHWAY sign
  interferon_PA <- bind_rows(
    filter_pathways("nterferon", pathways),
    filter_pathways("IFN", pathways)
  )
  
  # TNF PATHWAY sign
  tnf_PA <- filter_pathways("tnf", pathways)
  
  # COX PATHWAY sign
  cox_PA <- bind_rows(
    filter_pathways("COX", pathways),
    filter_pathways("cyclooxygenase", pathways),
    filter_pathways("prostagland", pathways)
  )
  
  # IDO PATHWAY sign
  ido_PA <- bind_rows(
    filter_pathways("IDO", pathways),
    filter_pathways("indoleamine", pathways)
  )
  
  # MHC PATHWAY sign
  mhc_PA <- filter_pathways("MHC", pathways)
  
  # EXTRACELLULAR PATHWAY sign
  extracellular_PA <- bind_rows(
    filter_pathways("xtracellular", pathways),
    filter_pathways("ECM", pathways),
    filter_pathways("tissue", pathways),
    filter_pathways("regeneration", pathways),
    filter_pathways("metallo", pathways),
    filter_pathways("membrane", pathways),
    filter_pathways("matrisome", pathways)
  )
  
  # COLLAGEN, LAMININ, and GLYCOPROTEINS PATHWAY signs
  collagen_PA <- filter_pathways("ollagen", pathways)
  laminin_PA <- filter_pathways("aminin", pathways)
  glycoprot_PA <- bind_rows(
    filter_pathways("glycoprot", pathways),
    filter_pathways("proteoglyc", pathways)
  )
  
  # STRESS, PROLIFERATION, GROWTH, CELL RESPONSE PATHWAY signs
  stress_PA <- bind_rows(
    filter_pathways("tress", pathways),
    filter_pathways("UPR", pathways),
    filter_pathways("ER", pathways)
  )
  
  proliferation_PA <- bind_rows(
    filter_pathways("cycle", pathways),
    filter_pathways("cyclins", pathways),
    filter_pathways("CDK", pathways)
  )
  
  growth_PA <- bind_rows(
    filter_pathways("RTK", pathways),
    filter_pathways("growth", pathways),
    filter_pathways("factor", pathways)
  )
  
  cellresponses_PA <- filter_pathways("CD", pathways)
  dna_repain_PA <- filter_pathways("repair", pathways)
  unfolded_prot_PA <- filter_pathways("unfolded", pathways)
  
  # Combine all pathways into one data frame
  PATHWAYS <- bind_rows(
    apoptosis_PA,
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
    unfolded_prot_PA
  )
  
  return(PATHWAYS)
}

addPathwayCategory <- function(pathways){
  
  #' Add Pathway Category Based on Keywords in Pathway Names
  #'
  #' This function adds a new column to the input pathways data frame, categorizing each pathway
  #' based on specific keywords found in the pathway name. The categories are determined using `grep` 
  #' with predefined patterns that are matched against pathway names (e.g., "apoptosis", "inflammatory", etc.).
  #' The resulting data frame will include an additional column with the corresponding pathway category.
  #'
  #' @param pathways A data frame containing pathway information. It must have at least one column:
  #'   - `pathway`: A character vector with pathway names or descriptions.
  #'
  #' @return A data frame containing the original pathways along with an additional column called `category`,
  #'   which is a categorical label assigned based on the pathway name. The categories are determined 
  #'   by matching keywords related to various biological processes (e.g., apoptosis, inflammation, immune response, etc.).
  
  pathways <- pathways %>%
  mutate(pathway_category = case_when(
    # Disease
    grepl("SARS-CoV|HIV|APOBEC3G|Leishmania|Disease|Viral|HCMV|Defective|Disorder|Drug resistance|Mutants|Cancer|Syndrome|Deficiency|Infection|Resistance|Impaired|Loss of function|Defects|Abnormalities|Congenital", pathway, ignore.case = TRUE) ~ "Disease",

    # Immune System
    grepl("immune|antigen|B Cell|BCR|TCR|ER-Phagosome|NK|TLR3|CLEC7A|antibody|NOD1/2 Signaling Pathway|NF-kappa-B|Growth hormone receptor signaling|PKR|interleukin|IFN|interferon|TNF|IL|Inflammasomes|Neutrophil|Receptors|Cascade", pathway, ignore.case = TRUE) ~ "Immune System",

    # Developmental Biology
    grepl("Developmental Biology|Cell Lineages|Gastrulation|Transition|Melanocyte", pathway, ignore.case = TRUE) ~ "Developmental Biology",

    # Signal Transduction
    grepl("Signal Transduction|RHO|RND2|NF-kB is activated|GTPase cycle|TGF-beta|relaxin|ERBB4|GPCR|MAPK|Phosphorylation|GTPase", pathway, ignore.case = TRUE) ~ "Signal Transduction",

    # Metabolism of proteins
    grepl("glycosilation|deubiquitination|GTP hydrolysis|Translation|codon|Protein|Folding|Ubiquitination|SUMO|N-glycan|O-glycosylation|Processing", pathway, ignore.case = TRUE) ~ "Metabolism of proteins",

    # Vesicle-mediated transport
    grepl("vesicle|Vesicle-mediated transport|Golgi|ER|Trafficking", pathway, ignore.case = TRUE) ~ "Vesicle-mediated transport",

    # Metabolism
    grepl("Nicotina|fructuose|fatty acids|SREBP|chondroitin|CoA|Biosynthesis|Catabolism|Synthesis|Acids|Glucose|Glycogen|Lipids|Nucleotide|Vitamin|Oxidations", pathway, ignore.case = TRUE) ~ "Metabolism",

    # Cell Cycle
    grepl("cell cycle|mitosis|G1|G2|mitotic|Synthesis of DNA|meio|checkpoint|chromosome|Telomeres|Anaphase|Metaphase|Telophase|Cytokinesis", pathway, ignore.case = TRUE) ~ "Cell Cycle",

    # Chromatin organization
    grepl("chromatin", pathway, ignore.case = TRUE) ~ "Chromatin organization",

    # DNA Repair
    grepl("repair|DNA Damage|Excision Repair|Mismatch Repair|Break|Homologous Recombination|Resolution", pathway, ignore.case = TRUE) ~ "DNA Repair",

    # Cellular responses to stimuli
    grepl("cellular response|Metallothioneins|(NLR) signaling pathways|chaperone|senescence|Response|Stress", pathway, ignore.case = TRUE) ~ "Cellular responses to stimuli",

    # Metabolism of RNA
    grepl("rRNA|mRNA|tRNA|RNA|Metabolism of RNA|Decay|Splicing|Editing|Transport", pathway, ignore.case = TRUE) ~ "Metabolism of RNA",

    # Neuronal System
    grepl("neural|Neuronal System|Channels|Synapses|Transmission|Potassium|Acetylcholine|GABA", pathway, ignore.case = TRUE) ~ "Neuronal System",

    # Gene expression (Transcription)
    grepl("RNA Polymerase|DNA methylation|RUNX|Transcript|Gene expression|Transcription|Epigenetic|Activates", pathway, ignore.case = TRUE) ~ "Gene expression (Transcription)",

    # Transport of small molecules
    grepl("plasma|Transport of small molecules|Transporters|Uptake|Exchange|Lipoprotein|Ion|Sodium|Zinc", pathway, ignore.case = TRUE) ~ "Transport of small molecules",

    # Extracellular matrix organization
    grepl("ECM|membrane|extracellular|glycoprotein|proteoglycan|collagen|ROBO|kerat|cell surface|cell junction|syndecan|gap junction|Extracellular matrix", pathway, ignore.case = TRUE) ~ "Extracellular matrix organization",

    # Autophagy
    grepl("autophagy|mitophagy|selective", pathway, ignore.case = TRUE) ~ "Autophagy",

    # Hemostasis
    grepl("platelet|PECAM1|Clot|Hemostasis|Homeostasis", pathway, ignore.case = TRUE) ~ "Hemostasis",

    # Digestion and absorption
    grepl("Digestion and absorption|Absorption|Intestinal", pathway, ignore.case = TRUE) ~ "Digestion and absorption",

    # Drug ADME
    grepl("ADME|Drug", pathway, ignore.case = TRUE) ~ "Drug ADME",

    # Programmed Cell Death
    grepl("Programmed Cell Death|Cell Death|Apoptosis|Necrosis", pathway, ignore.case = TRUE) ~ "Programmed Cell Death",

    # Sensory Perception
    grepl("Sensory Perception|Perception|Taste|Sound|Visual|Phototransduction", pathway, ignore.case = TRUE) ~ "Sensory Perception",

    # Circadian Clock
    grepl("Circadian Clock|Circadian", pathway, ignore.case = TRUE) ~ "Circadian Clock",

    # Organelle biogenesis and maintenance
    grepl("Biogenesis|Maintenance", pathway, ignore.case = TRUE) ~ "Organelle biogenesis and maintenance",

    # Reproduction
    grepl("Reproduction|Fertilization", pathway, ignore.case = TRUE) ~ "Reproduction",

    # Muscle contraction
    grepl("Contraction|Cardiac", pathway, ignore.case = TRUE) ~ "Muscle contraction",

    # Cell-Cell communication
    grepl("Cell-Cell communication|Junction|Adhesion|Regulation", pathway, ignore.case = TRUE) ~ "Cell-Cell communication",

    TRUE ~ "Other"
  ))

  pathways <- pathways %>%
    mutate(pathway_subcategory = case_when(
      grepl("Inflam", pathway, ignore.case = TRUE) ~ "Inflammatory Response",
      grepl("Innate", pathway, ignore.case = TRUE) ~ "Innate Immunity",
      grepl("interleukin", pathway, ignore.case = TRUE) ~ "IL",
      grepl("IL", pathway, ignore.case = FALSE) ~ "IL",
      grepl("IFN|interferon", pathway, ignore.case = TRUE) ~ "IFN",
      grepl("TNF", pathway, ignore.case = TRUE) ~ "TNF",
      grepl("Immun", pathway_category, ignore.case = TRUE) ~ "Immune System",
      TRUE ~ "Other"
    ))
  
  return(pathways)
}

preparePathwayData <- function(pathways, group, n_top, n_tail, per_cluster){
  
  #' Prepare Pathway Data for Plotting
  #'
  #' This function prepares the pathways data for further analysis or plotting. It filters the input pathways data frame
  #' according to specific parameters such as the group(s), top `n` pathways, tail `n` pathways, and whether to separate 
  #' the data by clusters. Additionally, the function reorders the pathways by the sum of their NES values across all treatment 
  #' groups and returns a filtered pathways data frame that can be used for subsequent analysis or visualizations (e.g., heatmaps).
  #'
  #' @param pathways A data frame containing pathway information with columns:
  #'   - `Group`: A factor or character vector with group names or treatment labels.
  #'   - `NES`: Numeric values representing the normalized enrichment score for each pathway.
  #'   - `pathway`: A character vector of pathway names.
  #' @param group A character vector of group names that you wish to filter the pathways by. If `per_cluster` is `TRUE`, 
  #'   pathways will be filtered by group.
  #' @param n_top An integer specifying the number of top pathways to retain for each group based on the NES score.
  #' @param n_tail An integer specifying the number of bottom pathways to retain for each group based on the NES score.
  #' @param per_cluster A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will filter and process pathways per 
  #'   each cluster/group. If `FALSE`, all pathways are processed together without cluster grouping.
  #'
  #' @return A data frame containing the filtered pathways data with the following columns:
  #'   - `pathway`: Pathway name.
  #'   - `NES`: Normalized enrichment score for the pathway.
  #'   - `Group`: The treatment or group associated with the pathway.
  #'   - `order`: A new column providing the rank order of pathways.
  #'   - `sum_nes`: The sum of NES values for each pathway across the groups (used for ordering pathways based on NES).
  #' 
  #' The data is filtered to retain the top `n_top` pathways with the highest NES values and the bottom `n_tail` pathways,
  #' and then it is ordered based on the sum of NES values across all treatment groups (useful for plotting or heatmap visualization).
  
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

segregateIntoSeveralColumns <- function(pathways, facetParam, n_columns){
  
  #' Segregate Pathways into Multiple Columns for Faceting
  #'
  #' This function takes the pathways data frame and creates an additional column called `Column` to facilitate
  #' the division of a plot into multiple columns for visualization using `facet_grid()` in `ggplot2`. The `Column`
  #' is calculated by dividing the values of the specified `facetParam` into `n_columns` groups in such a way as to
  #' balance the frequencies of the groups, minimizing the difference in size between columns.
  #'
  #' @param pathways A data frame containing pathways information, which includes the `facetParam` column.
  #' @param facetParam A string specifying the column name in the `pathways` data frame that will be used to divide the data into columns.
  #' @param n_columns An integer specifying the number of columns to divide the data into for faceting.
  #' 
  #' @return A data frame that includes all the original columns of the `pathways` data frame along with a new column `Column`,
  #'   which indicates the group assigned to each pathway based on the division of `facetParam` into `n_columns`.
  #'
  #' @details
  #' The function works by:
  #' 1. Determining the unique values of the `facetParam` and their frequencies.
  #' 2. Greedily assigning each value to one of the `n_columns` in a way that minimizes the imbalance in the group sizes.
  #' 3. Adding a `Column` column to the pathways data frame to indicate the group assignment.
  #'
  #' This is useful for visualizations where you want to split a large number of facet levels into several columns, 
  #' ensuring that the columns are balanced in size for better plot aesthetics.
  
  # 0. Delete the Column variable if already exists
  pathways <- pathways %>% dplyr::select(any_of(
    colnames(pathways)[!grepl("Column", colnames(pathways))]
  ))
  
  # 1. Select the unique values for facetParam field
  xi <- unique(pathways %>% dplyr::select(!!sym(facetParam)) %>% arrange(!!sym(facetParam)))
  
  # 2. Get frequencies for each value
  fi <- table(pathways[,facetParam])
  
  # 3. Remove values with freq = 0
  fi <- subset(table(pathways[,facetParam]), subset = fi > 0)
  fi <- as.data.frame(fi)
  fi <- fi[order(-fi$Freq), ]  # Sort by frequency in descending order
  
  # 4. Initialize the sum for each group
  sum_group <- numeric(n_columns)  # Initialize a vector for the sums of n_columns
  
  # 5. Greedily assign elements to the groups - add column 'group' to the fi DF
  for (i in 1:nrow(fi)) {
    # Find the group with the least sum
    group_index <- which.min(sum_group)
    
    # Assign the current element to this group
    sum_group[group_index] <- sum_group[group_index] + fi$Freq[i]
    fi$group[i] <- as.character(group_index)  # Assign the group index (1 to n_columns)
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

complete_DE_results <- function(current_de_results, group1, group2, subset_var = NULL, subset_value = NULL, cluster = NULL, clust_res = NULL) {
  
  ## Completes the results DF with metadata columns
  ## This is for any comparability so the subset parameter is passed in subset_var
  ## Function inputs:
  #' @param current_de_results
  #' @param group1 - metadata columns
  #' @param group2 - metadata columns
  #' @param subset_var - metadata columns
  #' @param subset_value - metadata columns
  #' @param cluster - metadata columns
  #' @param clust_res - metadata columns
  #' @return completed DF
  
  current_de_results <- current_de_results[, !(colnames(current_de_results) %in% c("logCPM", "F"))]
  
  current_de_results[["group1"]] <- group1
  current_de_results[["group2"]] <- group2
  current_de_results[subset_var] <- subset_value
  current_de_results[["group"]] <- paste0(current_de_results[["group2"]], "_vs_", current_de_results[["group1"]])
  
  if (!is.null(cluster) && !is.null(clust_res)) {
    current_de_results[["cluster"]] <- cluster
    current_de_results[["clust_res"]] <- clust_res
  }
  
  # current_de_results <- current_de_results[complete.cases(current_de_results), ]
  # Many pathways won't have pathway_level for the moment and that is OK
  
  if ("pathway_level" %in% names(current_de_results)) {
    cols_to_check <- setdiff(names(current_de_results), "pathway_level")
  } else {
    cols_to_check <- names(current_de_results)
  }
  
  current_de_results <- current_de_results[complete.cases(current_de_results[, cols_to_check, drop = FALSE]), ]
  
  return(current_de_results)
}

update_DE_results_df <- function(de_results_df, current_de_results) {
  
  #' Update Global Differential Expression Results
  #'
  #' This function appends the results from a new differential expression analysis (`current_de_results`) 
  #' to an existing global results data frame (`de_results_df`). If the global results data frame 
  #' is empty, it initializes the global results with the `current_de_results`.
  #'
  #' @param de_results_df A data frame containing the existing global differential expression results.
  #'   This data frame will be updated by appending `current_de_results` to it.
  #' @param current_de_results A data frame containing the results of the current differential expression analysis
  #'   that will be added to the global results.
  #' 
  #' @return A data frame representing the updated global results, containing both the old and the newly added results.
  #' 
  #' @details
  #' The function checks if the `de_results_df` is empty (i.e., contains no rows). If it is empty, 
  #' it directly assigns the `current_de_results` as the global results. Otherwise, it appends the `current_de_results` 
  #' to the existing global results, preserving the existing data and adding new rows.
  #'
  #' This function is useful when accumulating results from multiple analyses into a single data frame for further processing or visualization.
  
  if (nrow(de_results_df) == 0) {
    de_results_df <- current_de_results
  } else {
    de_results_df <- rbind(de_results_df, current_de_results)
  }
  return(de_results_df)
}

post_process_and_save_results <- function(de_results_df, params_for_post_process) {
  
  #' Post-process and Save Differential Expression Results
  #'
  #' This function processes and saves the global differential expression results (`de_results_df`) 
  #' according to the specified parameters in `params_for_post_process`. It handles column removal, empty data frame creation, 
  #' and saving the results to a file.
  #'
  #' @param de_results_df A data frame containing the global differential expression results that need to be processed and saved.
  #' @param params_for_post_process A named list containing the parameters required for post-processing and saving the results. 
  #'   It must include the following:
  #'   \itemize{
  #'     \item \strong{ncols_for_empty_df}: The number of columns to create in case the results data frame is empty.
  #'     \item \strong{colnames_for_empty_df}: A character vector specifying the column names for the empty data frame.
  #'     \item \strong{results_file_name}: The name (and path) of the file where the results will be saved.
  #'   }
  #' 
  #' @return The same `de_results_df`, potentially modified and saved to a file.
  #' 
  #' @details
  #' This function performs the following tasks:
  #' \itemize{
  #'   \item It removes unnecessary columns (`logCPM` and `F`) from the `de_results_df`.
  #'   \item It checks whether the data frame is empty. If it is, it creates a new empty data frame with the specified number of columns and column names.
  #'   \item It saves the results to a file using the file name provided in `params_for_post_process$results_file_name`.
  #' }
  #'
  #' If any of the required parameters (`ncols_for_empty_df`, `colnames_for_empty_df`, `results_file_name`) are missing from 
  #' the `params_for_post_process` list, an error is raised.
  
  if (!all(c("ncols_for_empty_df", "colnames_for_empty_df", "results_file_name") %in% names(params_for_post_process))) {
    rlang::abort("Please, provide all required parameters within the list: ncol_for_empty_df, colnames_for_empty_df, results_file_name.")
  }
  
  # Remove unnecessary columns
  de_results_df <- de_results_df[, !(colnames(de_results_df) %in% c("logCPM", "F"))]
  
  # Check if resulting dataframe is empty
  if (nrow(de_results_df) == 0) {
    de_results_df <- data.frame(matrix(vector(), ncol = params_for_post_process$ncol_for_empty_df))
    colnames(de_results_df) <- params_for_post_process$colnames_for_empty_df
    de_metrics_list <- list()
  }
  
  # save results to file
  write.table(de_results_df, file = params_for_post_process$results_file_name)
  return(de_results_df)
}
