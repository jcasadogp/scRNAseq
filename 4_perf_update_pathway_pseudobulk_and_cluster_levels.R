source("~/functions/functions_dge_pathways.R")

analysis_directory <- "velsera_object/"
pathway_results_directory <- paste0("~/Results/", analysis_directory, "pathway_results/")

# Load the pathway results dataframes
all_results_dataframe <- read.table(paste0(pathway_results_directory, "All_path_results_total_pseudobulk.tsv"))
all_results_dataframe_cluster <- read.table(paste0(pathway_results_directory, "All_path_results_per_cluster.tsv"))
main_results_dataframe <- read.table(paste0(pathway_results_directory, "Main_path_results_total_pseudobulk.tsv"))
main_results_dataframe_cluster <- read.table(paste0(pathway_results_directory, "Main_path_results_per_cluster.tsv"))
parent_results_dataframe <- read.table(paste0(pathway_results_directory, "Parent_path_results_total_pseudobulk.tsv"))
parent_results_dataframe_cluster <- read.table(paste0(pathway_results_directory, "Parent_path_results_per_cluster.tsv"))

all_results_dataframe <- addManualPathwayLevels(all_results_dataframe)
all_results_dataframe_cluster <- addManualPathwayLevels(all_results_dataframe_cluster)
main_results_dataframe <- addManualPathwayLevels(main_results_dataframe)
main_results_dataframe_cluster <- addManualPathwayLevels(main_results_dataframe_cluster)
parent_results_dataframe <- addManualPathwayLevels(parent_results_dataframe)
parent_results_dataframe_cluster <- addManualPathwayLevels(parent_results_dataframe_cluster)

# Save the pathway results dataframes with the pathway levels updated
write.table(all_results_dataframe, paste0(pathway_results_directory, "All_path_results_total_pseudobulk.tsv"))
write.table(all_results_dataframe_cluster, paste0(pathway_results_directory, "All_path_results_per_cluster.tsv"))
write.table(main_results_dataframe, paste0(pathway_results_directory, "Main_path_results_total_pseudobulk.tsv"))
write.table(main_results_dataframe_cluster, paste0(pathway_results_directory, "Main_path_results_per_cluster.tsv"))
write.table(parent_results_dataframe, paste0(pathway_results_directory, "Parent_path_results_total_pseudobulk.tsv"))
write.table(parent_results_dataframe_cluster, paste0(pathway_results_directory, "Parent_path_results_per_cluster.tsv"))


## Check if the manual categories excel should be updated ##
pathways_levels_with_na <- c(
  all_results_dataframe[all_results_dataframe$database == "Reactome" & is.na(all_results_dataframe$pathway_level), "pathway"],
  all_results_dataframe_cluster[all_results_dataframe_cluster$database == "Reactome" & is.na(all_results_dataframe_cluster$pathway_level), "pathway"],
  main_results_dataframe[main_results_dataframe$database == "Reactome" & is.na(main_results_dataframe$pathway_level), "pathway"],
  main_results_dataframe_cluster[main_results_dataframe_cluster$database == "Reactome" & is.na(main_results_dataframe_cluster$pathway_level), "pathway"],
  parent_results_dataframe[parent_results_dataframe$database == "Reactome" & is.na(parent_results_dataframe$pathway_level), "pathway"],
  parent_results_dataframe_cluster[parent_results_dataframe_cluster$database == "Reactome" & is.na(parent_results_dataframe_cluster$pathway_level), "pathway"]
)

pathways_levels_with_na <- unique(pathways_levels_with_na)

manual <- read.xlsx("~/Data/Reactome_Files/manual_categories_excel.xlsx", sheet = "all", colNames = TRUE) %>% dplyr::select(pathway, pathway_level)
manual$pathway <- gsub("\u00A0", " ", manual$pathway)

if(!all(pathways_levels_with_na %in% manual$pathway)) {
  cat("-- There are some pathways that are not in the manual categories excel => it should be updated!\n")
  
  pathways_levels_with_na[!pathways_levels_with_na %in% manual$pathway]
} else {
  cat("-- No need to updated the excel => good to go!\n")
}