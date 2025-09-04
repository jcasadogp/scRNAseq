myDimPlot <- function(obj, g, s = NULL, ncol = 1, colors=NULL, title, label = FALSE, legend = TRUE){
  
  #' Custom Dimensional Reduction Plot
  #'
  #' This function creates a customized UMAP plot using Seurat's DimPlot function. 
  #' It allows grouping by a specific metadata field, optional splitting into subplots, 
  #' custom color palettes, title setting, label display, and legend visibility control.
  #'
  #' @param obj A Seurat object containing dimensional reduction results.
  #' @param g A string indicating the metadata column to group cells by (used in `group.by`).
  #' @param s Optional. A string specifying a metadata field to split the plot by (used in `split.by`). Default is NULL.
  #' @param ncol An integer specifying the number of columns when splitting the plot. Default is 1.
  #' @param colors Optional. A vector of colors to manually specify the palette for groups. Default is NULL.
  #' @param title A string for the plot title.
  #' @param label Logical. If TRUE, labels for each group are added to the plot. Default is FALSE.
  #' @param legend Logical. If FALSE, the plot will be returned without a legend. Default is TRUE.
  #'
  #' @return A ggplot object representing the UMAP dimensional reduction plot with the specified customizations.
  
  plot <- DimPlot(obj, reduction = "umap", group.by = g, split.by = s, cols = colors, 
                  label = label, repel = label, 
                  raster = FALSE, shuffle = TRUE, ncol= ncol) +
    ggtitle(title) +
    theme(strip.text.x = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10))
  
  if (!legend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}


plot_percentage_cc <- function(obj, metadata_field, position = "fill", subtitle = NULL) {
  load("~/Data/color_vectors.Rdata")
  
  df <- as.data.frame(table(obj$Phase, obj[[metadata_field]][,1]))
  
  colnames(df) <- c("Phase", metadata_field, "Counts")
  df$Phase <- factor(df$Phase, levels = c("G1", "S", "G2M"))
  
  if(metadata_field == "mnf_step") {
    df$mnf_step <- factor(df$mnf_step, levels = c("SVF", "MCS", "FDS")) 
  }
  
  df <- df %>%
    group_by(!!sym(metadata_field)) %>%
    mutate(Percentage = (Counts / sum(Counts)) * 100)
  
  if(position == "fill") {
    
    plot <- ggplot(df, aes(x = !!sym(metadata_field), y = Counts, fill = Phase)) +
      geom_bar(stat = "identity", position = "fill") +
      geom_text(
        data = subset(df, Percentage >= 5), 
        aes(label = paste0(round(Percentage, 1), "%")), 
        position = position_fill(vjust = 0.5), 
        color = "white", size = 5)
    
  } else if(position == "stack") {
    
    plot <- ggplot(df, aes(x = factor(!!sym(metadata_field)), y = Counts, fill = Phase)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(data = subset(df, Percentage >= 5), 
                aes(label = paste0(round(Percentage, 1), "%")), 
                # stat = "identity", 
                position = position_stack(vjust = 0.5), 
                color = "white", size = 5)
    
  }
  
  plot <- plot +
    scale_fill_manual(values = phase_colors) +
    labs(x = NULL, y = NULL, fill = NULL, 
         title = "Percentage of cells on each Cell Cycle Phase", subtitle = subtitle) +
    theme_linedraw() +
    theme(legend.position = "top", 
          plot.title = element_text(face = "bold"))
  
  return(list(plot = plot, df = df))
  
}

plot_percentage_cluster <- function(df, metadata_field, position = "fill", subtitle = NULL) {
  
  if(position == "fill") {
    
    plot <- ggplot(df, aes(x = !!sym(metadata_field), y = Counts, fill = Phase)) +
      geom_bar(stat = "identity", position = "fill") +
      geom_text(
        data = subset(df, Percentage >= 5), 
        aes(label = paste0(round(Percentage, 1), "%")), 
        position = position_fill(vjust = 0.5), 
        color = "white", size = 5)
    
  } else if(position == "stack") {
    
    plot <- ggplot(df, aes(x = factor(!!sym(metadata_field)), y = Counts, fill = Phase)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(data = subset(df, Percentage >= 5), 
                aes(label = paste0(round(Percentage, 1), "%")), 
                # stat = "identity", 
                position = position_stack(vjust = 0.5), 
                color = "white", size = 5)
    
  }
  
  plot <- plot +
    labs(x = NULL, y = NULL, fill = NULL, 
         title = "Percentage of cells on each Cluster", subtitle = subtitle) +
    theme_linedraw() +
    theme(plot.title = element_text(face = "bold"))
  
  return(list(plot = plot, df = df))
  
}


create_sina_plot <- function(data, x_var = "group", y_var = "NES", fill_var = "group", facet_var = "pathway_level", ncol, title = NULL, colors = NULL) {
  
  p <- ggplot(data, aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(fill_var))) + 
    ggforce::geom_sina(alpha = 0.7, shape = 21, color = "black") +
    labs(x = NULL, y = y_var, fill = NULL, title = title) +
    geom_hline(yintercept = 0) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top", 
      legend.text = element_text(size = 12),
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank()
    )
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(vars(!!sym(facet_var)), ncol = ncol)
  }
  
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  
  return(p)
}
