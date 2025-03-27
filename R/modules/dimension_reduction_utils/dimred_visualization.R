# R/modules/dimension_reduction_utils/dimred_visualization.R

#' @title Extend Color Palette
#' @description Extends a color palette for potentially large number of clusters or categories.
#'   Uses RColorBrewer to generate colors but extends the palette for larger numbers.
#' @param n_colors Integer number of colors needed
#' @return Character vector of hex color codes
#' @keywords internal
extend_color_palette <- function(n_colors) {
  if (n_colors <= 8) {
    return(RColorBrewer::brewer.pal(8, "Set2")[1:n_colors])
  } else {
    return(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_colors))
  }
}

#' @title Save Elbow Plot to File
#' @description Saves a PCA elbow plot to file with a vertical line indicating the suggested
#'   number of dimensions to use.
#' @param file File path to save to
#' @param seurat_obj Seurat object with PCA
#' @param suggested_dims Integer suggested number of dimensions
#' @param width Integer plot width in pixels
#' @param height Integer plot height in pixels
#' @param dpi Integer plot resolution
#' @return Invisible NULL, called for side effect
#' @export
save_elbow_plot <- function(file, seurat_obj, suggested_dims, width = 3000, height = 1800, dpi = 300) {
  png(file, width = width, height = height, res = dpi)
  
  print(ElbowPlot(seurat_obj, 
                  ndims = ncol(Embeddings(seurat_obj, "pca"))) +
          geom_vline(xintercept = suggested_dims, 
                     color = "red", 
                     linetype = "dashed"))
  
  dev.off()
}

#' @title Get Available Coloring Options
#' @description Determines available options for coloring UMAPs based on the metadata
#'   columns in a Seurat object. Prioritizes sample, cluster, and condition columns.
#' @param seurat_obj Seurat object
#' @return Named character vector of coloring options
#' @export
get_coloring_options <- function(seurat_obj) {
  # Base options
  options <- c("cluster" = "cluster")
  
  # Add sample if available
  if (!is.null(seurat_obj) && "sample" %in% colnames(seurat_obj@meta.data)) {
    options <- c("sample" = "sample", options)
  }
  
  # Add condition options if available
  meta_cols <- colnames(seurat_obj@meta.data)
  condition_cols <- grep("condition|treatment|group|genotype|timepoint|characteristics", 
                         meta_cols, value = TRUE, ignore.case = TRUE)
  
  if (length(condition_cols) > 0) {
    condition_options <- setNames(condition_cols, condition_cols)
    options <- c(options, condition_options)
  }
  
  # Always add gene option
  options <- c(options, "gene" = "gene")
  
  return(options)
}

#' @title Search Genes in Seurat Object
#' @description Searches for genes in a Seurat object by gene symbol or Ensembl ID.
#'   Returns a data frame of matching genes for display.
#' @param seurat_obj Seurat object
#' @param query String gene name or ID query
#' @return Data frame of matching genes or NULL if none found
#' @export
search_genes <- function(seurat_obj, query) {
  if (is.null(seurat_obj) || trimws(query) == "") {
    return(NULL)
  }
  
  # Convert query to lowercase for case-insensitive search
  query <- tolower(trimws(query))
  
  # Get gene mapping
  gene_mapping <- seurat_obj@misc$gene_mapping
  
  # Initialize empty results
  results <- NULL
  
  if (!is.null(gene_mapping)) {
    # Create a data frame of gene mappings for easier searching
    gene_df <- data.frame(
      ensembl_id = names(gene_mapping),
      gene_symbol = as.character(gene_mapping),
      stringsAsFactors = FALSE
    )
    
    # Remove any rows with NA values
    gene_df <- gene_df[!is.na(gene_df$ensembl_id) & !is.na(gene_df$gene_symbol), ]
    
    # Search by gene symbol (prioritize this)
    symbol_matches <- grep(query, tolower(gene_df$gene_symbol), value = FALSE)
    
    # Search by Ensembl ID if needed
    ensembl_matches <- grep(query, tolower(gene_df$ensembl_id), value = FALSE)
    
    # Combine matches (prioritize gene symbol matches)
    all_matches <- unique(c(symbol_matches, ensembl_matches))
    
    if (length(all_matches) > 0) {
      # Create results data frame
      results <- gene_df[all_matches, , drop = FALSE]
      
      # Only keep genes that are in the dataset
      results <- results[results$ensembl_id %in% rownames(seurat_obj), , drop = FALSE]
      
      # Sort by gene symbol
      if (nrow(results) > 0) {
        results <- results[order(results$gene_symbol), , drop = FALSE]
      } else {
        results <- NULL
      }
    }
  } else {
    # If no gene mapping is available, search directly in rownames
    matches <- grep(query, tolower(rownames(seurat_obj)), value = TRUE)
    
    if (length(matches) > 0) {
      results <- data.frame(
        ensembl_id = matches,
        gene_symbol = matches,
        stringsAsFactors = FALSE
      )
      results <- results[order(results$gene_symbol), , drop = FALSE]
    }
  }
  
  return(results)
}

#' @title Create 2D UMAP Plot
#' @description Creates a standard 2D UMAP visualization with customizable coloring and filtering.
#'   Can color by cluster, sample, gene expression, or any metadata column.
#' @param seurat_obj Seurat object
#' @param color_by String attribute to color by ("cluster", "sample", "gene", or metadata column)
#' @param gene_id String gene ID (only used if color_by = "gene")
#' @param reduction String name of reduction to use
#' @param label Logical whether to show labels
#' @param pt_size Numeric point size
#' @param active_items Optional vector of items to include (for filtering)
#' @param ... Additional parameters to pass to DimPlot or FeaturePlot
#' @return ggplot object
#' @export
create_2d_umap_plot <- function(seurat_obj, color_by = "cluster", gene_id = NULL,
                                reduction = "umap2d", label = TRUE, pt_size = 1, 
                                active_items = NULL, ...) {
  # Handle different coloring options
  if (color_by == "cluster" && "seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    # Color by cluster with optional filtering
    if (!is.null(active_items) && length(active_items) > 0) {
      # Filter to only show active clusters
      cells_to_keep <- seurat_obj$seurat_clusters %in% active_items
      plot_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
    } else {
      plot_obj <- seurat_obj
    }
    
    return(DimPlot(plot_obj, 
                   reduction = reduction, 
                   label = label,
                   pt.size = pt_size,
                   ...))
    
  } else if (color_by == "sample" && "sample" %in% colnames(seurat_obj@meta.data)) {
    # Color by sample with optional filtering
    if (!is.null(active_items) && length(active_items) > 0) {
      # Filter to only show active samples
      cells_to_keep <- seurat_obj$sample %in% active_items
      plot_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
    } else {
      plot_obj <- seurat_obj
    }
    
    return(DimPlot(plot_obj, 
                   reduction = reduction, 
                   group.by = "sample",
                   pt.size = pt_size,
                   ...))
    
  } else if (color_by == "gene" && !is.null(gene_id) && gene_id %in% rownames(seurat_obj)) {
    # For gene expression, we use the original object but apply cell highlighting if needed
    return(FeaturePlot(seurat_obj, 
                       features = gene_id, 
                       reduction = reduction,
                       pt.size = pt_size,
                       ...))
    
  } else if (color_by %in% colnames(seurat_obj@meta.data)) {
    # Color by any metadata column with optional filtering
    if (!is.null(active_items) && length(active_items) > 0) {
      # Filter to only show active values for the selected column
      cells_to_keep <- seurat_obj@meta.data[[color_by]] %in% active_items
      plot_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
    } else {
      plot_obj <- seurat_obj
    }
    
    return(DimPlot(plot_obj, 
                   reduction = reduction, 
                   group.by = color_by,
                   pt.size = pt_size,
                   ...))
  } else {
    # Default empty plot with error message
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "Please select a valid coloring option") + 
             theme_void())
  }
}

#' @title Create 3D UMAP Plot
#' @description Creates an interactive 3D UMAP visualization using plotly.
#'   Can color by cluster, sample, or any metadata column.
#' @param seurat_obj Seurat object
#' @param color_by String attribute to color by ("cluster", "sample", or metadata column)
#' @param reduction String name of reduction to use
#' @param point_size Numeric point size
#' @param opacity Numeric opacity (0-1)
#' @param active_items Optional vector of items to include (for filtering)
#' @return plotly object
#' @export
create_3d_umap_plot <- function(seurat_obj, color_by = "cluster", reduction = "umap3d",
                                point_size = 3, opacity = 0.7, active_items = NULL) {
  # Handle subsetting based on active items
  if (!is.null(active_items) && length(active_items) > 0) {
    if (color_by == "cluster" && "seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      cells_to_keep <- seurat_obj$seurat_clusters %in% active_items
      plot_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
    } else if (color_by == "sample" && "sample" %in% colnames(seurat_obj@meta.data)) {
      cells_to_keep <- seurat_obj$sample %in% active_items
      plot_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
    } else if (color_by %in% colnames(seurat_obj@meta.data)) {
      cells_to_keep <- seurat_obj@meta.data[[color_by]] %in% active_items
      plot_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
    } else {
      plot_obj <- seurat_obj
    }
  } else {
    plot_obj <- seurat_obj
  }
  
  # Check if we have cells to plot
  if (ncol(plot_obj) == 0) {
    return(plotly_empty() %>% 
             layout(title = "No cells to display with current filter settings"))
  }
  
  # Extract UMAP coordinates
  umap_data <- Embeddings(plot_obj[[reduction]])
  
  # Handle different coloring options
  if (color_by == "cluster" && "seurat_clusters" %in% colnames(plot_obj@meta.data)) {
    color_var <- as.factor(plot_obj$seurat_clusters)
    plot_title <- "3D UMAP - Colored by Clusters"
    hover_text <- paste("Cluster:", color_var)
  } else if (color_by == "sample" && "sample" %in% colnames(plot_obj@meta.data)) {
    color_var <- as.factor(plot_obj$sample)
    plot_title <- "3D UMAP - Colored by Sample"
    hover_text <- paste("Sample:", color_var)
  } else if (color_by %in% colnames(plot_obj@meta.data)) {
    color_var <- as.factor(plot_obj@meta.data[[color_by]])
    plot_title <- paste0("3D UMAP - Colored by ", color_by)
    hover_text <- paste(color_by, ":", color_var)
  } else {
    # Default if no valid coloring option
    color_var <- rep("Sample", ncol(plot_obj))
    plot_title <- "3D UMAP"
    hover_text <- rep("Cell", length(color_var))
  }
  
  # Create custom color palette
  if (length(unique(color_var)) > 1) {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(color_var)))
  } else {
    colors <- "blue"
  }
  
  # Create 3D plot
  p <- plot_ly(x = umap_data[,1], 
               y = umap_data[,2], 
               z = umap_data[,3],
               type = "scatter3d",
               mode = "markers",
               color = color_var,
               colors = colors,
               text = hover_text,
               hoverinfo = "text",
               marker = list(size = point_size, opacity = opacity)) %>%
    layout(title = plot_title,
           scene = list(xaxis = list(title = "UMAP 1"),
                        yaxis = list(title = "UMAP 2"),
                        zaxis = list(title = "UMAP 3")))
  
  # Custom legend with larger markers
  p$x$data <- lapply(p$x$data, function(d) {
    if (!is.null(d$marker)) {
      d$showlegend <- TRUE
      d$marker$sizeref <- 0.2
    }
    return(d)
  })
  
  # Adjust legend position and size
  p <- p %>% layout(
    margin = list(r = 120),
    legend = list(
      font = list(size = 12),
      itemsizing = "constant",
      y = 0.5
    )
  )
  
  return(p)
}

#' @title Create 3D Gene Expression UMAP Plot
#' @description Creates an interactive 3D UMAP visualization colored by gene expression.
#' @param seurat_obj Seurat object
#' @param gene_id String gene ID
#' @param reduction String name of reduction to use
#' @param point_size Numeric point size
#' @param opacity Numeric opacity (0-1)
#' @return plotly object
#' @export
create_3d_gene_umap_plot <- function(seurat_obj, gene_id, reduction = "umap3d",
                                     point_size = 4, opacity = 0.7) {
  # Check if gene exists in the dataset
  if (!(gene_id %in% rownames(seurat_obj))) {
    return(plotly_empty() %>% 
             layout(title = paste("Gene", gene_id, "not found in dataset")))
  }
  
  # Get gene symbol if available
  gene_symbol <- if (!is.null(seurat_obj@misc$gene_mapping) && 
                     gene_id %in% names(seurat_obj@misc$gene_mapping) &&
                     !is.na(seurat_obj@misc$gene_mapping[gene_id])) {
    seurat_obj@misc$gene_mapping[gene_id]
  } else {
    gene_id
  }
  
  # Get expression values
  expr_values <- GetAssayData(seurat_obj, slot = "data")[gene_id, ]
  
  # Extract UMAP coordinates
  umap_data <- Embeddings(seurat_obj[[reduction]])
  
  # Extract cluster information for hover text
  clusters <- seurat_obj$seurat_clusters
  
  # Determine color scale for expression values
  # Normalize expression to a 0-1 scale for better color gradients
  expr_min <- min(expr_values)
  expr_max <- max(expr_values)
  
  # Create normalized expression values for coloring (avoid division by zero)
  if (expr_min == expr_max) {
    normalized_expr <- rep(0.5, length(expr_values))
  } else {
    normalized_expr <- (expr_values - expr_min) / (expr_max - expr_min)
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(
    umap1 = umap_data[, 1],
    umap2 = umap_data[, 2],
    umap3 = umap_data[, 3],
    expr = expr_values,
    norm_expr = normalized_expr,
    cluster = as.factor(clusters)
  )
  
  # Create 3D plot
  p <- plot_ly(
    data = plot_data,
    x = ~umap1, 
    y = ~umap2, 
    z = ~umap3,
    color = ~expr,
    colors = colorRamp(c("lightgrey", "red")),
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = point_size,
      opacity = opacity,
      colorbar = list(
        title = "Expression"
      )
    ),
    hoverinfo = "text",
    text = ~paste("Cluster:", cluster, "<br>Expression:", round(expr, 3))
  ) %>%
    layout(
      title = paste("Expression of", gene_symbol, "in UMAP"),
      scene = list(
        xaxis = list(title = "UMAP 1"),
        yaxis = list(title = "UMAP 2"),
        zaxis = list(title = "UMAP 3")
      )
    )
  
  return(p)
}

#' @title Save a Plotly Object to File
#' @description Saves a plotly visualization to an HTML file.
#' @param plot plotly object
#' @param file String file path
#' @param title String title for the plot
#' @return Invisible NULL, called for side effect
#' @export
save_plotly <- function(plot, file, title = "UMAP Plot") {
  # Add a title if not already present
  if (!"title" %in% names(plot$x$layout)) {
    plot <- plot %>% layout(title = title)
  }
  
  # Save widget
  htmlwidgets::saveWidget(
    widget = plot,
    file = file,
    selfcontained = TRUE
  )
}

#' @title Save a ggplot Object to File
#' @description Saves a ggplot visualization to a file.
#' @param plot ggplot object
#' @param file String file path
#' @param width Numeric width in inches
#' @param height Numeric height in inches
#' @param dpi Integer resolution
#' @param device String device type
#' @return Invisible NULL, called for side effect
#' @export
save_ggplot <- function(plot, file, width = 8, height = 8, dpi = 300, device = "png") {
  ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = device
  )
}

#' @title Generate UMAP Report
#' @description Generates a comprehensive report on UMAP visualization and clustering results.
#' @param seurat_obj Seurat object
#' @return List containing plots and statistics
#' @export
generate_umap_report <- function(seurat_obj) {
  # Validate input
  if (!("seurat_clusters" %in% colnames(seurat_obj@meta.data))) {
    return(list(
      error = "No clustering data available"
    ))
  }
  
  if (!("umap" %in% names(seurat_obj@reductions))) {
    return(list(
      error = "No UMAP data available"
    ))
  }
  
  # Generate plots and statistics
  cluster_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
  
  # Add sample information if available
  sample_plot <- NULL
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    sample_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "sample")
  }
  
  # Get cluster statistics
  cluster_stats <- table(seurat_obj$seurat_clusters)
  cluster_pct <- prop.table(cluster_stats) * 100
  
  # Create summary data frame
  summary_df <- data.frame(
    Cluster = names(cluster_stats),
    CellCount = as.numeric(cluster_stats),
    Percentage = as.numeric(cluster_pct)
  )
  
  # Sample distribution if available
  sample_dist <- NULL
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    sample_dist <- table(seurat_obj$seurat_clusters, seurat_obj$sample)
  }
  
  return(list(
    cluster_plot = cluster_plot,
    sample_plot = sample_plot,
    summary = summary_df,
    sample_distribution = sample_dist
  ))
}

#' @title Create Split UMAP
#' @description Creates a UMAP visualization split by a metadata factor for comparison.
#' @param seurat_obj Seurat object
#' @param split_by String metadata column to split by
#' @param reduction String reduction to use
#' @param ... Additional parameters to pass to DimPlot
#' @return ggplot object
#' @export
create_split_umap <- function(seurat_obj, split_by, reduction = "umap", ...) {
  if (!(split_by %in% colnames(seurat_obj@meta.data))) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("Column", split_by, "not found in metadata")) + 
             theme_void())
  }
  
  DimPlot(seurat_obj, 
          reduction = reduction, 
          group.by = "seurat_clusters", 
          split.by = split_by,
          ...)
}

#' @title Create Multi-Feature UMAP
#' @description Creates a multi-panel UMAP visualization showing expression of multiple genes.
#' @param seurat_obj Seurat object
#' @param features Vector of gene IDs to visualize
#' @param reduction String reduction to use
#' @param ncol Integer number of columns in the plot grid
#' @param ... Additional parameters to pass to FeaturePlot
#' @return ggplot object
#' @export
create_multi_feature_umap <- function(seurat_obj, features, reduction = "umap", ncol = 2, ...) {
  # Filter features to those available in the dataset
  valid_features <- features[features %in% rownames(seurat_obj)]
  
  if (length(valid_features) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "None of the requested features are in the dataset") + 
             theme_void())
  }
  
  FeaturePlot(seurat_obj, 
              features = valid_features, 
              reduction = reduction, 
              ncol = ncol,
              ...)
}