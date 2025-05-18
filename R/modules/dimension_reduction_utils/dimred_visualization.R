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

#' @title Create 2D UMAP Plot with Improved Filtering
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
  using_condition_labels <- FALSE
  
  # First handle filtering based on active_items if provided
  if (!is.null(active_items) && length(active_items) > 0) {
    using_condition_labels <- FALSE
    if (color_by != "cluster" && color_by != "sample" && color_by != "gene" &&
        "condition_label" %in% colnames(seurat_obj@meta.data) && 
        color_by %in% colnames(seurat_obj@meta.data)) {
      using_condition_labels <- TRUE
    }
    if (color_by == "cluster" && "seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      # Extract cluster numbers from formatted active items (e.g., "0:Cluster 0" -> "0")
      numeric_clusters <- sapply(strsplit(active_items, ":"), function(x) as.numeric(x[1]))
      cells_to_keep <- seurat_obj$seurat_clusters %in% numeric_clusters
      
      # Only filter if we have matches
      if (any(cells_to_keep)) {
        seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
      } else {
        # Return message if no cells match
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "No cells match the selected clusters") + 
                 theme_void())
      }
    } else if (color_by == "sample" && "sample" %in% colnames(seurat_obj@meta.data)) {
      # For sample filtering
      # Check if we need to extract sample names from formatted items
      if (any(grepl(":", active_items))) {
        # Items are formatted like "SampleName:Label"
        sample_ids <- sapply(strsplit(active_items, ":"), function(x) x[1])
      } else {
        sample_ids <- active_items
      }
      
      cells_to_keep <- seurat_obj$sample %in% sample_ids
      
      if (any(cells_to_keep)) {
        seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
      } else {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "No cells match the selected samples") + 
                 theme_void())
      }
    } else if (color_by %in% colnames(seurat_obj@meta.data)) {
      # Extract condition IDs from active_items if they're in "id:name" format
      if (any(grepl(":", active_items))) {
        # Split by the first colon only, in case there are multiple colons
        condition_ids <- sapply(strsplit(active_items, ":", fixed = TRUE), function(x) x[1])
      } else {
        condition_ids <- active_items
      }
      
      # Use original values (not labels) for filtering
      cells_to_keep <- seurat_obj@meta.data[[color_by]] %in% condition_ids
    }
  }
  
  # Then handle different coloring options
  if (color_by == "cluster" && "seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    # Check if we have a cluster_label column to use instead
    if ("cluster_label" %in% colnames(seurat_obj@meta.data)) {
      # Store original cluster labels (they might be in 0:name format already)
      original_labels <- seurat_obj$cluster_label
      
      # Create a column with just numbers for labeling the plot
      seurat_obj$cluster_number <- seurat_obj$seurat_clusters
      
      # Create a column with shortened legend labels (0:abc format)
      seurat_obj$shortened_labels <- sapply(original_labels, function(label) {
        parts <- strsplit(label, ":")[[1]]
        if (length(parts) == 2) {
          # Label is already in number:name format
          cluster_num <- parts[1]
          full_name <- parts[2]
          
          # Get first 3 letters of the name
          short_name <- if (nchar(full_name) >= 3) {
            substr(full_name, 1, 3)
          } else {
            full_name
          }
          
          paste0(cluster_num, ":", short_name)
        } else {
          # Label is not in number:name format yet
          cluster_num <- as.character(seurat_obj$seurat_clusters[seurat_obj$cluster_label == label][1])
          
          # Get first 3 letters of the name
          short_name <- if (nchar(label) >= 3) {
            substr(label, 1, 3)
          } else {
            label
          }
          
          paste0(cluster_num, ":", short_name)
        }
      })
      
      # Create ordered factor for plotting
      cluster_ids <- as.character(seurat_obj$seurat_clusters)
      unique_clusters <- sort(as.numeric(unique(cluster_ids)))
      
      # Create an ordering map to maintain numeric order
      ordering_map <- data.frame(
        cluster_id = as.character(unique_clusters),
        label = NA,
        stringsAsFactors = FALSE
      )
      
      # Map the shortened labels to the ordering
      for (i in 1:nrow(ordering_map)) {
        cluster_id <- ordering_map$cluster_id[i]
        # Find the first matching shortened label
        matching_idx <- which(cluster_ids == cluster_id)[1]
        if (!is.na(matching_idx)) {
          ordering_map$label[i] <- seurat_obj$shortened_labels[matching_idx]
        }
      }
      
      # Set the factor levels in numerical order
      seurat_obj$shortened_labels <- factor(seurat_obj$shortened_labels, 
                                            levels = ordering_map$label[!is.na(ordering_map$label)])
      
      # Create plot with customized labels
      if (label) {
        # If labeling is requested, use DimPlot with numbers as labels
        p <- DimPlot(seurat_obj, 
                     reduction = reduction, 
                     group.by = "shortened_labels",  # Group by the shortened labels for legend
                     label = TRUE,  # Show labels
                     label.size = 4,  # Adjust size as needed
                     repel = TRUE,  # Avoid overlapping
                     pt.size = pt_size,
                     ...)
        
        # Replace labels with just the cluster numbers
        p$labels <- as.character(seurat_obj$seurat_clusters[match(p$labels, seurat_obj$shortened_labels)])
        
        # Customize legend title
        p <- p + labs(color = "Cluster")
      } else {
        # If no labeling, just use the shortened labels for coloring
        p <- DimPlot(seurat_obj, 
                     reduction = reduction, 
                     group.by = "shortened_labels",
                     label = FALSE,
                     pt.size = pt_size,
                     ...)
        
        # Customize legend title
        p <- p + labs(color = "Cluster")
      }
      
      return(p)
    } 
  } else if (color_by == "sample" && "sample" %in% colnames(seurat_obj@meta.data)) {
    # Check if we have a sample_label column to use instead
    if ("sample_label" %in% colnames(seurat_obj@meta.data)) {
      # Use the sample_label column for visualization  
      return(DimPlot(seurat_obj, 
                     reduction = reduction,
                     group.by = "sample_label", 
                     pt.size = pt_size,
                     ...))
    }
    
    # Otherwise, use standard sample coloring (already filtered above if needed)
    return(DimPlot(seurat_obj, 
                   reduction = reduction, 
                   group.by = "sample",
                   pt.size = pt_size,
                   ...))
    
  } else if (color_by == "gene" && !is.null(gene_id) && gene_id %in% rownames(seurat_obj)) {
    # For gene expression, we use the object (filtered above if needed)
    return(FeaturePlot(seurat_obj, 
                       features = gene_id, 
                       reduction = reduction,
                       pt.size = pt_size,
                       ...))
    
  } else if (color_by %in% colnames(seurat_obj@meta.data)) {
    # Color by any metadata column (already filtered above if needed)
    if (using_condition_labels) {
      return(DimPlot(seurat_obj, 
                     reduction = reduction, 
                     group.by = "condition_label",
                     pt.size = pt_size,
                     ...))
    } else {
      return(DimPlot(seurat_obj, 
                     reduction = reduction, 
                     group.by = color_by,
                     pt.size = pt_size,
                     ...))
    }
  } else {
    # Default empty plot with error message
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "Please select a valid coloring option") + 
             theme_void())
  }
}


create_3d_umap_plot <- function(seurat_obj, color_by = "cluster", reduction = "umap3d",
                                point_size = 3, opacity = 0.7, active_items = NULL) {
  
  # Handle subsetting based on active items
  if (!is.null(active_items) && length(active_items) > 0) {
    if (color_by == "cluster" && "seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      # For clusters, we can handle numeric values directly
      numeric_clusters <- as.numeric(active_items)
      cells_to_keep <- seurat_obj$seurat_clusters %in% numeric_clusters
    } else if (color_by == "sample" && "sample" %in% colnames(seurat_obj@meta.data)) {
      # For samples, use direct IDs stored in active_items
      cells_to_keep <- seurat_obj$sample %in% active_items
    } else if (color_by %in% colnames(seurat_obj@meta.data)) {
      # For condition columns, straight matching by original values
      cells_to_keep <- seurat_obj@meta.data[[color_by]] %in% active_items
    } else {
      cells_to_keep <- rep(TRUE, ncol(seurat_obj))
    }
    
    # Check if any cells match the filter
    if (!any(cells_to_keep)) {
      print("WARNING: No cells match the filter criteria")
      return(plotly_empty() %>% 
               layout(title = "No cells match the current filter settings"))
    }
    
    # Subset to matching cells
    plot_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
  } else {
    # No active items filtering - use all cells
    plot_obj <- seurat_obj
    print("Using all cells (no active item filtering)")
  }
    
  # Extract UMAP coordinates
  umap_data <- Embeddings(plot_obj[[reduction]])
  
  # Handle different coloring options
  if (color_by == "cluster" && "seurat_clusters" %in% colnames(plot_obj@meta.data)) {
    # Check if cluster_label column exists
    if ("cluster_label" %in% colnames(plot_obj@meta.data)) {
      # Get the cluster labels 
      cluster_ids <- as.character(plot_obj$seurat_clusters)
      cluster_labels <- plot_obj$cluster_label
      
      # Handle both formats: with or without colon
      short_labels <- sapply(cluster_labels, function(label) {
        if (grepl(":", label)) {
          # If label contains colon, extract cluster number and first 3 letters of name
          parts <- strsplit(label, ":")[[1]]
          cluster_num <- parts[1]
          name_part <- if (length(parts) > 1 && nchar(parts[2]) >= 3) {
            substr(parts[2], 1, 3)
          } else if (length(parts) > 1) {
            parts[2]
          } else {
            "Cluster"
          }
          paste0(cluster_num, ":", name_part)
        } else {
          # Otherwise, create a label with cluster number and first 3 letters
          cluster_id <- cluster_ids[which(cluster_labels == label)[1]]
          name_part <- if (nchar(label) >= 3) {
            substr(label, 1, 3)
          } else {
            label
          }
          paste0(cluster_id, ":", name_part)
        }
      })
      
      # Create ordered factor for short labels
      unique_short_labels <- unique(short_labels)
      # Extract numeric cluster IDs from the labels
      cluster_numbers <- sapply(strsplit(unique_short_labels, ":"), function(x) {
        as.numeric(x[1])
      })
      # Order by numeric cluster ID
      ordered_indices <- order(cluster_numbers)
      ordered_short_labels <- unique_short_labels[ordered_indices]
      
      color_var <- factor(short_labels, levels = ordered_short_labels)
      plot_title <- "3D UMAP - Colored by Clusters"
      hover_text <- paste("Cluster:", cluster_ids, "<br>Label:", cluster_labels)
    } else {
      # Regular cluster coloring without custom labels
      color_var <- as.factor(plot_obj$seurat_clusters)
      plot_title <- "3D UMAP - Colored by Clusters"
      hover_text <- paste("Cluster:", color_var)
    }
  } else if (color_by == "sample" && "sample" %in% colnames(plot_obj@meta.data)) {
    # Check if we have a sample_label column to use instead
    if ("sample_label" %in% colnames(plot_obj@meta.data)) {
      color_var <- as.factor(plot_obj$sample_label)
      plot_title <- "3D UMAP - Colored by Sample"
      hover_text <- paste("Sample:", plot_obj$sample, "<br>Label:", color_var)
    } else {
      color_var <- as.factor(plot_obj$sample)
      plot_title <- "3D UMAP - Colored by Sample"
      hover_text <- paste("Sample:", color_var)
    }
  } else if (color_by == "condition_label" && "condition_label" %in% colnames(plot_obj@meta.data)) {
    # Handle condition label coloring
    color_var <- as.factor(plot_obj$condition_label)
    condition_column <- names(plot_obj@meta.data)[grep("condition", names(plot_obj@meta.data), ignore.case = TRUE)[1]]
    original_value <- plot_obj@meta.data[[condition_column]]
    plot_title <- paste0("3D UMAP - Colored by ", condition_column)
    hover_text <- paste(condition_column, ":", original_value, "<br>Label:", color_var)
  } else if (color_by %in% colnames(plot_obj@meta.data)) {
    # Check if we have condition labels for this column
    condition_column <- color_by
    if ("condition_label" %in% colnames(plot_obj@meta.data) && 
        plot_obj@meta.data[[condition_column]][1] == plot_obj@meta.data[[color_by]][1]) {
      # Use the condition_label column if it matches our current color_by
      color_var <- as.factor(plot_obj$condition_label)
      plot_title <- paste0("3D UMAP - Colored by ", color_by)
      hover_text <- paste(color_by, ":", plot_obj@meta.data[[color_by]], "<br>Label:", color_var)
    } else {
      color_var <- as.factor(plot_obj@meta.data[[color_by]])
      plot_title <- paste0("3D UMAP - Colored by ", color_by)
      hover_text <- paste(color_by, ":", color_var)
    }
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

#' @title Extract PCA Data for Download
#' @description Extracts PCA eigenvalues, variance explained, and coordinates for download
#' @param seurat_obj Seurat object containing PCA results
#' @return List containing multiple data frames with PCA information
#' @export
extractPCAData <- function(seurat_obj) {
  # Check if PCA exists
  if (!("pca" %in% names(seurat_obj@reductions))) {
    return(NULL)
  }
  
  # Extract eigenvalues and variance explained
  pca_obj <- seurat_obj[["pca"]]
  eigenvalues <- (Stdev(pca_obj))^2
  variance_explained <- eigenvalues / sum(eigenvalues)
  cumulative_variance <- cumsum(variance_explained)
  
  # Create summary data frame
  pca_summary <- data.frame(
    PC = 1:length(eigenvalues),
    Eigenvalue = eigenvalues,
    VarianceExplained = variance_explained,
    CumulativeVariance = cumulative_variance
  )
  
  # Extract cell embeddings (PCA coordinates)
  pca_embeddings <- as.data.frame(Embeddings(pca_obj))
  
  # Add cell IDs
  pca_embeddings$cell_id <- rownames(pca_embeddings)
  
  # Reorder columns to put cell_id first
  pca_embeddings <- pca_embeddings[, c("cell_id", setdiff(colnames(pca_embeddings), "cell_id"))]
  
  # Get sample information if available
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    sample_info <- seurat_obj@meta.data[, "sample", drop = FALSE]
    
    # Add to embeddings
    pca_embeddings$sample <- sample_info[rownames(pca_embeddings), "sample"]
  }
  
  # Get cluster information if available
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cluster_info <- seurat_obj@meta.data[, "seurat_clusters", drop = FALSE]
    
    # Add to embeddings
    pca_embeddings$cluster <- cluster_info[rownames(pca_embeddings), "seurat_clusters"]
  }
  
  # Return as a list
  return(list(
    summary = pca_summary,
    embeddings = pca_embeddings
  ))
}