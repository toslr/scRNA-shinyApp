#' @title Create Volcano Plot
#' @description Creates a volcano plot visualization for differential expression results,
#'   highlighting significantly differentially expressed genes.
#' @param results Data frame containing differential expression results with p-values and fold changes
#' @return A ggplot object with the volcano plot visualization
#' @export
createVolcanoPlot <- function(results) {
  ggplot(results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = abs(avg_log2FC) > 0.25 & p_val_adj < 0.05)) +
    scale_color_manual(values = c("grey", "red")) +
    theme_classic() +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = unique(results$comparison),
         color = "Significant") +
    theme(legend.position = "bottom")
}

#' @title Create Expression Heatmap
#' @description Creates a heatmap of gene expression across clusters for a set of genes.
#'   The heatmap shows scaled expression values to highlight differential patterns.
#' @param seurat_obj Seurat object containing the data
#' @param top_genes Vector of gene IDs to include in the heatmap
#' @param cluster_labels Optional named vector of cluster labels
#' @param active_clusters Optional vector of active cluster IDs to include
#' @return A pheatmap object visualizing the expression patterns
#' @export
createExpressionHeatmap <- function(seurat_obj, top_genes, cluster_labels = NULL, active_clusters = NULL) {
  # Safety check for empty gene list
  if (length(top_genes) == 0) {
    # Create an empty plot with a message
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No genes to display") + 
             theme_void())
  }
  
  # Ensure cluster_labels is not a function (reactive value)
  if (is.function(cluster_labels)) {
    tryCatch({
      cluster_labels <- cluster_labels()
    }, error = function(e) {
      cluster_labels <- NULL
    })
  }
  
  # Create default cluster labels if not provided
  if (is.null(cluster_labels)) {
    unique_clusters <- sort(unique(seurat_obj$seurat_clusters))
    cluster_labels <- setNames(
      paste("Cluster", unique_clusters),
      as.character(unique_clusters)
    )
  }
  
  # Subset to active clusters if provided
  if (!is.null(active_clusters) && length(active_clusters) > 0) {
    active_cells <- seurat_obj$seurat_clusters %in% active_clusters
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
  }
  
  # Ensure all top genes are in the Seurat object
  genes_present <- intersect(top_genes, rownames(seurat_obj))
  if (length(genes_present) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "None of the selected genes are present in the dataset") + 
             theme_void())
  }
  
  # Get expression data for genes that exist in the dataset
  expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_present, , drop = FALSE]
  
  # Calculate cluster means
  clusters <- seurat_obj$seurat_clusters
  unique_clusters <- sort(unique(clusters))
  
  # Check if we have clusters
  if (length(unique_clusters) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No clusters to display") + 
             theme_void())
  }
  
  cluster_means <- sapply(unique_clusters, function(clust) {
    cluster_cells <- which(clusters == clust)
    if (length(cluster_cells) > 0) {
      rowMeans(expr_data[, cluster_cells, drop = FALSE])
    } else {
      rep(NA, nrow(expr_data))
    }
  })
  
  # Set column names using cluster labels
  col_names <- sapply(unique_clusters, function(x) {
    cluster_key <- as.character(x)
    if (cluster_key %in% names(cluster_labels)) {
      cluster_labels[[cluster_key]]
    } else {
      paste("Cluster", x)
    }
  })
  colnames(cluster_means) <- col_names
  
  # Scale data
  if (is.matrix(cluster_means) && nrow(cluster_means) > 1 && ncol(cluster_means) > 1) {
    # Only scale if we have enough data
    scaled_data <- t(scale(t(cluster_means)))
  } else {
    # Otherwise just use the original data
    scaled_data <- cluster_means
  }
  
  # Get gene labels
  gene_mapping <- seurat_obj@misc$gene_mapping
  gene_labels <- if (!is.null(gene_mapping) && all(rownames(scaled_data) %in% names(gene_mapping))) {
    gene_mapping[rownames(scaled_data)]
  } else {
    rownames(scaled_data)
  }
  gene_labels[is.na(gene_labels)] <- rownames(scaled_data)[is.na(gene_labels)]
  
  # Create heatmap
  tryCatch({
    pheatmap(scaled_data,
             labels_row = gene_labels,
             labels_col = colnames(cluster_means),
             main = paste("Expression Heatmap:", length(genes_present), "Genes"),
             angle_col = 45,
             fontsize_row = 10,
             cluster_cols = FALSE,
             treeheight_row = 0,
             draw = TRUE)
  }, error = function(e) {
    # If pheatmap fails, return a ggplot error message
    ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Error creating heatmap:", e$message)) + 
      theme_void()
  })
}

#' @title Create General Heatmap
#' @description Creates a cluster-specific gene expression heatmap that highlights marker genes
#'   for each cluster. Unlike the standard expression heatmap, this preserves gene order by cluster
#'   to help visualize cluster-specific expression patterns.
#' @param seurat_obj Seurat object containing the data
#' @param genes Vector of gene IDs to include in the heatmap, typically ordered by cluster
#' @param cluster_labels Optional named vector of cluster labels
#' @param active_clusters Optional vector of active cluster IDs to include
#' @param cluster_order Optional vector specifying the ordering of clusters in the visualization
#' @return A pheatmap object visualizing the cluster-specific gene expression
#' @export
createGeneralHeatmap <- function(seurat_obj, genes, cluster_labels = NULL, active_clusters = NULL, cluster_order = NULL) {
  # Safety check for empty gene list
  if (length(genes) == 0) {
    # Create an empty plot with a message
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No genes to display") + 
             theme_void())
  }
  
  # Ensure cluster_labels is not a function (reactive value)
  if (is.function(cluster_labels)) {
    tryCatch({
      cluster_labels <- cluster_labels()
    }, error = function(e) {
      cluster_labels <- NULL
    })
  }
  
  # Create default cluster labels if not provided
  if (is.null(cluster_labels)) {
    unique_clusters <- sort(unique(seurat_obj$seurat_clusters))
    cluster_labels <- setNames(
      paste("Cluster", unique_clusters),
      as.character(unique_clusters)
    )
  }
  
  # Subset to active clusters if provided
  if (!is.null(active_clusters) && length(active_clusters) > 0) {
    existing_active_clusters <- intersect(active_clusters, unique(seurat_obj$seurat_clusters))
    
    if (length(existing_active_clusters) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "None of the active clusters exist in the current dataset") + 
               theme_void())
    }
    
    # Only use active clusters that actually exist
    active_cells <- seurat_obj$seurat_clusters %in% existing_active_clusters
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
  }
  
  # Ensure all genes are in the Seurat object
  genes_present <- intersect(genes, rownames(seurat_obj))
  if (length(genes_present) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "None of the selected genes are present in the dataset") + 
             theme_void())
  }
  
  # Get expression data for selected genes
  expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_present, , drop = FALSE]
  
  # Calculate cluster means
  clusters <- seurat_obj$seurat_clusters
  unique_clusters <- sort(unique(clusters))
  
  # Check if we have clusters
  if (length(unique_clusters) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No clusters to display") + 
             theme_void())
  }
  
  # Order clusters if provided
  if (!is.null(cluster_order) && length(cluster_order) > 0) {
    # Convert to numeric
    ordered_clusters <- as.numeric(cluster_order)
    # Only use clusters that exist in our data
    ordered_clusters <- intersect(ordered_clusters, unique_clusters)
    # If no matching clusters (completely new set), use the default ordering
    if (length(ordered_clusters) == 0) {
      ordered_clusters <- unique_clusters
    } else {
      # Add any clusters that weren't included in the order
      ordered_clusters <- c(ordered_clusters, setdiff(unique_clusters, ordered_clusters))
    }
    
    # Replace unique_clusters with ordered version
    unique_clusters <- ordered_clusters
  }
  
  cluster_means <- sapply(unique_clusters, function(clust) {
    cluster_cells <- which(clusters == clust)
    if (length(cluster_cells) > 0) {
      rowMeans(expr_data[, cluster_cells, drop = FALSE])
    } else {
      rep(NA, nrow(expr_data))
    }
  })
  
  # Set column names using cluster labels
  col_names <- sapply(unique_clusters, function(x) {
    cluster_key <- as.character(x)
    if (cluster_key %in% names(cluster_labels)) {
      cluster_labels[[cluster_key]]
    } else {
      paste("Cluster", x)
    }
  })
  colnames(cluster_means) <- col_names
  
  # Scale data
  if (is.matrix(cluster_means) && nrow(cluster_means) > 1 && ncol(cluster_means) > 1) {
    # Only scale if we have enough data
    scaled_data <- t(scale(t(cluster_means)))
  } else {
    # Otherwise just use the original data
    scaled_data <- cluster_means
  }
  
  # Get gene labels
  gene_mapping <- seurat_obj@misc$gene_mapping
  gene_labels <- if (!is.null(gene_mapping) && all(rownames(scaled_data) %in% names(gene_mapping))) {
    gene_mapping[rownames(scaled_data)]
  } else {
    rownames(scaled_data)
  }
  gene_labels[is.na(gene_labels)] <- rownames(scaled_data)[is.na(gene_labels)]
  
  # Create heatmap with ordered rows (genes already ordered by cluster)
  tryCatch({
    pheatmap(scaled_data,
             labels_row = gene_labels,
             labels_col = colnames(cluster_means),
             main = "Top Cluster-Specific Genes",
             angle_col = 45,
             fontsize_row = 8,
             cluster_cols = FALSE,
             cluster_rows = FALSE,  # Don't cluster rows to maintain diagonal pattern
             treeheight_row = 0,
             draw = TRUE)
  }, error = function(e) {
    # If pheatmap fails, return a ggplot error message
    ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Error creating heatmap:", e$message)) + 
      theme_void()
  })
}

#' @title Create Gene Expression Boxplot with Statistics
#' @description Creates a boxplot of gene expression across clusters or conditions with statistical
#'   significance annotations between groups.
#' @param seurat_obj Seurat object containing the data
#' @param gene_id String ID of the gene to visualize
#' @param group_by String variable to group by (default: "seurat_clusters")
#' @param comparisons List of pairs for statistical comparison (default: NULL)
#' @param active_groups Optional vector of active groups to include (default: NULL)
#' @param cluster_labels Optional named vector of labels for groups (default: NULL)
#' @return A ggplot object with the boxplot visualization
#' @export
createExpressionBoxplot <- function(seurat_obj, gene_id, group_by = "seurat_clusters",
                                    comparisons = NULL, active_groups = NULL, 
                                    cluster_labels = NULL) {
  # Safety check for gene existence in the dataset
  if (is.null(gene_id) || gene_id == "" || !(gene_id %in% rownames(seurat_obj))) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("Gene not found in dataset. Please try another gene.")) + 
             theme_void())
  }
  
  # Get gene symbol if available
  gene_symbol <- if (!is.null(seurat_obj@misc$gene_mapping) && 
                     gene_id %in% names(seurat_obj@misc$gene_mapping) &&
                     !is.na(seurat_obj@misc$gene_mapping[gene_id])) {
    seurat_obj@misc$gene_mapping[gene_id]
  } else {
    gene_id
  }
  
  # Handle group variable
  if (!(group_by %in% colnames(seurat_obj@meta.data))) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("Group variable", group_by, "not found in metadata")) + 
             theme_void())
  }
  
  # Ensure cluster_labels is not a function (reactive value)
  if (is.function(cluster_labels)) {
    tryCatch({
      cluster_labels <- cluster_labels()
    }, error = function(e) {
      cluster_labels <- NULL
    })
  }
  
  # Get available groups properly
  available_groups <- NULL
  tryCatch({
    if (is.data.frame(seurat_obj[[group_by]])) {
      available_groups <- unlist(seurat_obj[[group_by]])
    } else {
      available_groups <- seurat_obj[[group_by]]
    }
    available_groups <- as.character(unique(available_groups))
  }, error = function(e) {
    available_groups <- as.character(unique(seurat_obj@meta.data[[group_by]]))
  })
  
  if (is.null(available_groups) || length(available_groups) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "Could not determine available groups") + 
             theme_void())
  }
  
  # Process active_groups
  if (is.null(active_groups) || length(active_groups) == 0) {
    active_groups <- available_groups
  } else {
    active_groups <- as.character(active_groups)
    valid_active_groups <- intersect(active_groups, available_groups)
    
    if (length(valid_active_groups) == 0 && group_by == "seurat_clusters") {
      # Try zero-based to one-based adjustment for cluster numbers
      active_plus_one <- as.character(as.numeric(active_groups) + 1)
      valid_plus_one <- intersect(active_plus_one, available_groups)
      
      active_minus_one <- as.character(as.numeric(active_groups) - 1)
      valid_minus_one <- intersect(active_minus_one, available_groups)
      
      if (length(valid_plus_one) > length(valid_minus_one)) {
        valid_active_groups <- valid_plus_one
      } else if (length(valid_minus_one) > 0) {
        valid_active_groups <- valid_minus_one
      }
    }
    
    if (length(valid_active_groups) == 0) {
      valid_active_groups <- available_groups
    }
    
    active_groups <- valid_active_groups
  }
  
  # Filter for cells matching our active groups
  cells_to_keep <- seurat_obj@meta.data[[group_by]] %in% active_groups
  
  # Check if any cells match
  if (sum(cells_to_keep) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No cells match the active groups selection") + 
             theme_void())
  }
  
  # Subset to matching cells
  filtered_seurat <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
  
  # Get expression data
  expr_data <- GetAssayData(filtered_seurat, slot = "data")[gene_id, ]
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Expression = expr_data,
    Group = filtered_seurat@meta.data[[group_by]],
    stringsAsFactors = FALSE
  )
  
  # Handle factor conversion and labels for plotting
  if (group_by == "seurat_clusters" && !is.null(cluster_labels)) {
    plot_data$Group <- as.character(plot_data$Group)
    plot_data$GroupLabel <- sapply(plot_data$Group, function(x) {
      if (x %in% names(cluster_labels)) {
        return(cluster_labels[[x]])
      } else {
        return(paste("Cluster", x))
      }
    })
    
    # Create ordered factor based on numeric cluster IDs
    # First extract the unique groups and their order
    unique_groups <- unique(plot_data$Group)
    group_order <- order(as.numeric(unique_groups))
    ordered_levels <- unique_groups[group_order]
    
    # Map the ordered clusters to their labels
    ordered_labels <- sapply(ordered_levels, function(x) {
      if (x %in% names(cluster_labels)) {
        return(cluster_labels[[x]])
      } else {
        return(paste("Cluster", x))
      }
    })
    
    # Set factor levels for ordered plotting
    plot_data$GroupLabel <- factor(plot_data$GroupLabel, levels = ordered_labels)
    group_var_to_plot <- "GroupLabel"
  } else {
    # For non-cluster variables, just order numerically if possible
    if (all(grepl("^[0-9]+$", plot_data$Group))) {
      plot_data$Group <- factor(plot_data$Group, 
                                levels = unique(plot_data$Group)[order(as.numeric(unique(plot_data$Group)))])
    } else {
      plot_data$Group <- factor(plot_data$Group)
    }
    group_var_to_plot <- "Group"
  }
  
  # Create base boxplot
  p <- ggplot(plot_data, aes_string(x = group_var_to_plot, y = "Expression")) +
    geom_boxplot(outlier.size = 0.5, fill = "lightblue") +
    geom_jitter(width = 0.2, height = 0, size = 0.5, alpha = 0.3) +
    labs(title = paste("Expression of", gene_symbol),
         y = "Normalized expression",
         x = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12))
  
  # Add statistical comparisons if requested
  if (!is.null(comparisons) && requireNamespace("ggpubr", quietly = TRUE)) {
    # Get the levels from the factor for consistent ordering
    actual_groups <- levels(plot_data[[group_var_to_plot]])
    
    if (length(actual_groups) > 1) {
      if (length(actual_groups) <= 3) {
        # For few groups, compare all pairs
        valid_comparisons <- combn(actual_groups, 2, simplify = FALSE)
      } else {
        # For many groups, just compare adjacent pairs
        valid_comparisons <- lapply(1:(length(actual_groups)-1), function(i) {
          c(actual_groups[i], actual_groups[i+1])
        })
      }
      
      # Add statistical comparisons
      p <- p + ggpubr::stat_compare_means(
        comparisons = valid_comparisons,
        method = "wilcox.test",
        label = "p.format"
      )
    }
  }
  
  return(p)
}