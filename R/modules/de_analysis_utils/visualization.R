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
#' @param group_by Optional string indicating which metadata column to group by (default is "seurat_clusters")
#' @return A pheatmap object visualizing the expression patterns
#' @export
createExpressionHeatmap <- function(seurat_obj, top_genes, cluster_labels = NULL, 
                                    active_clusters = NULL, group_by = NULL) {
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
  
  # Determine grouping to use
  if (is.null(group_by)) {
    # Try to infer from the input data or use default
    if ("comparison" %in% colnames(seurat_obj@meta.data)) {
      # If this is a result from comparison analysis
      comparison_text <- unique(seurat_obj@meta.data$comparison)[1]
      if (grepl("Sample", comparison_text, ignore.case = TRUE)) {
        group_by <- "sample"
      } else if (any(grepl("condition", colnames(seurat_obj@meta.data), ignore.case = TRUE))) {
        # Try to find a condition column
        condition_cols <- grep("condition", colnames(seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)
        if (length(condition_cols) > 0) {
          group_by <- condition_cols[1]
        } else {
          # Default to clusters
          group_by <- "seurat_clusters"
        }
      } else {
        # Default to clusters
        group_by <- "seurat_clusters"
      }
    } else {
      # Default to clusters
      group_by <- "seurat_clusters"
    }
  }
  
  # Create default group labels if not provided
  if (is.null(cluster_labels)) {
    # Get the unique groups
    if (group_by == "seurat_clusters") {
      unique_groups <- sort(unique(seurat_obj$seurat_clusters))
      default_prefix <- "Cluster"
    } else if (group_by == "sample") {
      unique_groups <- unique(seurat_obj$sample)
      default_prefix <- "Sample"
    } else if (group_by %in% colnames(seurat_obj@meta.data)) {
      unique_groups <- unique(seurat_obj@meta.data[[group_by]])
      # Create a nicer display name for the column
      default_prefix <- gsub("_", " ", group_by)
      default_prefix <- paste0(toupper(substr(default_prefix, 1, 1)), 
                               substr(default_prefix, 2, nchar(default_prefix)))
    } else {
      # Fallback
      unique_groups <- unique(as.character(seurat_obj$seurat_clusters))
      default_prefix <- "Group"
    }
    
    # Create the default labels
    cluster_labels <- setNames(
      paste(default_prefix, unique_groups),
      as.character(unique_groups)
    )
  }
  
  # Handle filtering based on group_by
  if (group_by == "seurat_clusters" && !is.null(active_clusters) && length(active_clusters) > 0) {
    # Filter by active clusters
    active_cells <- seurat_obj$seurat_clusters %in% active_clusters
    if (any(active_cells)) {
      seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
    }
  } else if (group_by == "sample" && !is.null(active_clusters) && length(active_clusters) > 0) {
    # Reinterpret active_clusters as active samples in this context
    active_cells <- seurat_obj$sample %in% active_clusters
    if (any(active_cells)) {
      seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
    }
  } else if (group_by %in% colnames(seurat_obj@meta.data) && 
             !is.null(active_clusters) && length(active_clusters) > 0) {
    # Reinterpret active_clusters as active values for this metadata column
    active_cells <- seurat_obj@meta.data[[group_by]] %in% active_clusters
    if (any(active_cells)) {
      seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
    }
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
  
  # Calculate means by group
  if (group_by == "seurat_clusters") {
    groups <- seurat_obj$seurat_clusters
    unique_groups <- sort(unique(groups))
  } else if (group_by == "sample") {
    groups <- seurat_obj$sample
    unique_groups <- unique(groups)
  } else if (group_by %in% colnames(seurat_obj@meta.data)) {
    groups <- seurat_obj@meta.data[[group_by]]
    unique_groups <- unique(groups)
  } else {
    # Fallback to clusters
    groups <- seurat_obj$seurat_clusters
    unique_groups <- sort(unique(groups))
  }
  
  # Check if we have groups
  if (length(unique_groups) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No groups to display") + 
             theme_void())
  }
  
  # Calculate means by group
  group_means <- sapply(unique_groups, function(group) {
    group_cells <- which(groups == group)
    if (length(group_cells) > 0) {
      rowMeans(expr_data[, group_cells, drop = FALSE])
    } else {
      rep(NA, nrow(expr_data))
    }
  })
  
  # Set column names using provided labels
  col_names <- sapply(unique_groups, function(x) {
    group_key <- as.character(x)
    if (group_key %in% names(cluster_labels)) {
      cluster_labels[[group_key]]
    } else {
      if (group_by == "seurat_clusters") {
        paste("Cluster", x)
      } else if (group_by == "sample") {
        paste("Sample", x)
      } else {
        paste("Group", x)
      }
    }
  })
  colnames(group_means) <- col_names
  
  # Scale data
  if (is.matrix(group_means) && nrow(group_means) > 1 && ncol(group_means) > 1) {
    # Only scale if we have enough data
    scaled_data <- t(scale(t(group_means)))
  } else {
    # Otherwise just use the original data
    scaled_data <- group_means
  }
  
  # Get gene labels
  gene_mapping <- seurat_obj@misc$gene_mapping
  gene_labels <- if (!is.null(gene_mapping) && all(rownames(scaled_data) %in% names(gene_mapping))) {
    gene_mapping[rownames(scaled_data)]
  } else {
    rownames(scaled_data)
  }
  gene_labels[is.na(gene_labels)] <- rownames(scaled_data)[is.na(gene_labels)]
  
  # Create heatmap title based on grouping
  if (group_by == "seurat_clusters") {
    title_prefix <- "Cluster"
  } else if (group_by == "sample") {
    title_prefix <- "Sample"
  } else {
    title_prefix <- gsub("_", " ", group_by)
    title_prefix <- paste0(toupper(substr(title_prefix, 1, 1)), 
                           substr(title_prefix, 2, nchar(title_prefix)))
  }
  
  # Create heatmap
  tryCatch({
    pheatmap(scaled_data,
             labels_row = gene_labels,
             labels_col = colnames(group_means),
             main = paste(title_prefix, "Expression Heatmap:", length(genes_present), "Genes"),
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
#' @description Creates a heatmap of cluster-specific marker genes.
#' @param seurat_obj Seurat object containing the data
#' @param genes Vector of genes to include in the heatmap
#' @param cluster_labels Named vector of cluster labels
#' @param active_clusters Vector of active cluster IDs
#' @param cluster_order Optional vector specifying the order of clusters
#' @return A pheatmap object showing cluster-specific gene expression
#' @keywords internal
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
  
  # Determine what grouping is being used
  current_ident <- DefaultAssay(seurat_obj)
  
  # Create default cluster labels if not provided
  if (is.null(cluster_labels)) {
    # Check what identity is being used
    unique_groups <- sort(unique(Idents(seurat_obj)))
    ident_name <- deparse(substitute(Idents(seurat_obj)))
    
    # Create default labels based on current identity
    if (ident_name == "sample") {
      group_name <- "Sample"
    } else if (ident_name == "seurat_clusters") {
      group_name <- "Cluster"
    } else {
      group_name <- "Group"
    }
    
    cluster_labels <- setNames(
      paste(group_name, unique_groups),
      as.character(unique_groups)
    )
  }
  
  # Subset to active clusters if provided
  if (!is.null(active_clusters) && length(active_clusters) > 0) {
    existing_active_clusters <- intersect(as.character(active_clusters), as.character(unique(Idents(seurat_obj))))
    
    if (length(existing_active_clusters) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "None of the active items exist in the current dataset") + 
               theme_void())
    }
    
    # Only use active groups that actually exist
    active_cells <- Idents(seurat_obj) %in% existing_active_clusters
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
  
  # Calculate means by identity group
  current_groups <- Idents(seurat_obj)
  unique_groups <- sort(unique(current_groups))
  
  # Check if we have groups
  if (length(unique_groups) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No groups to display") + 
             theme_void())
  }
  
  # Order groups if provided
  if (!is.null(cluster_order) && length(cluster_order) > 0) {
    # Convert to same type (character)
    ordered_groups <- as.character(cluster_order)
    # Only use groups that exist in our data
    ordered_groups <- intersect(ordered_groups, as.character(unique_groups))
    # If no matching groups (completely new set), use the default ordering
    if (length(ordered_groups) == 0) {
      ordered_groups <- as.character(unique_groups)
    } else {
      # Add any groups that weren't included in the order
      ordered_groups <- c(ordered_groups, setdiff(as.character(unique_groups), ordered_groups))
    }
    
    # Replace unique_groups with ordered version (convert to same type as unique_groups)
    unique_groups <- ordered_groups
    if (is.numeric(current_groups)) {
      unique_groups <- as.numeric(unique_groups)
    } else if (is.factor(current_groups)) {
      unique_groups <- factor(unique_groups, levels = levels(current_groups))
    }
  }
  
  # Calculate means for each group
  group_means <- sapply(unique_groups, function(group) {
    group_cells <- which(current_groups == group)
    if (length(group_cells) > 0) {
      rowMeans(expr_data[, group_cells, drop = FALSE])
    } else {
      rep(NA, nrow(expr_data))
    }
  })
  
  # Set column names using provided labels
  col_names <- sapply(unique_groups, function(x) {
    group_key <- as.character(x)
    if (group_key %in% names(cluster_labels)) {
      cluster_labels[[group_key]]
    } else {
      paste("Group", x)
    }
  })
  colnames(group_means) <- col_names
  
  # Scale data
  if (is.matrix(group_means) && nrow(group_means) > 1 && ncol(group_means) > 1) {
    # Only scale if we have enough data
    scaled_data <- t(scale(t(group_means)))
  } else {
    # Otherwise just use the original data
    scaled_data <- group_means
  }
  
  # Get gene labels
  gene_mapping <- seurat_obj@misc$gene_mapping
  gene_labels <- if (!is.null(gene_mapping) && all(rownames(scaled_data) %in% names(gene_mapping))) {
    gene_mapping[rownames(scaled_data)]
  } else {
    rownames(scaled_data)
  }
  gene_labels[is.na(gene_labels)] <- rownames(scaled_data)[is.na(gene_labels)]
  
  # Create heatmap with ordered rows (genes already ordered by group)
  tryCatch({
    # Get the identity name for the title
    if (identical(Idents(seurat_obj), seurat_obj$sample)) {
      ident_title <- "Sample"
    } else if (identical(Idents(seurat_obj), seurat_obj$seurat_clusters)) {
      ident_title <- "Cluster"
    } else {
      # Try to get a nice title from the column name
      ident_cols <- colnames(seurat_obj@meta.data)
      ident_match <- NULL
      for (col in ident_cols) {
        if (all(Idents(seurat_obj) == seurat_obj@meta.data[[col]])) {
          ident_match <- col
          break
        }
      }
      ident_title <- if (!is.null(ident_match)) {
        # Make a nice title from the column name
        gsub("_", " ", ident_match)
      } else {
        "Group"
      }
    }
    
    pheatmap(scaled_data,
             labels_row = gene_labels,
             labels_col = colnames(group_means),
             main = paste("Top", ident_title, "Specific Genes"),
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

#' @title Extract Volcano Plot Data
#' @description Extracts data used in a volcano plot for download
#' @param results Data frame containing differential expression results
#' @return Data frame formatted for volcano plot
#' @export
extractVolcanoData <- function(results) {
  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }
  
  # Create a data frame with the relevant columns
  volcano_data <- data.frame(
    gene_id = rownames(results),
    gene_symbol = results$gene,
    avg_log2FC = results$avg_log2FC,
    p_value = results$p_val,
    adjusted_p_value = results$p_val_adj,
    significant = (abs(results$avg_log2FC) > 0.25 & results$p_val_adj < 0.05),
    comparison = results$comparison
  )
  
  # Add other useful columns if they exist
  if ("pct.1" %in% colnames(results)) {
    volcano_data$pct_1 = results$pct.1
  }
  if ("pct.2" %in% colnames(results)) {
    volcano_data$pct_2 = results$pct.2
  }
  
  return(volcano_data)
}

#' @title Extract Heatmap Expression Data
#' @description Extracts expression data used in a heatmap for download
#' @param seurat_obj Seurat object containing the data
#' @param genes Vector of genes included in the heatmap
#' @param group_by Metadata column to group by (e.g., "seurat_clusters", "sample")
#' @param cluster_labels Optional named vector of cluster labels
#' @param active_groups Optional vector of active groups to include
#' @return List containing scaled and raw expression data
#' @export
extractHeatmapData <- function(seurat_obj, genes, group_by = "seurat_clusters", 
                               cluster_labels = NULL, active_groups = NULL) {
  if (is.null(genes) || length(genes) == 0) {
    return(NULL)
  }
  
  # Ensure all genes are in the Seurat object
  genes_present <- intersect(genes, rownames(seurat_obj))
  if (length(genes_present) == 0) {
    return(NULL)
  }
  
  # Get expression data for genes that exist in the dataset
  expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_present, , drop = FALSE]
  
  # Get grouping information
  if (group_by == "seurat_clusters") {
    groups <- seurat_obj$seurat_clusters
    unique_groups <- sort(unique(groups))
  } else if (group_by == "sample") {
    groups <- seurat_obj$sample
    unique_groups <- unique(groups)
  } else if (group_by %in% colnames(seurat_obj@meta.data)) {
    groups <- seurat_obj@meta.data[[group_by]]
    unique_groups <- unique(groups)
  } else {
    # Fallback to clusters
    group_by <- "seurat_clusters"
    groups <- seurat_obj$seurat_clusters
    unique_groups <- sort(unique(groups))
  }
  
  # Filter by active groups if provided
  if (!is.null(active_groups) && length(active_groups) > 0) {
    unique_groups <- intersect(as.character(unique_groups), as.character(active_groups))
  }
  
  # Calculate means by group
  group_means <- sapply(unique_groups, function(group) {
    group_cells <- which(groups == group)
    if (length(group_cells) > 0) {
      rowMeans(expr_data[, group_cells, drop = FALSE])
    } else {
      rep(NA, nrow(expr_data))
    }
  })
  
  # Set column names using provided labels
  col_names <- sapply(unique_groups, function(x) {
    group_key <- as.character(x)
    if (!is.null(cluster_labels) && group_key %in% names(cluster_labels)) {
      cluster_labels[[group_key]]
    } else {
      if (group_by == "seurat_clusters") {
        paste("Cluster", x)
      } else if (group_by == "sample") {
        paste("Sample", x)
      } else {
        paste("Group", x)
      }
    }
  })
  colnames(group_means) <- col_names
  
  # Calculate scaled data
  scaled_data <- NULL
  if (is.matrix(group_means) && nrow(group_means) > 1 && ncol(group_means) > 1) {
    scaled_data <- t(scale(t(group_means)))
  } else {
    scaled_data <- group_means
  }
  
  # Get gene symbols if available
  gene_mapping <- seurat_obj@misc$gene_mapping
  gene_symbols <- NULL
  if (!is.null(gene_mapping)) {
    gene_symbols <- gene_mapping[genes_present]
    gene_symbols[is.na(gene_symbols)] <- genes_present[is.na(gene_symbols)]
  } else {
    gene_symbols <- genes_present
  }
  
  # Create result object
  result <- list(
    raw_data = as.data.frame(group_means),
    scaled_data = as.data.frame(scaled_data),
    genes = genes_present,
    gene_symbols = gene_symbols,
    group_by = group_by,
    groups = unique_groups
  )
  
  return(result)
}

#' @title Extract Gene Expression Boxplot Data
#' @description Extracts gene expression data used in a boxplot for download
#' @param seurat_obj Seurat object containing the data
#' @param gene_id String ID of the gene being visualized
#' @param group_by String metadata column to group by
#' @param active_groups Optional vector of active groups to include
#' @return A data frame with expression values and grouping information
#' @export
extractBoxplotData <- function(seurat_obj, gene_id, group_by = "seurat_clusters", active_groups = NULL) {
  # Check if gene exists
  if (!(gene_id %in% rownames(seurat_obj))) {
    return(NULL)
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
  
  # Get grouping information
  if (!(group_by %in% colnames(seurat_obj@meta.data))) {
    return(NULL)
  }
  
  # Extract group values for each cell
  group_values <- seurat_obj@meta.data[[group_by]]
  
  # Create data frame with expression and group information
  plot_data <- data.frame(
    cell_id = names(expr_values),
    expression = as.numeric(expr_values),
    group = group_values,
    stringsAsFactors = FALSE
  )
  
  # Filter by active groups if provided
  if (!is.null(active_groups) && length(active_groups) > 0) {
    active_groups <- as.character(active_groups)
    plot_data <- plot_data[plot_data$group %in% active_groups, ]
  }
  
  # Add sample info if available
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    plot_data$sample <- seurat_obj@meta.data[plot_data$cell_id, "sample"]
  }
  
  # Add cluster info if available and not already included
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data) && group_by != "seurat_clusters") {
    plot_data$cluster <- seurat_obj@meta.data[plot_data$cell_id, "seurat_clusters"]
  }
  
  # Add metadata
  attr(plot_data, "gene_id") <- gene_id
  attr(plot_data, "gene_symbol") <- gene_symbol
  attr(plot_data, "group_by") <- group_by
  
  return(plot_data)
}