# R/modules/dimension_reduction_utils/dimred_utils.R

# ==================== Helper Functions for Dimension Reduction ====================

# Find elbow point for optimal dimensions
find_elbow <- function(x, y) {
  # Focus on the rate of change and look for where the rate of change stabilizes
  diffs <- diff(y) / diff(x)
  
  window_size <- 3 # Adjustable
  rolling_std <- sapply(1:(length(diffs) - window_size), function(i) {
    sd(diffs[i:(i + window_size)])
  })
  
  threshold <- mean(rolling_std) * 0.1
  elbow_idx <- which(rolling_std < threshold)[1]
  
  return(x[elbow_idx])
}

# Save elbow plot to file
saveElbowPlot <- function(file, seurat_obj, suggested_dims) {
  png(file, width = 3000, height = 1800, res = 300)
  
  print(ElbowPlot(seurat_obj, 
                  ndims = ncol(Embeddings(seurat_obj, "pca"))) +
          geom_vline(xintercept = suggested_dims, 
                     color = "red", 
                     linetype = "dashed"))
  
  dev.off()
}

# Extend color palette for potentially large number of clusters
extendColorPalette <- function(n_colors) {
  if (n_colors <= 8) {
    return(RColorBrewer::brewer.pal(8, "Set2")[1:n_colors])
  } else {
    return(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_colors))
  }
}

# Function to create gene heatmap for top marker genes per cluster
createClusterMarkerHeatmap <- function(seurat_obj, n_genes_per_cluster = 5) {
  # Check if clustering has been performed
  if (!("seurat_clusters" %in% colnames(seurat_obj@meta.data))) {
    return(NULL)
  }
  
  # List to store top markers
  top_markers <- list()
  all_markers <- c()
  
  # For each cluster, find markers
  clusters <- sort(unique(seurat_obj$seurat_clusters))
  
  withProgress(message = 'Finding cluster markers', value = 0, {
    for (i in seq_along(clusters)) {
      incProgress(1/length(clusters),
                  detail = paste("Processing cluster", clusters[i]))
      
      # Find markers for this cluster vs all others
      markers <- FindMarkers(
        seurat_obj,
        ident.1 = clusters[i],
        min.pct = 0.25,
        logfc.threshold = 0.25
      )
      
      # Sort by p-value
      markers <- markers[order(markers$p_val_adj), ]
      
      # Store top N genes
      top_n <- min(n_genes_per_cluster, nrow(markers))
      if (top_n > 0) {
        top_genes <- rownames(markers)[1:top_n]
        top_markers[[as.character(clusters[i])]] <- top_genes
        all_markers <- c(all_markers, top_genes)
      }
    }
  })
  
  # Remove duplicates in all_markers
  all_markers <- unique(all_markers)
  
  # If no markers found
  if (length(all_markers) == 0) {
    return(NULL)
  }
  
  # Create heatmap of expression across clusters
  # Get expression data
  expr_data <- AverageExpression(
    seurat_obj,
    features = all_markers,
    assays = "RNA",
    slot = "data",
    group.by = "seurat_clusters"
  )
  
  # Convert to matrix for heatmap
  expr_matrix <- expr_data$RNA
  
  # Scale the data
  expr_matrix_scaled <- t(scale(t(expr_matrix)))
  
  # Create heatmap
  pheatmap::pheatmap(
    expr_matrix_scaled,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 8,
    fontsize_col = 10,
    main = "Top Marker Genes by Cluster"
  )
}

# Function to generate a UMAP cluster report
generateClusterReport <- function(seurat_obj) {
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

# Function to find similar clusters
findSimilarClusters <- function(seurat_obj, resolution_range = seq(0.1, 1.0, by = 0.1)) {
  # Validate input
  if (!("pca" %in% names(seurat_obj@reductions))) {
    return(list(
      error = "No PCA data available for clustering"
    ))
  }
  
  results <- list()
  
  # Store original clustering if it exists
  original_clusters <- NULL
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    original_clusters <- seurat_obj$seurat_clusters
  }
  
  # Find neighbors using PCA
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  
  # Run clustering at different resolutions
  withProgress(message = 'Analyzing cluster stability', value = 0, {
    for (i in seq_along(resolution_range)) {
      res <- resolution_range[i]
      
      incProgress(1/length(resolution_range),
                  detail = paste("Clustering at resolution", res))
      
      seurat_obj <- FindClusters(seurat_obj, resolution = res)
      
      # Store results
      results[[paste0("res", res)]] <- list(
        resolution = res,
        n_clusters = length(unique(seurat_obj$seurat_clusters)),
        clusters = seurat_obj$seurat_clusters
      )
    }
  })
  
  # Restore original clustering if it existed
  if (!is.null(original_clusters)) {
    seurat_obj$seurat_clusters <- original_clusters
  }
  
  # Extract clustering information
  cluster_counts <- sapply(results, function(x) x$n_clusters)
  
  # Create summary data frame
  summary_df <- data.frame(
    Resolution = resolution_range,
    NumberOfClusters = cluster_counts
  )
  
  return(list(
    results = results,
    summary = summary_df
  ))
}