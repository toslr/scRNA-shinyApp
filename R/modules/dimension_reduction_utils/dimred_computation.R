# R/modules/dimension_reduction_utils/dimred_computation.R

#' @title Find Elbow Point in PCA Variance
#' @description Identifies the optimal dimension cutoff based on the rate of change 
#'   in variance explained by principal components.
#' @param x Vector of dimensions
#' @param y Vector of variance explained
#' @return Integer representing the optimal dimension cutoff
#' @keywords internal
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

#' @title Compute Suggested Dimensions from PCA
#' @description Automatically determines the optimal number of principal components
#'   to use for downstream analysis based on the elbow point of variance explained.
#' @param seurat_obj Seurat object with PCA run
#' @return Integer suggesting optimal number of PCs to use
#' @export
compute_suggested_dims <- function(seurat_obj) {
  # Ensure PCA has been run
  if (!("pca" %in% names(seurat_obj@reductions))) {
    stop("PCA has not been run on this Seurat object")
  }
  
  # Get variance explained
  pca_data <- Embeddings(seurat_obj, "pca")
  stdev <- Stdev(seurat_obj[["pca"]])
  var_explained <- stdev^2 / sum(stdev^2)
  
  # Find elbow point
  dims <- 1:length(var_explained)
  find_elbow(dims, var_explained)
}

#' @title Run UMAP with Specified Parameters
#' @description Runs UMAP on a Seurat object with the specified parameters, creating
#'   either 2D or 3D embeddings.
#' @param seurat_obj Seurat object with PCA run
#' @param n_dims Integer number of PCA dimensions to use
#' @param n_components Integer number of UMAP dimensions (2 or 3)
#' @param reduction.name String name for stored reduction
#' @param reduction.key String key prefix for dimension names
#' @param ... Additional parameters to pass to RunUMAP
#' @return Seurat object with added UMAP reduction
#' @keywords internal
run_umap <- function(seurat_obj, n_dims, n_components = 2, 
                     reduction.name = "umap", reduction.key = "UMAP_", ...) {
  
  # Validate inputs
  if (!("pca" %in% names(seurat_obj@reductions))) {
    stop("PCA has not been run on this Seurat object")
  }
  
  if (n_dims > ncol(Embeddings(seurat_obj, "pca"))) {
    warning("Requested more dimensions than available. Using maximum available.")
    n_dims <- ncol(Embeddings(seurat_obj, "pca"))
  }
  
  # Run UMAP
  RunUMAP(seurat_obj, 
          dims = 1:n_dims, 
          n.components = n_components, 
          reduction.name = reduction.name,
          reduction.key = reduction.key,
          ...)
}

#' @title Process Both 2D and 3D UMAPs
#' @description Creates both 2D and 3D UMAP embeddings for a Seurat object
#'   to enable different visualization options.
#' @param seurat_obj Seurat object with PCA run
#' @param n_dims Integer number of PCA dimensions to use
#' @param ... Additional parameters to pass to RunUMAP
#' @return Seurat object with both 2D and 3D UMAP reductions
#' @export
process_umaps <- function(seurat_obj, n_dims, ...) {
  # Run 2D UMAP
  seurat_2d <- run_umap(seurat_obj, n_dims, n_components = 2, 
                        reduction.name = "umap2d", reduction.key = "UMAP2D_", ...)
  
  # Run 3D UMAP
  seurat_3d <- run_umap(seurat_2d, n_dims, n_components = 3,
                        reduction.name = "umap3d", reduction.key = "UMAP3D_", ...)
  
  # Also add a standard umap reduction for compatibility
  if (!("umap" %in% names(seurat_3d@reductions))) {
    seurat_3d@reductions$umap <- seurat_3d@reductions$umap2d
  }
  
  return(seurat_3d)
}

#' @title Run Clustering on Seurat Object
#' @description Performs neighbor finding and clustering on a Seurat object with the 
#'   specified parameters.
#' @param seurat_obj Seurat object with PCA run
#' @param n_dims Integer number of PCA dimensions to use
#' @param resolution Numeric clustering resolution
#' @param ... Additional parameters to pass to FindClusters
#' @return Seurat object with clustering results
#' @export
run_clustering <- function(seurat_obj, n_dims, resolution = 0.5, ...) {
  # Find neighbors first
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims)
  
  # Find clusters
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, ...)
  
  return(seurat_obj)
}

#' @title Get Cluster Statistics
#' @description Extracts statistics about clusters from a Seurat object, including
#'   cell counts, percentages, and sample distribution.
#' @param seurat_obj Seurat object with clustering results
#' @return Data frame with cluster statistics
#' @export
get_cluster_stats <- function(seurat_obj) {
  # Check if clustering has been performed
  if (!("seurat_clusters" %in% colnames(seurat_obj@meta.data))) {
    return(data.frame(Error = "No clustering data available"))
  }
  
  # Get cluster identities
  clusters <- seurat_obj$seurat_clusters
  
  # Count cells per cluster
  cell_counts <- table(clusters)
  
  # Calculate percentage of cells per cluster
  cell_percentages <- prop.table(cell_counts) * 100
  
  # Create output data frame
  stats <- data.frame(
    Cluster = names(cell_counts),
    CellCount = as.numeric(cell_counts),
    Percentage = round(as.numeric(cell_percentages), 2)
  )
  
  # If sample information is available, calculate distribution by sample
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    # Create a contingency table of clusters by sample
    cluster_by_sample <- table(seurat_obj$seurat_clusters, seurat_obj$sample)
    
    # For each sample, add a column to the stats data frame
    for (sample_name in colnames(cluster_by_sample)) {
      sample_counts <- cluster_by_sample[, sample_name]
      stats[[paste0("Count_", sample_name)]] <- sample_counts
      
      # Calculate percentage within each cluster
      for (i in 1:nrow(stats)) {
        cluster_id <- stats$Cluster[i]
        cluster_total <- stats$CellCount[i]
        sample_count <- sample_counts[cluster_id]
        stats[[paste0("Pct_", sample_name)]][i] <- round((sample_count / cluster_total) * 100, 2)
      }
    }
  }
  
  return(stats)
}

#' @title Find Similar Clusters Across Resolutions
#' @description Runs clustering at different resolutions to help identify optimal 
#'   clustering parameters.
#' @param seurat_obj Seurat object with PCA run
#' @param resolution_range Vector of resolutions to test
#' @param dims Integer vector of PCA dimensions to use
#' @return List containing results and summary of clustering at different resolutions
#' @export
find_similar_clusters <- function(seurat_obj, resolution_range = seq(0.1, 1.0, by = 0.1), dims = 1:20) {
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
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  
  # Run clustering at different resolutions
  for (i in seq_along(resolution_range)) {
    res <- resolution_range[i]
    
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
    
    # Store results
    results[[paste0("res", res)]] <- list(
      resolution = res,
      n_clusters = length(unique(seurat_obj$seurat_clusters)),
      clusters = seurat_obj$seurat_clusters
    )
  }
  
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

#' @title Filter Seurat Object by Samples
#' @description Filters a Seurat object to include only cells from specified samples.
#' @param seurat_obj Seurat object
#' @param active_samples Vector of sample names to keep
#' @return Filtered Seurat object
#' @export
filter_by_samples <- function(seurat_obj, active_samples) {
  # Check if we have a valid Seurat object and active samples
  if (is.null(seurat_obj) || is.null(active_samples) || length(active_samples) == 0) {
    return(seurat_obj)
  }
  
  # Check if sample column exists
  if (!"sample" %in% colnames(seurat_obj@meta.data)) {
    warning("Sample column not found in Seurat object metadata")
    return(seurat_obj)
  }
  
  # Get cells that belong to active samples
  cells_to_keep <- seurat_obj$sample %in% active_samples
  
  # Only subset if necessary
  if (all(cells_to_keep)) {
    return(seurat_obj)  # All cells are already active
  } else if (!any(cells_to_keep)) {
    warning("No cells match the active samples selection")
    return(seurat_obj)  # Return original to avoid empty object
  }
  
  # Subset the Seurat object
  subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
}

#' @title Filter Seurat Object by Conditions
#' @description Filters a Seurat object to include only cells with specified condition values.
#' @param seurat_obj Seurat object
#' @param condition_column String name of the condition column
#' @param active_conditions Vector of condition values to keep
#' @return Filtered Seurat object
#' @export
filter_by_conditions <- function(seurat_obj, condition_column, active_conditions) {
  # Check if we have valid inputs
  if (is.null(seurat_obj) || is.null(condition_column) || 
      is.null(active_conditions) || length(active_conditions) == 0) {
    return(seurat_obj)
  }
  
  # Check if condition column exists
  if (!(condition_column %in% colnames(seurat_obj@meta.data))) {
    warning(paste("Condition column", condition_column, "not found in Seurat object metadata"))
    return(seurat_obj)
  }
  
  # Get cells that match active conditions
  cells_to_keep <- seurat_obj@meta.data[[condition_column]] %in% active_conditions
  
  # Only subset if necessary
  if (all(cells_to_keep)) {
    return(seurat_obj)  # All cells are already active
  } else if (!any(cells_to_keep)) {
    warning("No cells match the active conditions selection")
    return(seurat_obj)  # Return original to avoid empty object
  }
  
  # Subset the Seurat object
  subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
}

#' @title Create Cluster Marker Heatmap
#' @description Creates a heatmap of marker genes for each cluster to help
#'   with identification and interpretation of clusters.
#' @param seurat_obj Seurat object with clusters
#' @param n_genes_per_cluster Integer number of genes per cluster
#' @return A pheatmap object
#' @export
create_cluster_marker_heatmap <- function(seurat_obj, n_genes_per_cluster = 5) {
  # Check if clustering has been performed
  if (!("seurat_clusters" %in% colnames(seurat_obj@meta.data))) {
    return(NULL)
  }
  
  # List to store top markers
  top_markers <- list()
  all_markers <- c()
  
  # For each cluster, find markers
  clusters <- sort(unique(seurat_obj$seurat_clusters))
  
  # Find markers for each cluster
  for (i in seq_along(clusters)) {
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