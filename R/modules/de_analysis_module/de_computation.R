# R/modules/de_analysis_module/de_computation.R

performDEanalysis <- function(seurat_obj, ident.1, ident.2 = NULL, active_clusters = NULL) {
  # Handle subsetting for active clusters
  if (!is.null(active_clusters)) {
    active_cells <- seurat_obj$seurat_clusters %in% active_clusters
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
  }
  
  # Perform DE analysis
  de_results <- FindMarkers(seurat_obj,
                            ident.1 = ident.1,
                            ident.2 = ident.2,
                            min.pct = 0.25,
                            logfc.threshold = 0.25)
  
  # Sort by adjusted p-value
  de_results <- de_results[order(de_results$p_val_adj), ]
  
  return(addGeneNames(de_results, seurat_obj))
}

addGeneNames <- function(de_results, seurat_obj) {
  gene_mapping <- seurat_obj@misc$gene_mapping
  if (!is.null(gene_mapping)) {
    de_results$gene <- gene_mapping[rownames(de_results)]
    de_results$gene[is.na(de_results$gene)] <- rownames(de_results)[is.na(de_results$gene)]
  } else {
    warning("Gene mapping not found in Seurat object")
    de_results$gene <- rownames(de_results)
  }
  return(de_results)
}

# This function is not used directly anymore as we use a custom approach in the main module
computeGeneralDEGenes <- function(seurat_obj, active_clusters, genes_per_cluster) {
  # Check if we have any active clusters
  if (length(active_clusters) == 0) {
    warning("No active clusters provided")
    return(NULL)
  }
  
  # Initialize list to store results
  all_de_results <- list()
  
  # First, subset the Seurat object to only include cells from active clusters
  active_cells <- seurat_obj$seurat_clusters %in% active_clusters
  seurat_subset <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
  
  # Get the actual clusters present in the subset
  available_clusters <- sort(unique(seurat_subset$seurat_clusters))
  
  # Check if we have enough cells after subsetting
  if (length(available_clusters) == 0) {
    warning("No clusters available after subsetting")
    return(NULL)
  }
  
  # If only one cluster is active, use a different approach
  if (length(available_clusters) == 1) {
    # When only one cluster is active, find genes that are highly expressed in this cluster
    print("Only one cluster active - finding highly expressed genes")
    
    # Get normalized expression
    norm_data <- GetAssayData(seurat_subset, slot = "data")
    
    # Calculate average expression
    avg_expr <- rowMeans(norm_data)
    
    # Sort by expression and take top genes
    top_genes <- names(sort(avg_expr, decreasing = TRUE))[1:(genes_per_cluster*2)]
    
    return(top_genes)
  }
  
  # For each cluster, find markers against all other clusters
  withProgress(message = "Computing DE genes for clusters", detail = "Processing...", value = 0, {
    for(cluster in available_clusters) {
      incProgress(1/length(available_clusters), detail = paste("Processing cluster", cluster))
      
      # For each cluster, find markers against all other clusters
      de_results <- FindMarkers(seurat_subset,
                                ident.1 = cluster,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
      
      de_results <- addGeneNames(de_results, seurat_subset)
      
      # Handle case where fewer significant genes than requested
      num_genes_to_take <- min(genes_per_cluster, nrow(de_results))
      
      if (num_genes_to_take > 0) {
        top_genes <- rownames(de_results[order(de_results$p_val_adj), ])[1:num_genes_to_take]
        all_de_results[[as.character(cluster)]] <- top_genes
      } else {
        print(paste("No significant genes found for cluster", cluster))
      }
    }
  })
  
  # Combine all unique genes
  unique_genes <- unique(unlist(all_de_results))
  
  if (length(unique_genes) == 0) {
    warning("No differentially expressed genes found")
    # Return some genes to avoid empty heatmap
    norm_data <- GetAssayData(seurat_subset, slot = "data")
    avg_expr <- rowMeans(norm_data)
    unique_genes <- names(sort(avg_expr, decreasing = TRUE))[1:(genes_per_cluster*length(available_clusters))]
  }
  
  print(paste("Found", length(unique_genes), "unique genes across all clusters"))
  return(unique_genes)
}