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

computeGeneralDEGenes <- function(seurat_obj, active_clusters, genes_per_cluster) {
  # Initialize list to store results
  all_de_results <- list()
  
  # First, subset the Seurat object to only include cells from active clusters
  active_cells <- seurat_obj$seurat_clusters %in% active_clusters
  seurat_subset <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
  
  # Get the actual clusters present in the subset
  available_clusters <- sort(unique(seurat_subset$seurat_clusters))
  withProgress(message = "Computing DE genes for clusters", detail = "Processing...", value = 0, {
    for(cluster in available_clusters) {
      incProgress(1/length(available_clusters), detail = paste("Processing cluster", cluster))
      de_results <- FindMarkers(seurat_subset,
                                ident.1 = cluster,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
      de_results <- addGeneNames(de_results, seurat_subset)
      top_genes <- rownames(de_results[order(de_results$p_val_adj), ])[1:genes_per_cluster]
      all_de_results[[as.character(cluster)]] <- top_genes
    }
  })
  
  # Combine all unique genes
  unique_genes <- unique(unlist(all_de_results))
  print(paste("Found", length(unique_genes), "unique genes across all clusters"))
  return(unique_genes)
}