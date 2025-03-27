#' @title Perform Differential Expression Analysis
#' @description Performs differential expression analysis between clusters in a Seurat object.
#'   Can compare one cluster vs all others or one cluster vs a specific other cluster.
#' @param seurat_obj Seurat object containing clustered data
#' @param ident.1 Cluster ID to use as the reference group
#' @param ident.2 Optional cluster ID to compare against. If NULL, compares against all other clusters
#' @param active_clusters Optional vector of cluster IDs to include in the analysis
#' @return A data frame containing differential expression results with gene names
#' @export
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

#' @title Add Gene Names
#' @description Adds gene symbol names to differential expression results.
#'   Uses the gene mapping stored in the Seurat object to convert Ensembl IDs to gene symbols.
#' @param de_results Data frame containing differential expression results
#' @param seurat_obj Seurat object containing gene name mapping
#' @return The DE results data frame with an additional column for gene symbols
#' @export
addGeneNames <- function(de_results, seurat_obj) {
  gene_mapping <- seurat_obj@misc$gene_mapping
  if (!is.null(gene_mapping)) {
    de_results$gene <- gene_mapping[rownames(de_results)]
    de_results$gene[is.na(de_results$gene)] <- rownames(de_results)[is.na(de_results$gene)]
  } else {
    de_results$gene <- rownames(de_results)
  }
  return(de_results)
}