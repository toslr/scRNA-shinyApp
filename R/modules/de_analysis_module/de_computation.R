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