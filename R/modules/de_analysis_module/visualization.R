# R/modules/de_analysis_module/visualization.R

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

createExpressionHeatmap <- function(seurat_obj, top_genes, cluster_labels) {
  # Get expression data
  expr_data <- GetAssayData(seurat_obj, slot = "data")[top_genes, ]
  
  # Calculate cluster means
  clusters <- seurat_obj$seurat_clusters
  unique_clusters <- sort(unique(clusters))

  cluster_means <- sapply(unique_clusters, function(clust) {
    rowMeans(expr_data[, clusters == clust, drop = FALSE])
  })

  colnames(cluster_means) <- sapply(unique_clusters, function(x) cluster_labels[[as.character(x)]])
  
  # Scale data
  scaled_data <- t(scale(t(cluster_means)))
  
  # Get gene labels
  gene_mapping <- seurat_obj@misc$gene_mapping
  gene_labels <- if (!is.null(gene_mapping)) {
    gene_mapping[rownames(scaled_data)]
  } else {
    rownames(scaled_data)
  }
  gene_labels[is.na(gene_labels)] <- rownames(scaled_data)[is.na(gene_labels)]

  pheatmap(scaled_data,
           labels_row = gene_labels,
           labels_col = colnames(cluster_means),
           main = paste("Top", length(top_genes), "DE Genes"),
           angle_col = 45,
           fontsize_row = 10,
           cluster_cols = FALSE,
           treeheight_row = 0,
           draw = TRUE)
}