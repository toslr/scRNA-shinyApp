# R/observers.R

setupObservers <- function(steps_completed, seurat_data, processed_seurat, 
                           clustered_seurat, de_module) {
  # Track step completion status
  observeEvent(seurat_data(), {
    steps_completed$data_input <- !is.null(seurat_data())
  })
  
  observeEvent(processed_seurat(), {
    steps_completed$qc <- !is.null(processed_seurat())
  })
  
  observeEvent(clustered_seurat(), {
    if (!is.null(clustered_seurat())) {
      steps_completed$dimred <- "umap" %in% names(clustered_seurat()@reductions)
      steps_completed$clustering <- "seurat_clusters" %in% colnames(clustered_seurat()@meta.data)
    }
  })
  
  observeEvent(de_module$status(), {
    steps_completed$de <- !is.null(de_module$status()) && de_module$status() == "completed"
  })
}