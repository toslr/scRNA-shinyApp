# R/observers.R

setupObservers <- function(steps_completed, seurat_data, metadata_module, processed_seurat, 
                           clustered_seurat, de_module, sample_management) {
  
  # Track metadata completion
  observe({
    print("Checking metadata completion")
    steps_completed$metadata <- !is.null(metadata_module$getMetadata())
  })
  
  # Track data input completion
  observe({
    steps_completed$data_input <- !is.null(seurat_data())
  })
  
  # Track QC completion
  observe({
    steps_completed$qc <- !is.null(processed_seurat())
  })
  
  # Track dimension reduction and clustering completion
  observe({
    if (!is.null(clustered_seurat())) {
      steps_completed$dimred <- any(c("umap2d", "umap3d", "umap") %in% names(clustered_seurat()@reductions))
      steps_completed$clustering <- "seurat_clusters" %in% colnames(clustered_seurat()@meta.data)
    }
  })
  
  # Track DE completion
  observe({
    steps_completed$de <- !is.null(de_module$status()) && de_module$status() == "completed"
  })
  
  # Track sample management completion
  observe({
    sample_active_status <- sample_management$getActiveStatus()
    if (!is.null(sample_active_status)) {
      active_samples <- names(sample_active_status[sample_active_status == TRUE])
      steps_completed$sample_management <- length(active_samples) > 0
    } else {
      steps_completed$sample_management <- FALSE
    }
  })
  
}