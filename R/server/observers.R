# R/observers.R

#' @title Setup Module Completion Observers
#' @description Sets up reactive observers that track the completion status of all analysis steps
#'   in the single-cell RNA-seq application. These observers update the steps_completed reactive
#'   values which drive the application's navigation and UI state.
#' @param steps_completed Reactive values object tracking completion status of each step
#' @param seurat_data Reactive expression containing the initial Seurat object
#' @param metadata_module Metadata module instance
#' @param processed_seurat Reactive expression containing the QC-processed Seurat object
#' @param clustered_seurat Reactive expression containing the clustered Seurat object
#' @param de_module Differential expression module instance
#' @param sample_management Sample management module instance
#' @param condition_management Condition management module instance
#' @return None (used for its side effects of setting up observers)
#' @export
setupObservers <- function(steps_completed, seurat_data, metadata_module, processed_seurat, 
                           clustered_seurat, de_module, sample_management, condition_management) {
  
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
  
  # Track condition management completion
  observe({
    if (!is.null(condition_management)) {
      condition_column <- condition_management$getConditionColumn()
      condition_active_status <- condition_management$getActiveStatus()
      
      if (!is.null(condition_column) && !is.null(condition_active_status)) {
        active_conditions <- names(condition_active_status[condition_active_status == TRUE])
        steps_completed$condition_management <- length(active_conditions) > 0
      } else {
        steps_completed$condition_management <- FALSE
      }
    } else {
      steps_completed$condition_management <- FALSE
    }
  })
}