# R/server/restore_ui.R

#' @title Restore QC UI State
#' @description Restores the QC UI state from the loaded data, including triggering the QC processing
#' @param qc_params List of QC parameters from saved state
#' @param session Shiny session object
#' @return None (called for side effects)
#' @export
restoreQcUI <- function(qc_params, session) {
  if (is.null(qc_params)) return(FALSE)
  
  # Update QC input values
  tryCatch({
    if (!is.null(qc_params$minFeature)) {
      updateNumericInput(session, "qc-minFeature", value = qc_params$minFeature)
    }
    if (!is.null(qc_params$maxFeature)) {
      updateNumericInput(session, "qc-maxFeature", value = qc_params$maxFeature)
    }
    if (!is.null(qc_params$maxMT)) {
      updateNumericInput(session, "qc-maxMT", value = qc_params$maxMT)
    }
    
    # Return whether QC was processed, but DON'T trigger the button
    return(!is.null(qc_params$qc_processed) && qc_params$qc_processed)
  }, error = function(e) {
    print(paste("Error restoring QC UI:", e$message))
  })
  
  return(FALSE)
}

#' @title Restore PCA UI State
#' @description Restores the dimension reduction UI state from loaded data
#' @param pca_params List of PCA parameters from saved state
#' @param session Shiny session object
#' @return None (called for side effects)
#' @export
restorePcaUI <- function(pca_params, session) {
  if (is.null(pca_params)) return(FALSE)
  
  tryCatch({
    # Update dimension input
    if (!is.null(pca_params$nDims)) {
      updateNumericInput(session, "dimRed-nDims", value = pca_params$nDims)
    }
    
    # Return whether dimensions were confirmed, but DON'T trigger the button
    return(!is.null(pca_params$dims_confirmed) && pca_params$dims_confirmed)
  }, error = function(e) {
    print(paste("Error restoring PCA UI:", e$message))
  })
  
  return(FALSE)
}

#' @title Restore Clustering UI State
#' @description Restores the clustering UI state from loaded data
#' @param clustering_params List of clustering parameters from saved state
#' @param session Shiny session object
#' @return Boolean indicating if clustering was restored
restoreClusteringUI <- function(clustering_params, session) {
  if (is.null(clustering_params)) return(FALSE)
  
  tryCatch({
    # Update resolution input
    if (!is.null(clustering_params$resolution)) {
      updateNumericInput(session, "dimRed-resolution", value = clustering_params$resolution)
    }
    
    # Return whether clustering was done, but DON'T trigger the button
    return(!is.null(clustering_params$clustering_done) && clustering_params$clustering_done)
  }, error = function(e) {
    print(paste("Error restoring clustering UI:", e$message))
  })
  
  return(FALSE)
}

#' @title Restore DE Analysis UI State
#' @description Restores the DE analysis UI state from loaded data
#' @param de_params List of DE parameters from saved state
#' @param session Shiny session object
#' @return Boolean indicating if DE analysis was restored
restoreDEAnalysisUI <- function(de_params, session) {
  if (is.null(de_params)) return(FALSE)
  
  tryCatch({
    # Update DE inputs
    if (!is.null(de_params$target_cluster_all)) {
      updateSelectInput(session, "de-targetClusterAll", selected = de_params$target_cluster_all)
    }
    
    if (!is.null(de_params$target_cluster1)) {
      updateSelectInput(session, "de-targetCluster1", selected = de_params$target_cluster1)
    }
    
    if (!is.null(de_params$target_cluster2)) {
      updateSelectInput(session, "de-targetCluster2", selected = de_params$target_cluster2)
    }
    
    if (!is.null(de_params$genes_per_cluster)) {
      updateNumericInput(session, "de-genesPerCluster", value = de_params$genes_per_cluster)
    }
    
    # Return whether DE was done, but DON'T trigger the buttons
    return(!is.null(de_params$de_analysis_done) && de_params$de_analysis_done)
  }, error = function(e) {
    print(paste("Error restoring DE UI:", e$message))
  })
  
  return(FALSE)
}