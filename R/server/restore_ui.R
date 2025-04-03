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
    
    # If QC was processed, we need to trigger the processing button click
    if (!is.null(qc_params$qc_processed) && qc_params$qc_processed) {
      # Delay slightly to let the inputs update first
      shinyjs::delay(300, {
        # Programmatically click the process button
        shinyjs::click("qc-processSeurat")
      })
      return(TRUE)
    }
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
    
    # If dimensions were confirmed, trigger the confirmation
    if (!is.null(pca_params$dims_confirmed) && pca_params$dims_confirmed) {
      # Delay slightly to let the inputs update first
      shinyjs::delay(600, {
        # Programmatically click the confirm dimensions button
        shinyjs::click("dimRed-confirmDims")
      })
      return(TRUE)
    }
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
    
    # If clustering was done, trigger the button click
    if (!is.null(clustering_params$clustering_done) && clustering_params$clustering_done) {
      shinyjs::delay(300, {
        shinyjs::click("dimRed-runClustering")
      })
      return(TRUE)
    }
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
    
    # If DE analysis was done, trigger the appropriate button based on analysis type
    if (!is.null(de_params$de_analysis_done) && de_params$de_analysis_done) {
      if (!is.null(de_params$analysis_type)) {
        if (de_params$analysis_type == "one_vs_all") {
          shinyjs::delay(300, {
            shinyjs::click("de-runDEAll")
          })
          return(TRUE)
        } else if (de_params$analysis_type == "pairwise") {
          shinyjs::delay(300, {
            shinyjs::click("de-runDEPair")
          })
          return(TRUE)
        } else if (de_params$analysis_type == "general_heatmap") {
          shinyjs::delay(300, {
            shinyjs::click("de-runGeneralHeatmap")
          })
          return(TRUE)
        }
      }
    }
  }, error = function(e) {
    print(paste("Error restoring DE UI:", e$message))
  })
  
  return(FALSE)
}