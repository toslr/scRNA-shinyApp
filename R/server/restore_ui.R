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