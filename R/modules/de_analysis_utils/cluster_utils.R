#' @title Initialize Cluster Labels
#' @description Creates default labels for clusters with format "Cluster X" while 
#'   preserving any existing labels from the current_labels parameter.
#' @param clusters Numeric vector of cluster IDs
#' @param current_labels Optional named vector of existing cluster labels
#' @return A named vector mapping cluster IDs to their labels
#' @export
initializeClusterLabels <- function(clusters, current_labels = NULL) {
  new_labels <- setNames(
    paste("Cluster", clusters), 
    as.character(clusters)
  )
  
  # Merge with existing labels if they exist
  if (!is.null(current_labels)) {
    existing_clusters <- names(current_labels)
    new_labels[existing_clusters] <- current_labels[existing_clusters]
  }
  
  return(new_labels)
}

#' @title Initialize Active Status
#' @description Creates a named vector indicating which clusters are active (TRUE/FALSE),
#'   with the default being all clusters active. Preserves existing statuses if provided.
#' @param clusters Numeric vector of cluster IDs
#' @param current_active Optional named vector of existing active statuses
#' @return A named logical vector indicating which clusters are active
#' @export
initializeActiveStatus <- function(clusters, current_active = NULL) {
  new_active <- setNames(
    rep(TRUE, length(clusters)),
    as.character(clusters)
  )
  
  if (!is.null(current_active)) {
    existing_clusters <- names(current_active)
    new_active[existing_clusters] <- current_active[existing_clusters]
  }
  
  return(new_active)
}

#' @title Get Cluster Label
#' @description Retrieves the label for a specific cluster, handling reactive values
#'   and providing a fallback to "Cluster X" if no label is found.
#' @param cluster Cluster ID (numeric)
#' @param labels Named vector of cluster labels, can be a reactive value
#' @return A character string with the label for the specified cluster
#' @export
getClusterLabel <- function(cluster, labels) {
  # Safety check to handle reactive values
  if (is.function(labels)) {
    tryCatch({
      labels <- labels()
    }, error = function(e) {
      return(paste("Cluster", cluster))
    })
  }
  
  # Check if labels is NULL or not a list/vector
  if (is.null(labels) || !is.vector(labels)) {
    return(paste("Cluster", cluster))
  }
  
  cluster_key <- as.character(cluster)
  
  if (cluster_key %in% names(labels) && !is.null(labels[[cluster_key]])) {
    labels[[cluster_key]]
  } else {
    paste("Cluster", cluster)
  }
}