# R/modules/de_analysis_module/cluster_utils.R

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

getClusterLabel <- function(cluster, labels) {
  # Safety check to handle reactive values
  if (is.function(labels)) {
    tryCatch({
      labels <- labels()
      print("Evaluated reactive labels in getClusterLabel")
    }, error = function(e) {
      print(paste("Error evaluating labels in getClusterLabel:", e$message))
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