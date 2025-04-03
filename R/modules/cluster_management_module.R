#' @title Cluster Management Module UI
#' @description Creates the UI for the cluster management module which allows users
#'   to toggle clusters and edit cluster labels.
#' @param id The module ID
#' @return A Shiny UI element containing the cluster management interface
#' @export
clusterManagementUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("clusterRows"))
  )
}

#' @title Cluster Management Module Server
#' @description Server logic for the cluster management module which handles cluster
#'   activation/deactivation and label editing.
#' @param id The module ID
#' @param clustered_seurat Reactive expression containing the clustered Seurat object
#' @return A list of reactive expressions for accessing cluster labels and active status
#' @export
clusterManagementServer <- function(id, clustered_seurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize state using reactiveValues for more stable reactivity
    state <- reactiveValues(
      cluster_labels = NULL,
      active_clusters = NULL,
      all_clusters = NULL,
      input_values = list(),
      last_update = NULL
    )
    
    # Get clusters from Seurat object and initialize state
    observe({
      available_clusters <- getAvailableClusters(clustered_seurat())
      
      # Skip if no clusters available
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      # Update all_clusters
      state$all_clusters <- available_clusters
      
      # Initialize labels if needed
      if (is.null(state$cluster_labels) || 
          !all(as.character(available_clusters) %in% names(state$cluster_labels))) {
        new_labels <- initializeClusterLabels(available_clusters, state$cluster_labels)
        state$cluster_labels <- new_labels
      }
      
      # Initialize active status if needed
      if (is.null(state$active_clusters) || 
          !all(as.character(available_clusters) %in% names(state$active_clusters))) {
        state$active_clusters <- initializeActiveStatus(available_clusters, state$active_clusters)
      }
    })
    
    # Setup UI for cluster rows
    output$clusterRows <- renderUI({
      # Skip if no clusters
      if (is.null(state$all_clusters) || length(state$all_clusters) == 0) {
        return(div(
          class = "alert alert-info",
          "No clusters available. Please run the clustering step first."
        ))
      }
      
      # Generate UI controls
      tagList(
        lapply(state$all_clusters, function(cluster) {
          cluster_key <- as.character(cluster)
          is_active <- if (!is.null(state$active_clusters) && cluster_key %in% names(state$active_clusters)) {
            state$active_clusters[[cluster_key]]
          } else {
            TRUE  # Default to active
          }
          
          # Get current label with fallback options
          current_label <- NULL
          
          # First check if we have a stored input value
          input_id <- paste0("label_", cluster)
          if (!is.null(state$input_values[[input_id]])) {
            current_label <- state$input_values[[input_id]]
          } 
          # Then check saved labels
          else if (!is.null(state$cluster_labels) && cluster_key %in% names(state$cluster_labels)) {
            current_label <- state$cluster_labels[[cluster_key]]
          } 
          # Fallback to default cluster name
          else {
            current_label <- paste("Cluster", cluster)
          }
          
          div(
            style = paste0(
              "margin-bottom: 10px; padding: 8px; border-radius: 4px; ",
              if (is_active) "background-color: #f8f9fa;" else "background-color: #e9ecef; opacity: 0.8;"
            ),
            fluidRow(
              column(2,
                     checkboxInput(ns(paste0("active_", cluster)), 
                                   label = paste("", cluster),
                                   value = is_active)
              ),
              column(10, 
                     tags$div(
                       textInput(ns(paste0("label_", cluster)),
                                 label = NULL,
                                 value = current_label)
                     )
              )
            )
          )
        })
      )
    })
    
    # Track input values for each cluster's label
    observe({
      req(state$all_clusters)
      
      lapply(state$all_clusters, function(cluster) {
        input_id <- paste0("label_", cluster)
        
        # Update input_values when input changes
        if (!is.null(input[[input_id]])) {
          # This isolate prevents a circular dependency
          isolate({
            state$input_values[[input_id]] <- input[[input_id]]
          })
        }
      })
    })
    
    # Handle Select All checkbox
    observeEvent(input$selectAllClusters, {
      req(state$all_clusters)
      
      # Create a new active_clusters list with all clusters set to the checkbox value
      new_active <- setNames(
        rep(input$selectAllClusters, length(state$all_clusters)),
        as.character(state$all_clusters)
      )
      
      # Update state
      state$active_clusters <- new_active
      
      # Update individual checkboxes to match
      for (cluster in state$all_clusters) {
        updateCheckboxInput(session, paste0("active_", cluster), value = input$selectAllClusters)
      }
      
      state$last_update <- Sys.time()
    })
    
    # Handle individual cluster activation
    observe({
      req(state$all_clusters)
      
      for (cluster in state$all_clusters) {
        local({
          local_cluster <- cluster
          cluster_key <- as.character(local_cluster)
          input_id <- paste0("active_", local_cluster)
          
          if (!is.null(input[[input_id]])) {
            observeEvent(input[[input_id]], {
              old_value <- state$active_clusters[[cluster_key]]
              # Update active status for this cluster
              state$active_clusters[[cluster_key]] <- input[[input_id]]
              
              # Update "Select All" checkbox
              all_active <- all(unlist(state$active_clusters))
              updateCheckboxInput(session, "selectAllClusters", value = all_active)
              
              state$last_update <- Sys.time()
            }, ignoreInit = TRUE)
          }
        })
      }
    })
    
    # Update labels when the update button is clicked
    observeEvent(input$updateAllLabels, {
      req(state$all_clusters)
      
      # Collect all input values
      for (cluster in state$all_clusters) {
        cluster_key <- as.character(cluster)
        input_id <- paste0("label_", cluster)
        
        # Only update if we have a value
        if (!is.null(input[[input_id]])) {
          state$cluster_labels[[cluster_key]] <- input[[input_id]]
          state$input_values[[input_id]] <- input[[input_id]]
        }
      }
      
      showNotification("Cluster labels updated", type = "message")
    })
    
    # Return reactive expressions
    list(
      getClusterLabels = reactive({ state$cluster_labels }),
      getActiveStatus = reactive({ state$active_clusters }),
      getFullState = function() {
        list(
          cluster_labels = state$cluster_labels,
          active_clusters = state$active_clusters,
          all_clusters = state$all_clusters,
          input_values = state$input_values
        )
      },
      setFullState = function(saved_state) {
        if (!is.null(saved_state)) {
          if (!is.null(saved_state$cluster_labels)) {
            state$cluster_labels <- saved_state$cluster_labels
          }
          if (!is.null(saved_state$active_clusters)) {
            state$active_clusters <- saved_state$active_clusters
          }
          if (!is.null(saved_state$input_values)) {
            state$input_values <- saved_state$input_values
          }
        }
      },
      getActiveClusterIds = reactive({
        if (is.null(state$active_clusters)) return(NULL)
        as.numeric(names(state$active_clusters[unlist(state$active_clusters) == TRUE]))
      }),
      getActiveClusterList = reactive({
        if (is.null(state$active_clusters)) return(NULL)
        names(state$active_clusters[unlist(state$active_clusters) == TRUE])
      }),
      updateLabels = function(new_labels) {
        state$cluster_labels <- new_labels
        for (cluster_key in names(new_labels)) {
          input_id <- paste0("label_", gsub("^cluster_", "", cluster_key))
          state$input_values[[input_id]] <- new_labels[[cluster_key]]
        }
      },
      updateActiveStatus = function(new_active) {
        state$active_clusters <- new_active
      },
      last_update = reactive({ state$last_update })
    )
  })
}

#' @title Get Available Clusters
#' @description Helper function that extracts unique cluster IDs from a Seurat object.
#' @param seurat_obj Seurat object containing clustered data
#' @return Vector of unique cluster IDs, or NULL if no clusters are found
#' @keywords internal
getAvailableClusters <- function(seurat_obj) {
  if (is.null(seurat_obj) || !"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }
  sort(unique(seurat_obj$seurat_clusters))
}

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
    for (cluster in existing_clusters) {
      if (cluster %in% names(new_labels)) {
        new_labels[cluster] <- current_labels[cluster]
      }
    }
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
    for (cluster in existing_clusters) {
      if (cluster %in% names(new_active)) {
        new_active[cluster] <- current_active[cluster]
      }
    }
  }
  
  return(new_active)
}