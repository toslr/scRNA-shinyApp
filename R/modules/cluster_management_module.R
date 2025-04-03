# R/modules/cluster_management_module.R

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
    
    # Add a loading lock to prevent competing updates
    is_loading_state <- reactiveVal(FALSE)
    ignore_input_changes <- reactiveVal(FALSE)
    
    # Use a dedicated reactiveVal for stable label storage
    # This is crucial to prevent flickering
    stable_labels <- reactiveVal(list())
    
    # Initialize state using reactiveValues for more stable reactivity
    state <- reactiveValues(
      cluster_labels = NULL,
      active_clusters = NULL,
      all_clusters = NULL,
      input_values = list(),
      last_update = NULL
    )
    
    # When clusters are initially loaded or change, initialize the stable labels
    observe({
      # Skip if we're in loading state
      if (is_loading_state()) {
        return(NULL)
      }
      
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
        
        # Get current labels
        current_labels <- stable_labels()
        
        # Create default labels while preserving existing ones
        new_labels <- initializeClusterLabels(available_clusters, current_labels)
        
        # Update both storage mechanisms atomically
        state$cluster_labels <- new_labels
        stable_labels(new_labels)
      }
      
      # Initialize active status if needed
      if (is.null(state$active_clusters) || 
          !all(as.character(available_clusters) %in% names(state$active_clusters))) {
        state$active_clusters <- initializeActiveStatus(available_clusters, state$active_clusters)
      }
    })
    
    # Setup UI for cluster rows - completely redesigned to prevent flickering
    output$clusterRows <- renderUI({
      # Skip if no clusters
      if (is.null(state$all_clusters) || length(state$all_clusters) == 0) {
        return(div(
          class = "alert alert-info",
          "No clusters available. Please run the clustering step first."
        ))
      }
      
      # Get the current stable labels - this is crucial for UI stability
      current_labels <- stable_labels()
      
      # Generate UI controls
      tagList(
        lapply(state$all_clusters, function(cluster) {
          cluster_key <- as.character(cluster)
          
          # First check if active
          is_active <- if (!is.null(state$active_clusters) && cluster_key %in% names(state$active_clusters)) {
            state$active_clusters[[cluster_key]]
          } else {
            TRUE  # Default to active
          }
          
          # Get label from stable source - this prevents flickering
          current_label <- if (!is.null(current_labels) && cluster_key %in% names(current_labels)) {
            current_labels[[cluster_key]]
          } else {
            paste("Cluster", cluster)
          }
          
          div(
            id = paste0("cluster-row-", cluster_key),
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
    
    # Track label changes with a dedicated observer for each cluster
    observe({
      req(state$all_clusters)
      
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      # Setup observers for each cluster's label input
      lapply(state$all_clusters, function(cluster) {
        cluster_key <- as.character(cluster)
        input_id <- paste0("label_", cluster)
        
        # Create an observer for this specific input
        observeEvent(input[[input_id]], {
          # Skip if in loading state
          if (ignore_input_changes()) {
            return(NULL)
          }
          
          # Get current label value from input
          current_value <- input[[input_id]]
          
          # Get current stable labels
          current_labels <- stable_labels()
          
          # Get previous value if it exists
          previous_value <- if (!is.null(current_labels) && cluster_key %in% names(current_labels)) {
            current_labels[[cluster_key]]
          } else {
            paste("Cluster", cluster)
          }
          
          # Only update if the value actually changed
          if (current_value != previous_value) {
            # Create a new copy of labels to modify
            new_labels <- current_labels
            new_labels[[cluster_key]] <- current_value
            
            # Update the stable storage
            stable_labels(new_labels)
            
            # Also update state for consistency
            state$cluster_labels[[cluster_key]] <- current_value
            
            # Log the change
            print(paste("Label for cluster", cluster_key, "changed from", 
                        previous_value, "to", current_value))
          }
        }, ignoreInit = TRUE)
      })
    })
    
    # Handle Select All checkbox
    observeEvent(input$selectAllClusters, {
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
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
      # Skip observer if we're loading a state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      req(state$all_clusters)
      
      for (cluster in state$all_clusters) {
        local({
          local_cluster <- cluster
          cluster_key <- as.character(local_cluster)
          input_id <- paste0("active_", local_cluster)
          
          if (!is.null(input[[input_id]])) {
            observeEvent(input[[input_id]], {
              # Skip if we're in loading state
              if (ignore_input_changes()) {
                return(NULL)
              }
              
              # Get current value before updating
              old_value <- FALSE
              if (!is.null(state$active_clusters) && cluster_key %in% names(state$active_clusters)) {
                old_value <- state$active_clusters[[cluster_key]]
              }
              
              # Only update and log if the value actually changed
              if (old_value != input[[input_id]]) {
                # Update active status for this cluster
                state$active_clusters[[cluster_key]] <- input[[input_id]]
                
                print(paste("Cluster", cluster_key, "active status changed from", 
                            old_value, "to", input[[input_id]]))
                
                # Update "Select All" checkbox
                all_active <- all(unlist(state$active_clusters))
                updateCheckboxInput(session, "selectAllClusters", value = all_active)
                
                state$last_update <- Sys.time()
              }
            }, ignoreInit = TRUE)
          }
        })
      }
    })
    
    # Update labels when the update button is clicked
    observeEvent(input$updateAllLabels, {
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      req(state$all_clusters)
      
      # Get current stable labels
      current_labels <- stable_labels()
      
      # Collect all input values
      for (cluster in state$all_clusters) {
        cluster_key <- as.character(cluster)
        input_id <- paste0("label_", cluster)
        
        # Only update if we have a value
        if (!is.null(input[[input_id]])) {
          current_labels[[cluster_key]] <- input[[input_id]]
        }
      }
      
      # Update stable labels in a single atomic operation
      stable_labels(current_labels)
      
      # Also update state for consistency
      state$cluster_labels <- current_labels
      
      showNotification("Cluster labels updated", type = "message")
    })
    
    # Function to get a comprehensive state snapshot for saving
    getFullState = function() {
      # Get the latest stable labels
      current_labels <- stable_labels()
      
      list(
        cluster_labels = current_labels,
        active_clusters = state$active_clusters,
        all_clusters = state$all_clusters
      )
    }
    
    # Function to restore state from a saved analysis
    setFullState = function(saved_state) {
      if (!is.null(saved_state)) {
        # Set loading flags to block other observers
        is_loading_state(TRUE)
        ignore_input_changes(TRUE)
        
        print("Setting cluster management full state")
        
        # Block all other observers and UI updates
        shinyjs::delay(100, {
          # Restore all_clusters if provided
          if (!is.null(saved_state$all_clusters) && length(saved_state$all_clusters) > 0) {
            state$all_clusters <- saved_state$all_clusters
          }
          
          # Restore labels atomically
          if (!is.null(saved_state$cluster_labels)) {
            print(paste("Restoring", length(saved_state$cluster_labels), "cluster labels"))
            
            # First update the stable storage
            stable_labels(saved_state$cluster_labels)
            
            # Then update state for consistency
            state$cluster_labels <- saved_state$cluster_labels
            
            # Print all labels being restored
            for (cluster_key in names(saved_state$cluster_labels)) {
              print(paste("Label for cluster", cluster_key, "=", saved_state$cluster_labels[[cluster_key]]))
            }
          }
          
          # Restore active status
          if (!is.null(saved_state$active_clusters)) {
            print(paste("Restoring", length(saved_state$active_clusters), "active status values:"))
            for (cluster_key in names(saved_state$active_clusters)) {
              print(paste(cluster_key, "->", saved_state$active_clusters[[cluster_key]]))
            }
            
            # Set active clusters
            state$active_clusters <- saved_state$active_clusters
          }
          
          # Update UI after a short delay
          shinyjs::delay(200, {
            # Update label inputs
            if (!is.null(saved_state$cluster_labels)) {
              for (cluster_key in names(saved_state$cluster_labels)) {
                label_value <- saved_state$cluster_labels[[cluster_key]]
                updateTextInput(session, paste0("label_", cluster_key), value = label_value)
              }
            }
            
            # Update checkboxes
            if (!is.null(saved_state$active_clusters)) {
              for (cluster_key in names(saved_state$active_clusters)) {
                is_active <- saved_state$active_clusters[[cluster_key]]
                updateCheckboxInput(session, paste0("active_", cluster_key), value = is_active)
              }
              
              # Update "Select All" checkbox
              all_active <- all(unlist(saved_state$active_clusters))
              updateCheckboxInput(session, "selectAllClusters", value = all_active)
            }
            
            # Force UI refresh
            state$last_update <- Sys.time()
            
            # Release the lock after a significant delay
            shinyjs::delay(500, {
              ignore_input_changes(FALSE)
              shinyjs::delay(200, {
                is_loading_state(FALSE)
                print("Cluster state loading complete - released locks")
              })
            })
          })
        })
      }
    }
    
    # Return reactive expressions
    list(
      getClusterLabels = reactive({ stable_labels() }),  # Use stable_labels instead of state
      getActiveStatus = reactive({ state$active_clusters }),
      getActiveClusterIds = reactive({
        if (is.null(state$active_clusters)) return(NULL)
        as.numeric(names(state$active_clusters[unlist(state$active_clusters) == TRUE]))
      }),
      getActiveClusterList = reactive({
        if (is.null(state$active_clusters)) return(NULL)
        names(state$active_clusters[unlist(state$active_clusters) == TRUE])
      }),
      updateLabels = function(new_labels) {
        # Update stable labels atomically
        stable_labels(new_labels)
        # Also update state for consistency
        state$cluster_labels <- new_labels
      },
      updateActiveStatus = function(new_active) {
        state$active_clusters <- new_active
      },
      last_update = reactive({ state$last_update }),
      getFullState = getFullState,
      setFullState = setFullState
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

#' @title Get Cluster Label
#' @description Retrieves the label for a specific cluster, handling reactive values
#'   and providing a fallback to "Cluster X" if no label is found.
#' @param cluster Cluster ID (numeric)
#' @param labels Named vector of cluster labels, can be a reactive value
#' @return A character string with the label for the specified cluster
#' @export
getClusterLabel <- function(cluster, labels) {
  # Safety check to handle reactive values
  if (is.reactive(labels)) {
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