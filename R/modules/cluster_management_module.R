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
    
    # Add loading locks to prevent competing updates
    is_loading_state <- reactiveVal(FALSE)
    ignore_input_changes <- reactiveVal(FALSE)
    
    # Add flags to control update flow and prevent cascades
    updating_from_checkbox <- reactiveVal(FALSE)
    updating_select_all <- reactiveVal(FALSE)
    
    # Use a dedicated reactiveVal for stable label storage
    # This is crucial to prevent flickering
    stable_labels <- reactiveVal(list())
    
    # Initialize state using reactiveValues for more stable reactivity
    state <- reactiveValues(
      cluster_labels = NULL,
      active_clusters = NULL,
      all_clusters = NULL,
      temp_labels = NULL,
      input_values = list(),
      last_update = NULL
    )
    
    # Flag to control UI updates
    inhibit_ui_updates <- reactiveVal(FALSE)
    
    # Add debouncing for UI updates to prevent multiple rapid redraws
    observe({
      # This will throttle UI updates when typing
      if (inhibit_ui_updates()) {
        # If inhibited, don't trigger UI updates for a short time
        invalidateLater(300) # Wait 300ms before allowing updates again
        inhibit_ui_updates(FALSE)
      }
    })
    
    # Create a deep copy of a list
    deep_copy <- function(list_obj) {
      if (is.null(list_obj)) return(list())
      result <- list()
      for (name in names(list_obj)) {
        result[[name]] <- list_obj[[name]]
      }
      return(result)
    }
    
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
        
        # Also initialize temp_labels if needed
        if (is.null(state$temp_labels)) {
          state$temp_labels <- list()
        }
      }
      
      # Initialize active status if needed
      if (is.null(state$active_clusters) || 
          !all(as.character(available_clusters) %in% names(state$active_clusters))) {
        state$active_clusters <- initializeActiveStatus(available_clusters, state$active_clusters)
      }
    })
    
    # Add JavaScript handler for individual checkboxes
    session$sendCustomMessage(type = "initializeClusterCheckboxes", message = list())
    
    # Setup UI for cluster rows
    output$clusterRows <- renderUI({
      # Skip if no clusters
      if (is.null(state$all_clusters) || length(state$all_clusters) == 0) {
        return(div(
          class = "alert alert-info",
          "No clusters available. Please run the clustering step first."
        ))
      }
      
      # Use isolate here to prevent UI rebuilding during typing
      isolate({
        # Get the current stable labels - this is crucial for UI stability
        current_labels <- stable_labels()
        current_active <- state$active_clusters
        current_temp <- state$temp_labels
        
        # Generate UI controls - no select all button here
        div(
          style = "max-height: 400px; overflow-y: auto; border: 1px solid #e3e3e3; padding: 5px; border-radius: 4px;",
          lapply(state$all_clusters, function(cluster) {
            cluster_key <- as.character(cluster)
            
            # First check if active
            is_active <- if (!is.null(current_active) && cluster_key %in% names(current_active)) {
              current_active[[cluster_key]]
            } else {
              TRUE  # Default to active
            }
            
            # Get label value to display - use temp label if available, otherwise stable label
            display_value <- if (!is.null(current_temp) && cluster_key %in% names(current_temp)) {
              current_temp[[cluster_key]]
            } else if (!is.null(current_labels) && cluster_key %in% names(current_labels)) {
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
                       tags$div(
                         style = "margin-top: 5px;",
                         tags$div(
                           class = "checkbox",
                           tags$label(
                             tags$input(
                               type = "checkbox",
                               id = ns(paste0("active_", cluster)),
                               checked = if(is_active) "checked" else NULL,
                               class = "cluster-management-checkbox",
                               `data-cluster` = cluster  # Add data attribute for cluster ID
                             )
                           )
                         )
                       )
                ),
                column(10, 
                       tags$div(
                         textInput(ns(paste0("label_", cluster)),
                                   label = NULL,
                                   value = display_value,
                                   width = "100%")
                       )
                )
              )
            )
          })
        )
      })
    })
    
    # Track label changes with a dedicated observer
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
          
          # Store the input value in a state variable without triggering UI updates
          isolate({
            if (is.null(state$temp_labels)) {
              state$temp_labels <- list()
            }
            
            # Update the temporary value
            state$temp_labels[[cluster_key]] <- input[[input_id]]
            
            # Print for debugging
            print(paste("Temp label for cluster", cluster_key, "set to:", input[[input_id]]))
          })
          
        }, ignoreInit = TRUE, ignoreNULL = TRUE)
      })
    })
    
    # Handle Select All checkbox
    observeEvent(input$selectAllClusters, {
      # Skip if we're already updating from a checkbox change
      if (updating_from_checkbox()) {
        return(NULL)
      }
      
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      # Mark that we're updating from Select All
      updating_select_all(TRUE)
      
      req(state$all_clusters)
      
      # Create a new active_clusters list with all clusters set to the checkbox value
      new_active <- setNames(
        rep(input$selectAllClusters, length(state$all_clusters)),
        as.character(state$all_clusters)
      )
      
      # Block individual checkbox observers
      ignore_input_changes(TRUE)
      
      # Update state
      state$active_clusters <- new_active
      
      # Update individual checkboxes using JavaScript
      selected_clusters <- if(input$selectAllClusters) as.character(state$all_clusters) else character(0)
      session$sendCustomMessage(
        type = "updateClusterCheckboxes",
        message = list(clusters = selected_clusters)
      )
      
      # Now it's safe to allow changes - slight delay to ensure UI updates first
      shinyjs::delay(200, {
        ignore_input_changes(FALSE)
        updating_select_all(FALSE)
        
        # Trigger state update
        state$last_update <- Sys.time()
      })
    }, ignoreInit = TRUE, priority = 10)  # Higher priority
    
    # Handler for individual cluster checkbox updates from JavaScript
    observeEvent(input$cluster_checkbox_changed, {
      # Skip if in loading state or during Select All update
      if (ignore_input_changes() || updating_select_all()) {
        return(NULL)
      }
      
      # Extract the cluster ID and checked state from the input
      cluster_id <- input$cluster_checkbox_changed$cluster
      is_checked <- input$cluster_checkbox_changed$checked
      
      # Mark that we're updating from a checkbox
      updating_from_checkbox(TRUE)
      
      # Get current active status - create a complete copy to avoid reference issues
      current_active <- deep_copy(state$active_clusters)
      
      # Get previous value with default
      old_value <- if (!is.null(current_active) && cluster_id %in% names(current_active)) {
        current_active[[cluster_id]]
      } else {
        TRUE
      }
      
      # Check if value actually changed
      if (old_value != is_checked) {
        # Update our deep copy with the new value
        current_active[[cluster_id]] <- is_checked
        
        # Update stable active status
        state$active_clusters <- current_active
        
        # Update "Select All" checkbox state
        all_active <- all(unlist(current_active))
        updateCheckboxInput(session, "selectAllClusters", value = all_active)
        
        # Signal update
        state$last_update <- Sys.time()
      }
      
      # Clear updating flag after a delay
      shinyjs::delay(100, {
        updating_from_checkbox(FALSE)
      })
    })
    
    # Update labels when the update button is clicked
    observeEvent(input$updateAllLabels, {
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      req(state$all_clusters)
      
      # Only proceed if we have temporary labels to apply
      if (!is.null(state$temp_labels) && length(state$temp_labels) > 0) {
        # Start with current stable labels
        current_labels <- stable_labels()
        if (is.null(current_labels)) {
          current_labels <- list()
        }
        
        # Merge in temporary labels
        for (cluster_key in names(state$temp_labels)) {
          current_labels[[cluster_key]] <- state$temp_labels[[cluster_key]]
        }
        
        # Apply the update
        stable_labels(current_labels)
        
        # Also update state for consistency
        state$cluster_labels <- current_labels
        
        # Clear temp labels to avoid duplicate application
        state$temp_labels <- list()
        
        # Trigger update in other modules
        state$last_update <- Sys.time()
        
        # Show notification
        showNotification("Cluster labels updated", type = "message")
      }
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
            
            # Ensure temp_labels is initialized properly
            state$temp_labels <- list()
            
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
            
            # Update checkboxes using JavaScript for better consistency
            if (!is.null(saved_state$active_clusters)) {
              active_clusters <- names(saved_state$active_clusters[unlist(saved_state$active_clusters) == TRUE])
              session$sendCustomMessage(
                type = "updateClusterCheckboxes",
                message = list(clusters = active_clusters)
              )
              
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
    return(list(
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
        # Signal the update
        state$last_update <- Sys.time()
      },
      updateActiveStatus = function(new_active) {
        state$active_clusters <- new_active
        
        # Update UI to reflect the new active status
        active_clusters <- names(new_active[unlist(new_active) == TRUE])
        session$sendCustomMessage(
          type = "updateClusterCheckboxes",
          message = list(clusters = active_clusters)
        )
      },
      last_update = reactive({ state$last_update }),
      getFullState = getFullState,
      setFullState = setFullState,
      updateFromButton = function() {
        # Only proceed if we have temporary labels to apply
        if (!is.null(state$temp_labels) && length(state$temp_labels) > 0) {
          # Start with current stable labels
          current_labels <- stable_labels()
          if (is.null(current_labels)) {
            current_labels <- list()
          }
          
          # Merge in temporary labels
          for (cluster_key in names(state$temp_labels)) {
            current_labels[[cluster_key]] <- state$temp_labels[[cluster_key]]
          }
          
          # Apply the update
          stable_labels(current_labels)
          
          # Also update state for consistency
          state$cluster_labels <- current_labels
          
          # Clear temp labels to avoid duplicate application
          state$temp_labels <- list()
          
          # Trigger update in other modules
          state$last_update <- Sys.time()
          
          # Show notification
          showNotification("Cluster labels updated", type = "message")
        } else {
          showNotification("No label changes to save", type = "message")
        }
      }
    ))
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