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
    
    # Initialize state
    state <- list(
      cluster_labels = reactiveVal(NULL),
      active_clusters = reactiveVal(NULL),
      temp_labels = reactiveVal(NULL),
      label_inputs = reactiveValues()
    )
    
    # Get clusters from Seurat object
    observe({
      available_clusters <- getAvailableClusters(clustered_seurat())
      
      # Skip if no clusters available
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      # Initialize or update labels if needed
      current_labels <- state$cluster_labels()
      if (is.null(current_labels) || 
          !all(as.character(available_clusters) %in% names(current_labels))) {
        
        new_labels <- initializeClusterLabels(available_clusters, current_labels)
        
        state$cluster_labels(new_labels)
        state$temp_labels(new_labels)
      }
      
      # Initialize active status if needed
      current_active <- state$active_clusters()
      if (is.null(current_active) || 
          !all(as.character(available_clusters) %in% names(current_active))) {
        state$active_clusters(initializeActiveStatus(available_clusters, current_active))
      }
    })
    
    # Handle initial population when clusters change
    observeEvent(clustered_seurat(), {
      available_clusters <- getAvailableClusters(clustered_seurat())
      
      if (!is.null(available_clusters) && length(available_clusters) > 0) {
        # Initialize new labels and active status
        new_labels <- initializeClusterLabels(available_clusters)
        state$cluster_labels(new_labels)
        state$temp_labels(new_labels)
        
        new_active <- initializeActiveStatus(available_clusters)
        state$active_clusters(new_active)
      }
    }, ignoreInit = TRUE)
    
    # Setup UI for cluster controls
    output$clusterControls <- renderUI({
      available_clusters <- getAvailableClusters(clustered_seurat())
      
      # Show message if no clusters
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(div(
          class = "alert alert-info",
          "No clusters available. Please run the clustering step first."
        ))
      }
      
      # Get current labels and active status
      current_temp_labels <- state$temp_labels()
      current_active <- state$active_clusters()
      
      # Return UI
      createClusterControls(ns, available_clusters, current_temp_labels, current_active)
    })
    
    # Setup UI for clusters
    output$clusterRows <- renderUI({
      available_clusters <- getAvailableClusters(clustered_seurat())
      
      # Show message if no clusters
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(div(
          class = "alert alert-info",
          "No clusters available. Please run the clustering step first."
        ))
      }
      
      # Get current labels and active status
      current_temp_labels <- state$temp_labels()
      current_active <- state$active_clusters()
      
      # Return UI
      createClusterControls(ns, available_clusters, current_temp_labels, current_active)
    })
    
    # Handle Select All checkbox
    observeEvent(input$selectAllClusters, {
      available_clusters <- getAvailableClusters(clustered_seurat())
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      current_active <- state$active_clusters()
      for (cluster in available_clusters) {
        cluster_key <- as.character(cluster)
        current_active[cluster_key] <- input$selectAllClusters
      }
      state$active_clusters(current_active)
      for (cluster in available_clusters) {
        updateCheckboxInput(session, paste0("active_", cluster), value = input$selectAllClusters)
      }
    })
    
    # Handle active status updates for individual clusters
    observe({
      available_clusters <- getAvailableClusters(clustered_seurat())
      
      # Skip empty clusters
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      lapply(available_clusters, function(cluster) {
        input_id <- paste0("active_", cluster)
        
        # Only observe if this input exists
        if (!is.null(input[[input_id]])) {
          observeEvent(input[[input_id]], {
            current_active <- state$active_clusters()
            if (is.null(current_active)) {
              return(NULL)
            }
            
            cluster_key <- as.character(cluster)
            
            # Check if this cluster exists in our active_clusters
            if (cluster_key %in% names(current_active)) {
              current_active[cluster_key] <- input[[input_id]]
              state$active_clusters(current_active)
              
              # Update select all checkbox based on all cluster checkboxes
              all_selected <- all(unlist(current_active))
              updateCheckboxInput(session, "selectAllClusters", value = all_selected)
            }
          }, ignoreInit = TRUE)
        }
      })
    })
    
    # Update label_inputs when text changes
    observe({
      available_clusters <- getAvailableClusters(clustered_seurat())
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      for(cluster in available_clusters) {
        local({
          local_cluster <- cluster
          input_id <- paste0("label_", local_cluster)
          
          # Only observe if this input exists
          if (!is.null(input[[input_id]])) {
            observeEvent(input[[input_id]], {
              state$label_inputs[[input_id]] <- input[[input_id]]
            }, ignoreInit = TRUE)
          }
        })
      }
    })
    
    # Update temp_labels only when update button is clicked
    observeEvent(input$updateAllLabels, {
      available_clusters <- getAvailableClusters(clustered_seurat())
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      current_temp <- state$temp_labels()
      
      for(cluster in available_clusters) {
        input_id <- paste0("label_", cluster)
        if(!is.null(state$label_inputs[[input_id]])) {
          current_temp[as.character(cluster)] <- state$label_inputs[[input_id]]
        }
      }
      
      state$temp_labels(current_temp)
      state$cluster_labels(current_temp)
      
      showNotification("Cluster labels updated", type = "message")
    })
    
    # Return reactive expressions
    list(
      getClusterLabels = reactive({ state$cluster_labels() }),
      getActiveStatus = reactive({ state$active_clusters() }),
      getActiveClusterIds = reactive({
        current_active <- state$active_clusters()
        if (is.null(current_active)) return(NULL)
        as.numeric(names(current_active[current_active == TRUE]))
      }),
      getActiveClusterList = reactive({
        current_active <- state$active_clusters()
        if (is.null(current_active)) return(NULL)
        names(current_active[current_active == TRUE])
      }),
      updateLabels = function(new_labels) {
        state$cluster_labels(new_labels)
        state$temp_labels(new_labels)
      },
      updateActiveStatus = function(new_active) {
        state$active_clusters(new_active)
      }
    )
  })
}

#' @title Create Cluster Controls
#' @description Helper function that creates UI controls for each cluster,
#'   including active/inactive toggle and label editing.
#' @param ns Namespace function
#' @param available_clusters Vector of available cluster IDs
#' @param current_temp_labels Named vector of current cluster labels
#' @param current_active Named vector of current active status (TRUE/FALSE)
#' @return A UI element with controls for each cluster
#' @keywords internal
createClusterControls <- function(ns, available_clusters, current_temp_labels, current_active) {
  tagList(
    # Individual cluster controls
    lapply(available_clusters, function(cluster) {
      cluster_key <- as.character(cluster)
      is_active <- if (cluster_key %in% names(current_active)) {
        current_active[[cluster_key]]
      } else {
        TRUE  # Default to active if not found
      }
      
      current_label <- if (cluster_key %in% names(current_temp_labels)) {
        current_temp_labels[[cluster_key]]
      } else {
        paste("Cluster", cluster)
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