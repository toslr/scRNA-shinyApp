# R/modules/cluster_management_module.R

clusterManagementUI <- function(id) {
  ns <- NS(id)
  
  div(class = "sidebar-section collapsible-section",
      div(class = "section-header",
          h4("Cluster Management", style = "display: inline;"),
          tags$button(
            class = "btn btn-link btn-collapse",
            tags$i(class = "fa fa-chevron-down")
          )
      ),
      div(class = "section-content",
          div(
            style = "display: flex; align-items: center; margin-bottom: 10px;",
            checkboxInput(ns("selectAllClusters"), "Select All Clusters", value = TRUE),
            tags$div(style = "margin-left: 8px;",
                     actionButton(ns("updateAllLabels"), "Save Labels", 
                                  class = "btn-sm btn-primary")
            )
          ),
          div(id = ns("clusterManagement"),
              style = "max-height: 300px; overflow-y: auto;",
              uiOutput(ns("clusterControls"))
          )
      )
  )
}

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

# Helper function: Create UI controls for each cluster
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

# Helper function: Get available clusters from Seurat object
getAvailableClusters <- function(seurat_obj) {
  if (is.null(seurat_obj) || !"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }
  sort(unique(seurat_obj$seurat_clusters))
}