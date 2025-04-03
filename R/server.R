# R/server.R

#' @title Build Application Server
#' @description Constructs the complete server logic for the single-cell RNA-seq analysis application,
#'   initializing all modules and establishing reactive data flow between components.
#' @return A Shiny server function that defines the reactive behavior of the application
#' @export
buildServer <- function() {
  function(input, output, session) {
    print("Starting server initialization")
    
    # Initialize metadata module first
    metadata_module <- metadataServer("metadata")
    
    metadata_handler <- reactive({
      metadata_module
    })
    
    # Initialize data input module
    data_input <- dataInputServer("dataInput", metadata_module = metadata_module)
    
    # Extract Seurat object with proper reactive chain
    seurat_data <- reactive({
      req(data_input())
      data_input()
    })
    
    # Initialize management modules
    sample_management <- sampleManagementServer("sampleManagement", seurat_data)
    condition_management <- conditionManagementServer("conditionManagement", seurat_data, metadata_module)
    
    # Chain the reactive values through the analysis modules
    qc_module <- qcServer("qc", seurat_data, sample_management, condition_management)
    # Extract the reactive data from the module
    processed_seurat <- reactive({
      if (is.list(qc_module) && is.function(qc_module$data)) {
        qc_module$data()
      } else {
        qc_module()
      }
    })
    
    cluster_management <- clusterManagementServer("clusterManagement", clustered_seurat = reactive({
      # Get the reactive data from dimred_module
      if (!is.null(dimred_module) && is.list(dimred_module) && is.function(dimred_module$data)) {
        return(dimred_module$data())
      } else if (!is.null(clustered_seurat())) {
        return(clustered_seurat())
      }
      return(NULL)
    }))
    
    # Initialize the dimension reduction module
    dimred_module <- dimensionReductionServer("dimRed", 
                                              processed_seurat, 
                                              sample_management, 
                                              condition_management, 
                                              cluster_management)
    
    # Extract the reactive data from the dimred module
    clustered_seurat <- reactive({
      if (is.list(dimred_module) && is.function(dimred_module$data)) {
        dimred_module$data()
      } else {
        dimred_module()
      }
    })
    
    de_module <- deAnalysisServer("de", clustered_seurat, cluster_management, sample_management, condition_management)
    
    # Create reactive values to track completion of each step
    steps_completed <- reactiveValues(
      data_input = FALSE,
      metadata = FALSE,
      qc = FALSE,
      dimred = FALSE,
      clustering = FALSE,
      de = FALSE,
      sample_management = FALSE,
      condition_management = FALSE
    )
    
    # Initialize save/load module
    loaded_analysis <- saveLoadServer("saveLoad", 
                                      seurat_data, 
                                      metadata_module, 
                                      qc_module, 
                                      dimred_module, 
                                      de_module,
                                      steps_completed,
                                      session)
    
    # Handle loaded analysis data with enhanced UI restoration
    observe({
      data <- loaded_analysis()
      if (is.null(data)) return()
      
      # Verify we have some of the key components
      if (is.null(data$seurat_data) && is.null(data$processed_seurat)) {
        showNotification("Warning: Loaded analysis may be incomplete", type = "warning")
        return()
      }
      
      # Update data objects first
      if (!is.null(data$seurat_data)) {
        data_input(data$seurat_data)
      }
      
      # Main UI restoration workflow with sequential delays
      # This ensures each UI element is updated in the correct order
      shinyjs::delay(500, {
        withProgress(message = 'Restoring analysis state...', value = 0, {
          # First restore QC parameters and simulate button click if needed
          incProgress(0.3, detail = "Restoring QC parameters")
          qc_processed <- FALSE
          
          if (!is.null(data$ui_state) && !is.null(data$ui_state$qc_params)) {
            qc_processed <- restoreQcUI(data$ui_state$qc_params, session)
          }
          
          # If QC wasn't processed in the saved state, we're done
          if (!qc_processed) {
            incProgress(1.0, detail = "Restoration complete")
            return()
          }
          
          # After QC is processed, restore PCA/dimension reduction state
          # Use a longer delay to ensure QC completes first
          shinyjs::delay(1000, {
            incProgress(0.6, detail = "Restoring dimension reduction")
            
            if (!is.null(data$ui_state) && !is.null(data$ui_state$pca_params)) {
              restorePcaUI(data$ui_state$pca_params, session)
            }
            
            incProgress(1.0, detail = "Restoration complete")
            showNotification("Analysis state fully restored", type = "message")
          })
        })
      })
    })
    
    # Render cluster controls UI
    output$clusterControls <- renderUI({
      has_clustered <- !is.null(clustered_seurat())
      has_clusters <- FALSE
      
      if (has_clustered) {
        tryCatch({
          has_clusters <- "seurat_clusters" %in% colnames(clustered_seurat()@meta.data)
        }, error = function(e) {
          print(paste("Error checking for clusters:", e$message))
        })
      }
      
      if (!has_clustered || !has_clusters) {
        return(div(
          class = "alert alert-info",
          style = "margin-top: 10px; margin-bottom: 10px;",
          "Cluster controls will appear here after the clustering step is complete."
        ))
      }
      
      div(
        style = "margin-top: 10px;",
        div(
          style = "display: flex; align-items: center; margin-bottom: 10px;",
          checkboxInput("selectAllClusters", "Select All Clusters", value = TRUE),
          tags$div(style = "margin-left: 8px;",
                   actionButton("updateAllLabels", "Save Labels", 
                                class = "btn-sm btn-primary")
          )
        ),
        clusterManagementUI("clusterManagement")
      )
    })
    
    # Render sample controls UI
    output$sampleControls <- renderUI({
      has_samples <- !is.null(seurat_data())
      has_available_samples <- FALSE
      
      if (has_samples) {
        tryCatch({
          has_available_samples <- "sample" %in% colnames(seurat_data()@meta.data) && 
            length(unique(seurat_data()$sample)) > 0
        }, error = function(e) {
          print(paste("Error checking for samples:", e$message))
        })
      }
      
      if (!has_samples || !has_available_samples) {
        return(div(
          class = "alert alert-info",
          style = "margin-top: 10px; margin-bottom: 10px;",
          "Sample controls will appear here after data is loaded."
        ))
      }
      
      div(
        style = "margin-top: 10px;",
        div(
          style = "display: flex; align-items: center; margin-bottom: 10px;",
          checkboxInput("selectAllSamples", "Select All Samples", value = TRUE),
          tags$div(style = "margin-left: 8px;",
                   actionButton("updateSampleLabels", "Save Labels", 
                                class = "btn-sm btn-primary")
          )
        ),
        sampleManagementUI("sampleManagement")
      )
    })
    
    # Render condition controls UI
    output$conditionControls <- renderUI({
      has_data <- !is.null(seurat_data())
      
      if (!has_data) {
        return(div(
          class = "alert alert-info",
          style = "margin-top: 10px; margin-bottom: 10px;",
          "Condition controls will appear here after data is loaded."
        ))
      }
      
      div(
        style = "margin-top: 10px;",
        conditionManagementUI("conditionManagement")
      )
    })
    
    # Handle the sample management all samples checkbox
    observeEvent(input$selectAllSamples, {
      if (is.function(sample_management$updateActiveStatus)) {
        samples_object <- seurat_data()
        if (!is.null(samples_object) && "sample" %in% colnames(samples_object@meta.data)) {
          all_samples <- unique(samples_object$sample)
          new_active <- setNames(
            rep(input$selectAllSamples, length(all_samples)),
            all_samples
          )
          sample_management$updateActiveStatus(new_active)
        }
      }
    })
    
    # Handle sample label updates button
    observeEvent(input$updateSampleLabels, {
      if (is.function(sample_management$getSampleLabels)) {
        # This will trigger the module's internal update function
        # The module handles the actual label updates
      }
    })
    
    # Handle the condition management all conditions checkbox
    observeEvent(input$selectAllConditions, {
      if (is.function(condition_management$updateActiveStatus)) {
        condition_column <- condition_management$getConditionColumn()
        if (!is.null(condition_column) && !is.null(seurat_data())) {
          # Get all available conditions for the selected column
          condition_values <- unique(as.character(seurat_data()@meta.data[[condition_column]]))
          if (length(condition_values) > 0) {
            new_active <- setNames(
              rep(input$selectAllConditions, length(condition_values)),
              condition_values
            )
            condition_management$updateActiveStatus(new_active)
          }
        }
      }
    }, ignoreInit = TRUE)
    
    # Handle condition label updates button
    observeEvent(input$updateConditionLabels, {
      if (is.function(condition_management$getConditionLabels)) {
        # This will trigger the module's internal update function
        # The module handles the actual label updates
      }
    })
    
    # Handle the cluster management all clusters checkbox
    observeEvent(input$selectAllClusters, {
      if (is.function(cluster_management$updateActiveStatus)) {
        clustered_obj <- clustered_seurat()
        if (!is.null(clustered_obj) && "seurat_clusters" %in% colnames(clustered_obj@meta.data)) {
          all_clusters <- as.character(sort(unique(clustered_obj$seurat_clusters)))
          new_active <- setNames(
            rep(input$selectAllClusters, length(all_clusters)),
            all_clusters
          )
          cluster_management$updateActiveStatus(new_active)
        }
      }
    })
    
    # Handle cluster label updates button
    observeEvent(input$updateAllLabels, {
      if (is.function(cluster_management$getClusterLabels)) {
        # This will trigger the module's internal update function
        # The module handles the actual label updates
      }
    })
    
    # Setup observers for tracking module completion status
    setupObservers(steps_completed, seurat_data, metadata_module, qc_module, 
                   dimred_module, de_module, sample_management, condition_management)
    
    # Setup dynamic section rendering
    setupSections(input, output, seurat_data, metadata_handler, processed_seurat, 
                  clustered_seurat, session)
    
    # Setup navigation panel
    setupNavigation(output, steps_completed)
    
    print("Server initialization complete")
  }
}