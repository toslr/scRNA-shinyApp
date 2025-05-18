#' @title Quality Control Module UI
#' @description Creates the UI for the quality control module which allows users to 
#'   visualize and filter scRNA-seq data based on quality metrics.
#' @param id The module ID
#' @return A Shiny UI element containing the QC interface
#' @export
qcUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Species-specific QC suggestions
    uiOutput(ns("qcSuggestions")),
    
    # Container for plots
    div(class = "qc-plots",
        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
            h4(style = "margin: 0;", "Quality Control metrics"),
            downloadButton(ns("downloadQCPlot"), "Save Plot", 
                           class = "btn-sm btn-success"),
            downloadButton(ns("downloadQCData"), "Save Data", 
                           class = "btn-sm btn-info")
        ),
        plotOutput(ns("qcPlot"), height = "600px")
    ),
    
    # Container for controls
    div(class = "qc-controls",
        uiOutput(ns("filterControls"))
    ),
    
    # Add QC status message
    uiOutput(ns("qcStatusMessage"))
  )
}

#' @title Quality Control Module Server
#' @description Server logic for the QC module that processes data filtering and visualization.
#' @param id The module ID
#' @param seurat_data Reactive expression containing the Seurat object
#' @param sample_management Optional sample management module for filtering by samples
#' @param condition_management Optional condition management module for filtering by conditions
#' @return A reactive expression containing the processed and filtered Seurat object
#' @export
qcServer <- function(id, seurat_data, sample_management = NULL, condition_management = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values to store state
    values <- reactiveValues(
      filtered_data = NULL,
      filtered_samples = NULL,
      filtered_conditions = NULL,
      condition_column = NULL,
      plot_update_trigger = runif(1),
      min_feature = 500,
      max_feature = 5000,
      max_mt = 5,
      sample_labels = NULL
    )
    
    # Status message reactive value
    values$qc_status_message <- NULL
    
    # Render the QC status message
    output$qcStatusMessage <- renderUI({
      if (!is.null(values$qc_status_message)) {
        div(
          class = "alert alert-info",
          style = "margin-top: 10px; margin-bottom: 10px;",
          icon("info-circle"),
          values$qc_status_message
        )
      }
    })
    
    # Add flag to control UI updates
    inhibit_ui_updates <- reactiveVal(FALSE)
    
    # Add debouncing for UI updates to prevent multiple rapid redraws
    observe({
      if (inhibit_ui_updates()) {
        invalidateLater(300) # Wait 300ms before allowing updates again
        inhibit_ui_updates(FALSE)
      }
    })
    
    # Watch for changes in sample management active status
    observe({
      req(sample_management)
      
      active_samples <- sample_management$getActiveSampleIds()
      
      if (!identical(values$filtered_samples, active_samples)) {
        values$filtered_samples <- active_samples
        values$plot_update_trigger <- runif(1)
      }
    })
    
    # Watch for changes in sample management labels
    observe({
      req(sample_management)
      
      current_labels <- sample_management$getSampleLabels()
      
      if (!identical(values$sample_labels, current_labels)) {
        values$sample_labels <- current_labels
        values$plot_update_trigger <- runif(1)
      }
    })
    
    # Watch for changes in condition management active status
    observe({
      req(condition_management)
      
      active_conditions <- condition_management$getActiveConditions()
      condition_column <- condition_management$getConditionColumn()
      
      if (!identical(values$filtered_conditions, active_conditions) || 
          !identical(values$condition_column, condition_column)) {
        
        values$filtered_conditions <- active_conditions
        values$condition_column <- condition_column
        values$plot_update_trigger <- runif(1)
      }
    })
    
    # Handle sample selection for plot with filtering
    plot_data <- reactive({
      req(seurat_data())
      
      # Add this line to make the rendering reactive to the update trigger
      trigger <- values$plot_update_trigger
      
      # Apply filtering to the Seurat object
      filtered_seurat <- filter_seurat_for_qc(
        seurat_data(),
        values$filtered_samples,
        values$filtered_conditions,
        values$condition_column,
        values$sample_labels
      )
      
      return(filtered_seurat)
    })
    
    # Create violin plot as a reactive expression
    qc_plot <- reactive({
      req(plot_data())
      createQCPlot(plot_data())
    })
    
    # Render the plot
    output$qcPlot <- renderPlot({
      qc_plot()
    })
    
    output$qcSuggestions <- renderUI({
      req(seurat_data())
      
      # Get QC suggestions based on detected species
      species_info <- get_species_from_seurat(seurat_data())
      species <- species_info$species
      
      # Get species-specific suggestions
      suggestions <- get_qc_suggestions(species)
      
      # Update values reactiveValue
      values$min_feature <- suggestions$min_feature
      values$max_feature <- suggestions$max_feature
      values$max_mt <- suggestions$max_mt
      
      # Only show suggestions if we detected a species
      if (!is.null(species) && species != "auto") {
        div(
          class = "alert alert-info",
          icon("info-circle"),
          paste("QC suggestions for", species, "data:"),
          tags$ul(
            tags$li(paste("Min Features:", suggestions$min_feature)),
            tags$li(paste("Max Features:", suggestions$max_feature)),
            tags$li(paste("Max MT%:", suggestions$max_mt))
          ),
          "These are general guidelines - adjust based on your specific dataset."
        )
      }
    })
    
    # Download handler for the QC plot
    output$downloadQCPlot <- downloadHandler(
      filename = function() {
        paste("qc_metrics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        saveQCPlot(file, qc_plot())
      }
    )
    
    # Download handler for the QC data
    output$downloadQCData <- downloadHandler(
      filename = function() {
        paste("qc_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
      },
      content = function(file) {
        # Extract the QC data
        qc_data <- extractQCData(plot_data())
        
        # Write to CSV
        write.csv(qc_data, file, row.names = FALSE)
      }
    )
    
    # Render filter controls
    output$filterControls <- renderUI({
      req(seurat_data())
      # Use isolate to prevent UI rebuilding while typing
      isolate({
        createFilterControls(session, values)
      })
    })
    
    # Observers to update the stored filter parameters when inputs change
    # Use isolate to prevent unwanted reactivity
    observeEvent(input$minFeature, {
      isolate({
        values$min_feature <- input$minFeature
      })
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    
    observeEvent(input$maxFeature, {
      isolate({
        values$max_feature <- input$maxFeature
      })
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    
    observeEvent(input$maxMT, {
      isolate({
        values$max_mt <- input$maxMT
      })
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    
    # Process data when the button is clicked
    observeEvent(input$processSeurat, {
      req(seurat_data(), input$minFeature, input$maxFeature, input$maxMT)
      
      # Inhibit UI updates during processing
      inhibit_ui_updates(TRUE)
      
      # Store current parameter values
      values$min_feature <- input$minFeature
      values$max_feature <- input$maxFeature
      values$max_mt <- input$maxMT
      
      # Create a filtered Seurat object
      filtered_seurat <- seurat_data()
      num_cells_before <- ncol(filtered_seurat)
      
      # Apply sample and condition filtering
      filtered_seurat <- filter_seurat_for_processing(
        filtered_seurat,
        values$filtered_samples,
        values$filtered_conditions,
        values$condition_column
      )
      
      # Update status with filtering context
      status_components <- build_qc_status_message(
        filtered_seurat,
        values$filtered_samples,
        values$filtered_conditions,
        values$condition_column,
        input$minFeature,
        input$maxFeature,
        input$maxMT
      )
      
      values$qc_status_message <- status_components$initial_message
      
      # Now process the filtered Seurat object
      values$filtered_data <- processQCFiltering(
        filtered_seurat, 
        input$minFeature, 
        input$maxFeature, 
        input$maxMT
      )
      
      # Update status message with results
      num_cells_after <- ncol(values$filtered_data)
      cells_filtered <- num_cells_before - num_cells_after
      percent_retained <- round(num_cells_after/num_cells_before * 100, 1)
      
      values$qc_status_message <- paste0(
        "QC complete - Retained ", num_cells_after, "/", num_cells_before, " cells (",
        percent_retained, "%). Filtered out ", cells_filtered, " cells.", 
        status_components$sample_info, status_components$condition_info, "; ", status_components$filter_summary
      )
      
      # Release UI update inhibitor
      shinyjs::delay(200, {
        inhibit_ui_updates(FALSE)
      })
    })
    
    # Return the processed data
    processed_seurat <- reactive({
      values$filtered_data
    })
    
    # Process with specific parameters (for programmatic use)
    processWithParameters <- function(min_feature, max_feature, max_mt) {
      req(seurat_data())
      
      # Update UI inputs to match the provided parameters
      updateNumericInput(session, "minFeature", value = min_feature)
      updateNumericInput(session, "maxFeature", value = max_feature)
      updateNumericInput(session, "maxMT", value = max_mt)
      
      # Update internal values
      values$min_feature <- min_feature
      values$max_feature <- max_feature
      values$max_mt <- max_mt
      
      # Apply filtering and process
      filtered_seurat <- filter_seurat_for_processing(
        seurat_data(),
        values$filtered_samples,
        values$filtered_conditions,
        values$condition_column
      )
      
      # Process with the specified parameters
      values$filtered_data <- processQCFiltering(
        filtered_seurat, 
        min_feature, 
        max_feature, 
        max_mt
      )
      
      return(values$filtered_data)
    }
    
    # Set processed data (for loading saved state)
    setProcessedData <- function(new_data) {
      if (!is.null(new_data)) {
        values$filtered_data <- new_data
      }
    }
    
    # Return a list with the processed data and methods for state restoration
    return(list(
      data = processed_seurat,
      getFilterParams = function() {
        list(
          min_feature = values$min_feature,
          max_feature = values$max_feature,
          max_mt = values$max_mt
        )
      },
      processWithParameters = processWithParameters,
      setProcessedData = setProcessedData
    ))
  })
}