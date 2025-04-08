#' @title Quality Control Module UI
#' @description Creates the UI for the quality control module which allows users to 
#'   visualize and filter scRNA-seq data based on quality metrics.
#' @param id The module ID
#' @return A Shiny UI element containing the QC interface
#' @export
qcUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Container for plots
    div(class = "qc-plots",
        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
            h4(style = "margin: 0;", "Quality Control metrics"),
            downloadButton(ns("downloadQCPlot"), "Save Plot", 
                           class = "btn-sm btn-success")
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
      sample_labels = NULL  # Add this to track sample labels
    )
    
    # Add this near the beginning of the function
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
    
    # Watch for changes in sample management active status
    observe({
      req(sample_management)
      
      # Get active samples
      active_samples <- sample_management$getActiveSampleIds()
      
      # Update filtered samples
      if (!identical(values$filtered_samples, active_samples)) {
        values$filtered_samples <- active_samples
        
        # Force re-rendering of plots
        values$plot_update_trigger <- runif(1)
      }
    })
    
    # Watch for changes in sample management labels
    observe({
      req(sample_management)
      
      # Get current sample labels
      current_labels <- sample_management$getSampleLabels()
      
      # Update if they've changed
      if (!identical(values$sample_labels, current_labels)) {
        values$sample_labels <- current_labels
        
        # Force re-rendering of plots
        values$plot_update_trigger <- runif(1)
      }
    })
    
    # Watch for changes in condition management active status
    observe({
      req(condition_management)
      
      # Get active conditions and the condition column
      active_conditions <- condition_management$getActiveConditions()
      condition_column <- condition_management$getConditionColumn()
      
      # Update filtered conditions and condition column
      if (!identical(values$filtered_conditions, active_conditions) || 
          !identical(values$condition_column, condition_column)) {
        
        values$filtered_conditions <- active_conditions
        values$condition_column <- condition_column
        
        # Force re-rendering of plots
        values$plot_update_trigger <- runif(1)
      }
    })
    
    # Handle sample selection for plot with filtering
    plot_data <- reactive({
      req(seurat_data())
      
      # Add this line to make the rendering reactive to the update trigger
      trigger <- values$plot_update_trigger
      
      # Create a filtered Seurat object
      filtered_seurat <- seurat_data()
      
      # Apply sample filtering if available
      if (!is.null(values$filtered_samples) && length(values$filtered_samples) > 0) {
        # Filter to show only cells from active samples
        cells_to_keep <- filtered_seurat$sample %in% values$filtered_samples
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        }
      }
      
      # Apply condition filtering if available
      if (!is.null(values$condition_column) && !is.null(values$filtered_conditions) && 
          length(values$filtered_conditions) > 0 && 
          values$condition_column %in% colnames(filtered_seurat@meta.data)) {
        
        # Filter to show only cells from active conditions
        cells_to_keep <- filtered_seurat@meta.data[[values$condition_column]] %in% values$filtered_conditions
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        }
      }
      
      # Apply sample limit to prevent overcrowding (get first 5 samples)
      filtered_seurat <- getSamplesToPlot(filtered_seurat)
      
      # If we have sample labels, apply them to the plot data
      if (!is.null(values$sample_labels)) {
        # For each sample in the dataset, update its display name if we have a label
        current_samples <- unique(filtered_seurat$sample)
        for (sample in current_samples) {
          if (sample %in% names(values$sample_labels)) {
            # Create a temporary variable to use for plotting
            filtered_seurat$sample_label <- filtered_seurat$sample
            # Replace sample IDs with their labels
            filtered_seurat$sample_label[filtered_seurat$sample == sample] <- values$sample_labels[[sample]]
          }
        }
      }
      
      # Return the filtered data
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
    
    # The rest of your code remains the same...
    
    # Download handler for the QC plot
    output$downloadQCPlot <- downloadHandler(
      filename = function() {
        paste("qc_metrics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        saveQCPlot(file, qc_plot())
      }
    )
    
    # Render filter controls
    output$filterControls <- renderUI({
      req(seurat_data())
      createFilterControls(session, values)
    })
    
    # Observers to update the stored filter parameters when inputs change
    observeEvent(input$minFeature, {
      values$min_feature <- input$minFeature
    })
    observeEvent(input$maxFeature, {
      values$max_feature <- input$maxFeature
    })
    observeEvent(input$maxMT, {
      values$max_mt <- input$maxMT
    })
    
    # Process data when the button is clicked
    observeEvent(input$processSeurat, {
      req(seurat_data(), input$minFeature, input$maxFeature, input$maxMT)
      
      # Store current parameter values
      values$min_feature <- input$minFeature
      values$max_feature <- input$maxFeature
      values$max_mt <- input$maxMT
      
      # Create a filtered Seurat object
      filtered_seurat <- seurat_data()
      num_cells_before <- ncol(filtered_seurat)
      
      # Prepare message components
      sample_info <- ""
      condition_info <- ""
      filter_summary <- paste0("QC thresholds - Min features: ", input$minFeature, 
                               ", Max features: ", input$maxFeature, 
                               ", Max MT%: ", input$maxMT)
      
      # Apply sample filtering if available
      if (!is.null(values$filtered_samples) && length(values$filtered_samples) > 0) {
        # Filter to show only cells from active samples
        cells_to_keep <- filtered_seurat$sample %in% values$filtered_samples
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          showNotification("No cells match the active sample selection. Cannot process data.", 
                           type = "error")
          return(NULL)
        }
      }
      
      # Apply condition filtering if available
      if (!is.null(values$condition_column) && !is.null(values$filtered_conditions) && 
          length(values$filtered_conditions) > 0 && 
          values$condition_column %in% colnames(filtered_seurat@meta.data)) {
        
        # Filter to show only cells from active conditions
        cells_to_keep <- filtered_seurat@meta.data[[values$condition_column]] %in% values$filtered_conditions
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          showNotification("No cells match the active condition selection. Cannot process data.", 
                           type = "error")
          return(NULL)
        }
      }
      
      # Set the status message
      values$qc_status_message <- paste0(
        "Processing QC with ", ncol(filtered_seurat), " cells", 
        sample_info, condition_info, "; ", filter_summary
      )
      
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
      values$qc_status_message <- paste0(
        "QC complete - Retained ", num_cells_after, " cells (", 
        round(num_cells_after/num_cells_before * 100, 1), 
        "%). Filtered out ", cells_filtered, " cells.", 
        sample_info, condition_info, "; ", filter_summary
      )
    })
    
    # Return the processed data
    processed_seurat <- reactive({
      values$filtered_data
    })
    
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
      
      # Process data with these parameters
      filtered_seurat <- seurat_data()
      
      # Apply sample filtering if available
      if (!is.null(values$filtered_samples) && length(values$filtered_samples) > 0) {
        cells_to_keep <- filtered_seurat$sample %in% values$filtered_samples
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        }
      }
      
      # Apply condition filtering if available
      if (!is.null(values$condition_column) && !is.null(values$filtered_conditions) && 
          length(values$filtered_conditions) > 0 && 
          values$condition_column %in% colnames(filtered_seurat@meta.data)) {
        
        cells_to_keep <- filtered_seurat@meta.data[[values$condition_column]] %in% values$filtered_conditions
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        }
      }
      
      # Process with the specified parameters
      values$filtered_data <- processQCFiltering(
        filtered_seurat, 
        min_feature, 
        max_feature, 
        max_mt
      )
      
      return(values$filtered_data)
    }
    
    setProcessedData <- function(new_data) {
      if (!is.null(new_data)) {
        values$filtered_data <- new_data
        # Update internal parameters based on the loaded data if possible
        tryCatch({
          # You could add code here to extract parameters from the loaded data
          # or just trust that the UI update has set the correct parameters
        }, error = function(e) {
          print(paste("Error updating QC parameters:", e$message))
        })
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

#' @title Get Samples to Plot
#' @description Selects a subset of samples to display in QC plots (limits to first 5 if more exist).
#' @param seurat_obj Seurat object containing the data
#' @return A Seurat object with only the selected samples
#' @keywords internal
getSamplesToPlot <- function(seurat_obj) {
  all_samples <- unique(seurat_obj$sample)
  if (length(all_samples) > 5) {
    samples_to_plot <- all_samples[1:5]
    cells_to_keep <- seurat_obj$sample %in% samples_to_plot
    return(subset(seurat_obj, cells = cells_to_keep))
  }
  return(seurat_obj)
}

#' @title Create QC Plot
#' @description Creates violin plots for key QC metrics (feature count, UMI count, mitochondrial percentage).
#' @param seurat_obj Seurat object to visualize
#' @return A ggplot object with violin plots for QC metrics
#' @keywords internal
createQCPlot <- function(seurat_obj) {
  # Check if we have sample_label column
  if ("sample_label" %in% colnames(seurat_obj@meta.data)) {
    # Use sample_label instead of sample for grouping
    VlnPlot(seurat_obj, 
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
            group.by = "sample_label",
            ncol = 3)
  } else {
    # Fall back to original behavior
    VlnPlot(seurat_obj, 
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
            group.by = "sample",
            ncol = 3)
  }
}

#' @title Save QC Plot
#' @description Saves the QC plot to a file with error handling.
#' @param file File path to save to
#' @param plot The plot to save
#' @return None
#' @keywords internal
saveQCPlot <- function(file, plot) {
  tryCatch({
    ggsave(file, plot = plot, device = "png", width = 12, height = 8, dpi = 300)
  }, error = function(e) {
    # Fallback to base R graphics if ggsave fails
    png(file, width = 12, height = 8, units = "in", res = 300)
    print(plot)
    dev.off()
  })
}

#' @title Create Filter Controls
#' @description Creates UI elements for filtering the Seurat object based on QC metrics.
#' @param session The current Shiny session
#' @return A UI element with filter controls
#' @keywords internal
createFilterControls <- function(session, values) {
  ns <- session$ns
  tagList(
    tags$div(
      id = ns("filter_controls"),
      tags$br(),
      tags$strong("Please adjust filtering parameters:"),
      tags$br(),
      tags$br(),
      numericInput(ns("minFeature"), "Minimum Features:", values$min_feature),
      numericInput(ns("maxFeature"), "Maximum Features:", values$max_feature),
      numericInput(ns("maxMT"), "Maximum MT %:", values$max_mt),
      actionButton(ns("processSeurat"), "Filter and run PCA")
    )
  )
}

#' @title Process QC Filtering
#' @description Filters cells based on QC metrics and runs preprocessing steps (normalization, scaling, PCA).
#' @param seurat_obj The Seurat object to process
#' @param min_feature Minimum feature count threshold
#' @param max_feature Maximum feature count threshold
#' @param max_mt Maximum mitochondrial percentage threshold
#' @return A processed Seurat object
#' @keywords internal
processQCFiltering <- function(seurat_obj, min_feature, max_feature, max_mt) {
  withProgress(message = 'Processing data', value=0, {
    # Filter cells
    incProgress(0.05, detail = "Filtering cells")
    seurat <- subset(seurat_obj, subset = nFeature_RNA > min_feature &
                       nFeature_RNA < max_feature &
                       percent.mt < max_mt)
    
    # Normalize data
    incProgress(0.2, detail = "Normalizing data")
    seurat <- JoinLayers(seurat)
    seurat <- NormalizeData(seurat)
    
    # Find variable features
    incProgress(0.2, detail = "Finding variable features")
    seurat <- FindVariableFeatures(seurat)
    
    # Scale data
    incProgress(0.2, detail = "Scaling data")
    seurat <- ScaleData(seurat)
    
    # Run PCA
    incProgress(0.2, detail = "Running PCA")
    seurat <- RunPCA(seurat, npcs = 50)
    
    return(seurat)
  })
}