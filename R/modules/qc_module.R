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
    )
  )
}

#' @title Quality Control Module Server
#' @description Server logic for the QC module that processes data filtering and visualization.
#' @param id The module ID
#' @param seurat_data Reactive expression containing the Seurat object
#' @return A reactive expression containing the processed and filtered Seurat object
#' @export
qcServer <- function(id, seurat_data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values to store state
    values <- reactiveValues(
      filtered_data = NULL
    )
    
    # Handle sample selection for plot
    plot_data <- reactive({
      req(seurat_data())
      getSamplesToPlot(seurat_data())
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
      createFilterControls(session)
    })
    
    # Process data when the button is clicked
    observeEvent(input$processSeurat, {
      req(seurat_data(), input$minFeature, input$maxFeature, input$maxMT)
      
      # Process the Seurat object
      values$filtered_data <- processQCFiltering(
        seurat_data(), 
        input$minFeature, 
        input$maxFeature, 
        input$maxMT
      )
    })
    
    # Return the processed data
    processed_seurat <- reactive({
      values$filtered_data
    })
    
    return(processed_seurat)
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
  VlnPlot(seurat_obj, 
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          group.by = "sample",
          ncol = 3)
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
createFilterControls <- function(session) {
  ns <- session$ns
  tagList(
    tags$div(
      id = ns("filter_controls"),
      tags$br(),
      tags$strong("Please adjust filtering parameters:"),
      tags$br(),
      tags$br(),
      numericInput(ns("minFeature"), "Minimum Features:", 500),
      numericInput(ns("maxFeature"), "Maximum Features:", 5000),
      numericInput(ns("maxMT"), "Maximum MT %:", 5),
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