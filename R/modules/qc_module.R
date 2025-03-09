# R/modules/qc_module.R

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

# Helper Functions

getSamplesToPlot <- function(seurat_obj) {
  all_samples <- unique(seurat_obj$sample)
  if (length(all_samples) > 5) {
    samples_to_plot <- all_samples[1:5]
    cells_to_keep <- seurat_obj$sample %in% samples_to_plot
    return(subset(seurat_obj, cells = cells_to_keep))
  }
  return(seurat_obj)
}

createQCPlot <- function(seurat_obj) {
  VlnPlot(seurat_obj, 
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          group.by = "sample",
          ncol = 3)
}

saveQCPlot <- function(file, plot) {
  tryCatch({
    ggsave(file, plot = plot, device = "png", width = 12, height = 8, dpi = 300)
  }, error = function(e) {
    print(paste("Error using ggsave:", e$message))
    print("Trying alternative saving method...")
    
    png(file, width = 12, height = 8, units = "in", res = 300)
    print(plot)
    dev.off()
  })
}

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