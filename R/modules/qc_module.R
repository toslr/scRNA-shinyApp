# R/modules/qc_module.R

qcUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Container for plots always exists
    div(class = "qc-plots",
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
    
    # Function to get samples to plot
    get_plot_data <- function(seurat_obj) {
      all_samples <- unique(seurat_obj$sample)
      if (length(all_samples) > 5) {
        samples_to_plot <- all_samples[1:5]
        cells_to_keep <- seurat_obj$sample %in% samples_to_plot
        return(subset(seurat_obj, cells = cells_to_keep))
      }
      return(seurat_obj)
    }
    
    # Separate reactive for plot data
    plot_data <- reactive({
      req(seurat_data())
      get_plot_data(seurat_data())
    })
    
    # Violin plots
    output$qcPlot <- renderPlot({
      req(plot_data())
      VlnPlot(plot_data(), 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              group.by = "sample",
              ncol = 3)
    })
    
    # Filter controls
    output$filterControls <- renderUI({
      req(seurat_data())
      tagList(
        tags$div(
          id = session$ns("filter_controls"),
          tags$br(),
          tags$strong("Please adjust filtering parameters:"),
          tags$br(),
          tags$br(),
          numericInput(session$ns("minFeature"), "Minimum Features:", 500),
          numericInput(session$ns("maxFeature"), "Maximum Features:", 5000),
          numericInput(session$ns("maxMT"), "Maximum MT %:", 5),
          actionButton(session$ns("processSeurat"), "Filter and run PCA")
        )
      )
    })
    
    # Run PCA
    processed_seurat <- eventReactive(input$processSeurat, {
      req(seurat_data(), input$minFeature, input$maxFeature, input$maxMT)
      withProgress(message = 'Processing data', value=0, {
        seurat <- seurat_data()
        incProgress(0.05, detail = "Filtering cells")
        seurat <- subset(seurat, subset = nFeature_RNA > input$minFeature &
                           nFeature_RNA < input$maxFeature &
                           percent.mt < input$maxMT)
        
        incProgress(0.2, detail = "Normalizing data")
        seurat <- JoinLayers(seurat)
        seurat <- NormalizeData(seurat)
        
        incProgress(0.2, detail = "Finding variable features")
        seurat <- FindVariableFeatures(seurat)
        
        incProgress(0.2, detail = "Scaling data")
        seurat <- ScaleData(seurat)
        
        incProgress(0.2, detail = "Running PCA")
        seurat <- RunPCA(seurat, npcs = 50)
        seurat
      })
    })
    
    return(processed_seurat)
  })
}