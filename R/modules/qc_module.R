#R/modules/qc_module.R

qcUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("qcPlot"), height = "600px"),
    textOutput(ns("sampleInfo")),
    uiOutput(ns("filterControls"))
  )
}

qcServer <- function(id, seurat_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Function to get samples to plot
    get_plot_data <- function(seurat_obj) {
      all_samples <- unique(seurat_obj$sample)
      if (length(all_samples) > 5) {
        samples_to_plot <- all_samples[1:5]
        # Only subset if we have more than 5 samples
        cells_to_keep <- seurat_obj$sample %in% samples_to_plot
        return(subset(seurat_obj, cells = cells_to_keep))
      }
      # If 5 or fewer samples, return original object
      return(seurat_obj)
    }
    
    output$qcPlot <- renderPlot({
      req(seurat_data())
      seurat_obj <- seurat_data()
      
      plot_obj <- get_plot_data(seurat_data())
      VlnPlot(plot_obj, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              group.by = "sample",
              ncol = 3)
    })
    
    # Add information about sample selection
    output$sampleInfo <- renderText({
      req(seurat_data())
      total_samples <- length(unique(seurat_data()$sample))
      if (total_samples > 5) {
        return(paste("Showing QC plots for first 5 samples out of", total_samples, "total samples"))
      }
      return(paste("Showing QC plots for all", total_samples, "samples"))
    })
    
    output$filterControls <- renderUI({
      req(seurat_data())
      tagList(
        tags$div(
          id = ns("filter_controls"),
          renderPrint("Please adjust filtering parameters:"),
          numericInput(ns("minFeature"), "Minimum Features:", 500),
          numericInput(ns("maxFeature"), "Maximum Features:", 5000),
          numericInput(ns("maxMT"), "Maximum MT %:", 5),
          actionButton(ns("processSeurat"), "Filter and run PCA")
        )
      )
    })
    
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