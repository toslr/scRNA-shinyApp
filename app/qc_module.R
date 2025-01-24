qcUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("qcPlot"), height = "600px"),
    uiOutput(ns("filterControls"))
  )
}

qcServer <- function(id, seurat_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$qcPlot <- renderPlot({
      req(seurat_data())
      VlnPlot(seurat_data(), 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              group.by = "sample",
              ncol = 3)
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