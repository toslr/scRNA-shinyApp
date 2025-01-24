dimensionReductionUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("elbowPlot"), height = "400px"),
    uiOutput(ns("dimControls")),
    plotOutput(ns("umapPlot"), height = "400px"),
    uiOutput(ns("clusterControls")),
    plotOutput(ns("clusterPlot"), height = "400px")
  )
}

dimensionReductionServer <- function(id, processed_seurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$elbowPlot <- renderPlot({
      req(processed_seurat())
      ElbowPlot(processed_seurat(), 
                ndims = ncol(Embeddings(processed_seurat(), "pca")))
    })
    
    output$dimControls <- renderUI({
      req(processed_seurat())
      tagList(
        tags$div(
          id = ns("dimension_controls"),
          renderPrint("Please adjust the number of PC for reduction"),
          numericInput(ns("nDims"), "Number of dimensions for UMAP:", 10),
          actionButton(ns("runUMAP"), "Run UMAP")
        )
      )
    })
    
    seurat_with_umap <- eventReactive(input$runUMAP, {
      req(processed_seurat(), input$nDims)
      withProgress(message = 'Computing UMAP...', {
        seurat <- processed_seurat()
        seurat <- RunUMAP(seurat, dims = 1:input$nDims)
        seurat
      })
    })
    
    output$clusterControls <- renderUI({
      req(seurat_with_umap())
      tagList(
        tags$div(
          id = ns("clustering_controls"),
          renderPrint("Please adjust clustering resolution:"),
          numericInput(ns("resolution"), "Clustering Resolution:", 0.5, min = 0, max = 2, step = 0.01),
          actionButton(ns("runClustering"), "Run Clustering")
        )
      )
    })
    
    clustered_seurat <- eventReactive(input$runClustering, {
      req(seurat_with_umap(), input$nDims, input$resolution)
      withProgress(message = 'Clustering...', {
        seurat <- seurat_with_umap()
        seurat <- FindNeighbors(seurat, dims = 1:input$nDims)
        seurat <- FindClusters(seurat, resolution = input$resolution)
        seurat
      })
    })
    
    output$umapPlot <- renderPlot({
      req(seurat_with_umap())
      req("umap" %in% names(seurat_with_umap()@reductions))
      DimPlot(seurat_with_umap(), reduction = "umap")
    })
    
    output$clusterPlot <- renderPlot({
      req(clustered_seurat())
      req("umap" %in% names(clustered_seurat()@reductions))
      DimPlot(clustered_seurat(), reduction = "umap", label = TRUE)
    })
    
    return(clustered_seurat)
  })
}