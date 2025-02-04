# R/modules/dimension_reduction_module.R

dimensionReductionUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("elbowPlot"), height = "400px"),
    textOutput(ns("suggestedDims")),  # Added suggestion text
    uiOutput(ns("dimControls")),
    plotOutput(ns("umapPlot"), height = "400px"),
    uiOutput(ns("clusterControls")),
    plotOutput(ns("clusterPlot"), height = "400px")
  )
}

find_elbow <- function(x, y) {
  # Focus on the rate of change
  diffs <- diff(y) / diff(x)
  
  # Look for where the rate of change stabilizes
  # Use a smaller window to be more sensitive to early changes
  window_size <- 3
  rolling_std <- sapply(1:(length(diffs) - window_size), function(i) {
    sd(diffs[i:(i + window_size)])
  })
  
  # Find where the variability in slope significantly decreases
  threshold <- mean(rolling_std) * 0.1  # Adjust this threshold as needed
  elbow_idx <- which(rolling_std < threshold)[1]
  
  return(x[elbow_idx])
}

dimensionReductionServer <- function(id, processed_seurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Calculate optimal dimensions
    suggested_dims <- reactive({
      req(processed_seurat())
      seurat <- processed_seurat()
      
      # Extract eigenvalues from PCA
      pca_data <- Embeddings(seurat, "pca")
      stdev <- Stdev(seurat[["pca"]])
      var_explained <- stdev^2 / sum(stdev^2)
      
      # Find elbow point
      dims <- 1:length(var_explained)
      optimal_dims <- find_elbow(dims, var_explained)
      
      return(optimal_dims)
    })
    
    output$elbowPlot <- renderPlot({
      req(processed_seurat())
      ElbowPlot(processed_seurat(), 
                ndims = ncol(Embeddings(processed_seurat(), "pca"))) +
        geom_vline(xintercept = suggested_dims(), 
                   color = "red", 
                   linetype = "dashed")
    })
    
    output$suggestedDims <- renderText({
      req(suggested_dims())
      paste("Suggested number of dimensions (based on elbow point):", suggested_dims())
    })
    
    output$dimControls <- renderUI({
      req(processed_seurat(), suggested_dims())
      tagList(
        tags$div(
          id = ns("dimension_controls"),
          renderPrint("Please adjust the number of PC for reduction"),
          numericInput(ns("nDims"), 
                       "Number of dimensions for UMAP:", 
                       value = suggested_dims(),  # Use suggested dims as default
                       min = 2,
                       max = ncol(Embeddings(processed_seurat(), "pca"))),
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
          numericInput(ns("resolution"), 
                       "Clustering Resolution:", 
                       0.5, 
                       min = 0, 
                       max = 2, 
                       step = 0.01),
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