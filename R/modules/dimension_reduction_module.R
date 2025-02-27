# R/modules/dimension_reduction_module.R

dimensionReductionUI <- function(id) {
  ns <- NS(id)
  tagList(
    div(
      uiOutput(ns("elbowSaveButton")),
      plotOutput(ns("elbowPlot"), height = "400px"),
      textOutput(ns("suggestedDims")),
    ),
    uiOutput(ns("dimControls")),
    uiOutput(ns("umapSection")),
    uiOutput(ns("clusterControls")),
    uiOutput(ns("clusterSection"))
  )
}

find_elbow <- function(x, y) {
  # Focus on the rate of change and look for where the rate of change stabilizes
  diffs <- diff(y) / diff(x)
  
  window_size <- 3 # Adjustable
  rolling_std <- sapply(1:(length(diffs) - window_size), function(i) {
    sd(diffs[i:(i + window_size)])
  })
  
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
    
    # Create elbow plot as reactive expression
    elbow_plot <- reactive({
      req(processed_seurat(), suggested_dims())
      ElbowPlot(processed_seurat(), 
                ndims = ncol(Embeddings(processed_seurat(), "pca"))) +
        geom_vline(xintercept = suggested_dims(), 
                   color = "red", 
                   linetype = "dashed")
    })
    
    # Render elbow plot
    output$elbowPlot <- renderPlot({
      elbow_plot()
    })
    
    output$elbowSaveButton <- renderUI({
      req(processed_seurat(), suggested_dims())
      div(style = "margin-top: 10px; text-align: right;",
          downloadButton(ns("downloadElbowPlot"), "Save Plot", 
                         class = "btn-sm btn-success"))
    })
    
    # Download handler for elbow plot
    output$downloadElbowPlot <- downloadHandler(
      filename = function() {
        paste("pca_elbow_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        # Create high-resolution PNG device directly
        png(file, width = 3000, height = 1800, res = 300)
        
        # Generate the plot directly to the device
        print(ElbowPlot(processed_seurat(), 
                        ndims = ncol(Embeddings(processed_seurat(), "pca"))) +
                geom_vline(xintercept = suggested_dims(), 
                           color = "red", 
                           linetype = "dashed"))
        
        # Close the device to save the file
        dev.off()
      }
    )
    
    output$suggestedDims <- renderText({
      req(suggested_dims())
      paste("Suggested number of dimensions (based on elbow point):", suggested_dims(),
            ". Please adjust the number of PC for reduction if needed.")
    })
    
    output$dimControls <- renderUI({
      req(processed_seurat(), suggested_dims())
      tagList(
        tags$div(
          id = ns("dimension_controls"),
          numericInput(ns("nDims"), 
                       "Number of dimensions for UMAP:", 
                       value = suggested_dims(),
                       min = 2,
                       max = ncol(Embeddings(processed_seurat(), "pca"))),
          actionButton(ns("runUMAP"), "Run UMAP")
        )
      )
    })
    
    # Compute UMAP
    seurat_with_umap <- eventReactive(input$runUMAP, {
      req(processed_seurat(), input$nDims)
      withProgress(message = 'Computing UMAP...', {
        seurat <- processed_seurat()
        seurat <- RunUMAP(seurat, dims = 1:input$nDims)
        seurat
      })
    })
    
    # Clustering xontrols
    output$clusterControls <- renderUI({
      req(seurat_with_umap())
      tagList(
        tags$div(
          id = ns("clustering_controls"),
          numericInput(ns("resolution"), 
                       "Please adjust clustering resolution:", 
                       0.5, 
                       min = 0, 
                       max = 2, 
                       step = 0.01),
          actionButton(ns("runClustering"), "Run Clustering")
        )
      )
    })
    
    # Cluster the data
    clustered_seurat <- eventReactive(input$runClustering, {
      req(seurat_with_umap(), input$nDims, input$resolution)
      withProgress(message = 'Clustering...', {
        seurat <- seurat_with_umap()
        seurat <- FindNeighbors(seurat, dims = 1:input$nDims)
        seurat <- FindClusters(seurat, resolution = input$resolution)
        seurat
      })
    })
    
    output$umapSection <- renderUI({
      req(seurat_with_umap())
      req("umap" %in% names(seurat_with_umap()@reductions))
      
      div(
        div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
            h4(style = "margin: 0;", "UMAP Visualization"),
            downloadButton(ns("downloadUMAPPlot"), "Save Plot", 
                           class = "btn-sm btn-success")
        ),
        plotOutput(ns("umapPlot"), height = "400px")
      )
    })
    
    # UMAP plot as reactive expression
    umap_plot <- reactive({
      req(seurat_with_umap())
      req("umap" %in% names(seurat_with_umap()@reductions))
      DimPlot(seurat_with_umap(), reduction = "umap")
    })
    
    # Render UMAP plot
    output$umapPlot <- renderPlot({
      umap_plot()
    })
    
    # Download handler for UMAP plot
    output$downloadUMAPPlot <- downloadHandler(
      filename = function() {
        paste("umap_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = umap_plot(), device = "png", width = 8, height = 6, dpi = 300)
      }
    )
    
    output$clusterSection <- renderUI({
      req(clustered_seurat())
      req("umap" %in% names(clustered_seurat()@reductions))
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      
      div(
        div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
            h4(style = "margin: 0;", "Cluster Visualization"),
            downloadButton(ns("downloadClusterPlot"), "Save Plot", 
                           class = "btn-sm btn-success")
        ),
        plotOutput(ns("clusterPlot"), height = "400px")
      )
    })
    
    # Cluster plot as reactive expression
    cluster_plot <- reactive({
      req(clustered_seurat())
      req("umap" %in% names(clustered_seurat()@reductions))
      DimPlot(clustered_seurat(), reduction = "umap", label = TRUE)
    })
    
    # Render cluster plot
    output$clusterPlot <- renderPlot({
      cluster_plot()
    })
    
    # Download handler for cluster plot
    output$downloadClusterPlot <- downloadHandler(
      filename = function() {
        paste("cluster_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = cluster_plot(), device = "png", width = 8, height = 6, dpi = 300)
      }
    )
    
    return(clustered_seurat)
  })
}