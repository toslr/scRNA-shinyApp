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
    
    umap_type <- reactiveVal("2D")
    
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
          div(style = "display: flex; gap: 10px; margin-top: 10px;",
              actionButton(ns("run2DUMAP"), "Run 2D UMAP"),
              actionButton(ns("run3DUMAP"), "Run Interactive 3D UMAP")
          )
        )
      )
    })
    
    # Compute UMAP
    seurat_with_umap <- reactiveVal(NULL)
    
    observeEvent(input$run2DUMAP, {
      req(processed_seurat(), input$nDims)
      withProgress(message = 'Computing 2D UMAP...', {
        seurat <- processed_seurat()
        seurat <- RunUMAP(seurat, dims = 1:input$nDims, n.components=2)
        seurat_with_umap(seurat)
        umap_type("2D")
      })
    })
    
    # 3D UMAP handler
    observeEvent(input$run3DUMAP, {
      req(processed_seurat(), input$nDims)
      withProgress(message = 'Computing 3D UMAP...', {
        seurat <- processed_seurat()
        # Run 3D UMAP
        seurat <- RunUMAP(seurat, dims = 1:input$nDims, n.components = 3)
        seurat_with_umap(seurat)
        umap_type("3D")
      })
    })
    
    # Clustering controls
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
      
      if (umap_type()=="2D") {
        div(
          div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
              h4(style = "margin: 0;", "2D UMAP Visualization"),
              downloadButton(ns("downloadUMAPPlot"), "Save Plot", 
                             class = "btn-sm btn-success")
          ),
          plotOutput(ns("umapPlot"), height = "400px")
        )
      } else {
        div(
          plotlyOutput(ns("umap3DPlot"), height = "600px")
        )
      }
    })
    
    # UMAP plot as reactive expression
    umap_plot <- reactive({
      req(seurat_with_umap())
      req("umap" %in% names(seurat_with_umap()@reductions))
      req(umap_type() == "2D")
      DimPlot(seurat_with_umap(), reduction = "umap")#, label = TRUE) +
        #theme(text = element_text(face = "bold")) +
        #guides(color = guide_legend(override.aes = list(size = 5))) +
        #geom_text(aes(label = clustered_seurat()$seurat_clusters),
        #          size = 5, fontface = "bold", check_overlap = TRUE)
    })
    
    # 3D UMAP plot as reactive expression
    umap3d_plot <- reactive({
      req(seurat_with_umap())
      req("umap" %in% names(seurat_with_umap()@reductions))
      req(umap_type() == "3D")
      
      # Extract UMAP coordinates
      umap_data <- Embeddings(seurat_with_umap()[["umap"]])
      
      # Extract cluster information if available
      if ("seurat_clusters" %in% colnames(seurat_with_umap()@meta.data)) {
        clusters <- seurat_with_umap()$seurat_clusters
        color_var <- as.factor(clusters)
        plot_title <- "3D UMAP - Colored by Clusters"
      } else {
        # Default coloring if no clusters yet
        color_var <- "blue"
        plot_title <- "3D UMAP"
      }
      
      colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(color_var)))
      
      # Create plotly 3D scatter plot
      plot_ly(x = umap_data[,1], 
              y = umap_data[,2], 
              z = umap_data[,3],
              type = "scatter3d",
              mode = "markers",
              color = color_var,
              colors=colors,
              marker = list(size = 3, opacity = 0.7)) %>%
        layout(legend = list(
          itemsizing = "constant",  # Forces legend symbols to be the same size
          itemwidth = 50,           # Width of the legend items
        )) %>%
        layout(title = "3D UMAP - Clusters",
               scene = list(xaxis = list(title = "UMAP 1"),
                            yaxis = list(title = "UMAP 2"),
                            zaxis = list(title = "UMAP 3")))
    })
    
    # Render UMAP plot
    output$umapPlot <- renderPlot({
      req(umap_type() == "2D")
      umap_plot()
    })
    
    output$umap3DPlot <- renderPlotly({
      req(umap_type() == "3D")
      umap3d_plot()
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
    
    # Download handler for 3D UMAP plot
    output$download3DUMAPPlot <- downloadHandler(
      filename = function() {
        paste("umap3d_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        # Save the plotly object as a static PNG
        plotly_obj <- umap3d_plot()
        plotly::save_image(plotly_obj, file, width = 8, height = 6, scale = 3)
      }
    )
    
    output$clusterSection <- renderUI({
      req(clustered_seurat())
      req("umap" %in% names(clustered_seurat()@reductions))
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      
      if (umap_type() == "2D") {
        div(
          div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
              h4(style = "margin: 0;", "Cluster Visualization"),
              downloadButton(ns("downloadClusterPlot"), "Save Plot", 
                             class = "btn-sm btn-success")
          ),
          plotOutput(ns("clusterPlot"), height = "400px")
        )
      } else {
        div(
          div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
              h4(style = "margin: 0;", "Cluster Visualization (3D)"),
              downloadButton(ns("download3DClusterPlot"), "Save Plot", 
                             class = "btn-sm btn-success")
          ),
          plotlyOutput(ns("cluster3DPlot"), height = "600px")
        )
      }
    })
    
    # 2D Cluster plot as reactive expression
    cluster_plot <- reactive({
      req(clustered_seurat())
      req("umap" %in% names(clustered_seurat()@reductions))
      req(umap_type() == "2D")
      DimPlot(clustered_seurat(), reduction = "umap", label = TRUE)
    })
    
    # 3D Cluster plot as reactive expression
    cluster3d_plot <- reactive({
      req(clustered_seurat())
      req("umap" %in% names(clustered_seurat()@reductions))
      req(umap_type() == "3D")
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      
      # Extract UMAP coordinates
      umap_data <- Embeddings(clustered_seurat()[["umap"]])
      
      # Extract cluster information
      clusters <- clustered_seurat()$seurat_clusters
      
      # Create plotly 3D scatter plot with clusters
      plot_ly(x = umap_data[,1], 
              y = umap_data[,2], 
              z = umap_data[,3],
              type = "scatter3d",
              mode = "markers",
              color = as.factor(clusters),
              marker = list(size = 3, opacity = 0.7)) %>%
        layout(title = "3D UMAP - Clusters",
               scene = list(xaxis = list(title = "UMAP 1"),
                            yaxis = list(title = "UMAP 2"),
                            zaxis = list(title = "UMAP 3")))
    })
    
    # Render cluster plot
    output$clusterPlot <- renderPlot({
      req(umap_type() == "2D")
      cluster_plot()
    })
    
    output$cluster3DPlot <- renderPlotly({
      req(umap_type() == "3D")
      cluster3d_plot()
    })
    
    # Download handler for 2D cluster plot
    output$downloadClusterPlot <- downloadHandler(
      filename = function() {
        paste("cluster_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = cluster_plot(), device = "png", width = 8, height = 6, dpi = 300)
      }
    )
    
    # Download handler for 3D cluster plot
    output$download3DClusterPlot <- downloadHandler(
      filename = function() {
        paste("cluster3d_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        # Save the plotly object as a static PNG
        plotly_obj <- cluster3d_plot()
        plotly::save_image(plotly_obj, file, width = 8, height = 6, scale = 3)
      }
    )
    
    return(clustered_seurat)
  })
}