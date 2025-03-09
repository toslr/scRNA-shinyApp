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
    
    # Add a reactive value to track UMAP type (2D or 3D)
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
    
    # Compute UMAP - combined observer for both 2D and 3D
    seurat_with_umap <- reactiveVal(NULL)
    
    # 2D UMAP handler
    observeEvent(input$run2DUMAP, {
      req(processed_seurat(), input$nDims)
      withProgress(message = 'Computing 2D UMAP...', {
        seurat <- processed_seurat()
        # Run 2D UMAP
        seurat <- RunUMAP(seurat, dims = 1:input$nDims, n.components = 2)
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
      
      if (umap_type() == "2D") {
        div(
          div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
              h4(style = "margin: 0;", "2D UMAP Visualization"),
              downloadButton(ns("downloadUMAPPlot"), "Save Plot", 
                             class = "btn-sm btn-success")
          ),
          plotOutput(ns("umapPlot"), height = "800px")
        )
      } else {
        div(
          div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
              h4(style = "margin: 0;", "3D UMAP Visualization"),
              downloadButton(ns("download3DUMAPPlot"), "Save Screenshot", 
                             class = "btn-sm btn-success")
          ),
          plotlyOutput(ns("umap3DPlot"), height = "800px")
        )
      }
    })
    
    # 2D UMAP plot as reactive expression
    umap_plot <- reactive({
      req(seurat_with_umap())
      req("umap" %in% names(seurat_with_umap()@reductions))
      req(umap_type() == "2D")
      
      # Check whether to color by sample or just show a single color
      if ("sample" %in% colnames(seurat_with_umap()@meta.data)) {
        DimPlot(seurat_with_umap(), reduction = "umap", group.by = "sample")
      } else {
        DimPlot(seurat_with_umap(), reduction = "umap")
      }
    })
    
    # 3D UMAP plot as reactive expression
    umap3d_plot <- reactive({
      req(seurat_with_umap())
      req("umap" %in% names(seurat_with_umap()@reductions))
      req(umap_type() == "3D")
      
      # Extract UMAP coordinates
      umap_data <- Embeddings(seurat_with_umap()[["umap"]])
      
      # Extract sample or cluster information
      if ("seurat_clusters" %in% colnames(seurat_with_umap()@meta.data)) {
        # Color by clusters if available
        color_var <- as.factor(seurat_with_umap()$seurat_clusters)
        plot_title <- "3D UMAP - Colored by Clusters"
      } else if ("sample" %in% colnames(seurat_with_umap()@meta.data)) {
        # Color by sample ID if clusters aren't available
        color_var <- as.factor(seurat_with_umap()$sample)
        plot_title <- "3D UMAP - Colored by Sample"
      } else {
        # Default coloring if neither is available
        color_var <- rep("Sample", ncol(seurat_with_umap()))
        plot_title <- "3D UMAP"
      }
      
      # Create custom color palette with sufficient colors
      if (length(unique(color_var)) > 1) {
        colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(color_var)))
      } else {
        colors <- "blue" # Default color if only one group
      }
      
      # Create plotly 3D scatter plot with improved legend
      p <- plot_ly(x = umap_data[,1], 
                   y = umap_data[,2], 
                   z = umap_data[,3],
                   type = "scatter3d",
                   mode = "markers",
                   color = color_var,
                   colors = colors,
                   marker = list(size = 3, opacity = 0.7)) %>%
        layout(title = plot_title,
               scene = list(xaxis = list(title = "UMAP 1"),
                            yaxis = list(title = "UMAP 2"),
                            zaxis = list(title = "UMAP 3")))
      
      # Custom legend with larger markers
      p$x$data <- lapply(p$x$data, function(d) {
        if (!is.null(d$marker)) {
          # Make legend markers larger
          d$showlegend <- TRUE
          d$marker$sizeref <- 0.2  # Smaller value = larger legend markers
        }
        return(d)
      })
      
      # Adjust legend position and size
      p <- p %>% layout(
        margin = list(r = 120),  # Make room for legend
        legend = list(
          font = list(size = 12),
          itemsizing = "constant",
          y = 0.5
        )
      )
      
      return(p)
    })
    
    # Render UMAP plots based on type
    output$umapPlot <- renderPlot({
      req(umap_type() == "2D")
      umap_plot()
    })
    
    output$umap3DPlot <- renderPlotly({
      req(umap_type() == "3D")
      umap3d_plot()
    })
    
    # Download handler for 2D UMAP plot
    output$downloadUMAPPlot <- downloadHandler(
      filename = function() {
        paste("umap_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = umap_plot(), device = "png", width = 8, height = 8, dpi = 300)
      }
    )
    
    # Download handler for 3D UMAP plot
    output$download3DUMAPPlot <- downloadHandler(
      filename = function() {
        paste("umap3d_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        # Save the plotly object as a static PNG
        plotly::save_image(umap3d_plot(), file, width = 800, height = 800)
      }
    )
    
    output$clusterSection <- renderUI({
      req(clustered_seurat())
      req("umap" %in% names(clustered_seurat()@reductions))
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      
      if (umap_type() == "2D") {
        div(
          div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
              h4(style = "margin: 0;", "Cluster Visualization (2D)"),
              downloadButton(ns("downloadClusterPlot"), "Save Plot", 
                             class = "btn-sm btn-success")
          ),
          plotOutput(ns("clusterPlot"), height = "800px")
        )
      } else {
        div(
          div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
              h4(style = "margin: 0;", "Cluster Visualization (3D)"),
              downloadButton(ns("download3DClusterPlot"), "Save Screenshot", 
                             class = "btn-sm btn-success")
          ),
          plotlyOutput(ns("cluster3DPlot"), height = "800px")
        )
      }
    })
    
    # 2D Cluster plot as reactive expression
    cluster_plot <- reactive({
      req(clustered_seurat())
      req("umap" %in% names(clustered_seurat()@reductions))
      req(umap_type() == "2D")
      
      DimPlot(clustered_seurat(), reduction = "umap", label = TRUE)#, label.size = 8)
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
      clusters_factor <- as.factor(clusters)
      
      # Create custom color palette
      colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(levels(clusters_factor)))
      
      # Create plotly 3D scatter plot with clusters
      p <- plot_ly(x = umap_data[,1], 
                   y = umap_data[,2], 
                   z = umap_data[,3],
                   type = "scatter3d",
                   mode = "markers",
                   color = clusters_factor,
                   colors = colors,
                   marker = list(size = 3, opacity = 0.7)) %>%
        layout(title = "3D UMAP - Clusters",
               scene = list(xaxis = list(title = "UMAP 1"),
                            yaxis = list(title = "UMAP 2"),
                            zaxis = list(title = "UMAP 3")))
      
      # Custom legend with larger markers
      p$x$data <- lapply(p$x$data, function(d) {
        if (!is.null(d$marker)) {
          # Make legend markers larger
          d$showlegend <- TRUE
          d$marker$sizeref <- 0.2  # Smaller value = larger legend markers
        }
        return(d)
      })
      
      # Adjust legend position and size
      p <- p %>% layout(
        margin = list(r = 120),  # Make room for legend
        legend = list(
          font = list(size = 12),
          itemsizing = "constant",
          y = 0.5
        )
      )
      
      return(p)
    })
    
    # Render cluster plots based on type
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
        ggsave(file, plot = cluster_plot(), device = "png", width = 8, height = 8, dpi = 300)
      }
    )
    
    # Download handler for 3D cluster plot
    output$download3DClusterPlot <- downloadHandler(
      filename = function() {
        paste("cluster3d_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        # Save the plotly object as a static PNG
        plotly::save_image(cluster3d_plot(), file, width = 800, height = 800)
      }
    )
    
    return(clustered_seurat)
  })
}