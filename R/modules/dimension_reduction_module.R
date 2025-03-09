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

dimensionReductionServer <- function(id, processed_seurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # State management
    values <- reactiveValues(
      umap_type = "2D",
      seurat_with_umap = NULL,
      clustered_seurat = NULL
    )
    
    # Find elbow point for optimal dimensions
    suggested_dims <- reactive({
      req(processed_seurat())
      seurat <- processed_seurat()
      
      pca_data <- Embeddings(seurat, "pca")
      stdev <- Stdev(seurat[["pca"]])
      var_explained <- stdev^2 / sum(stdev^2)
      
      dims <- 1:length(var_explained)
      find_elbow(dims, var_explained)
    })
    
    # Create elbow plot
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
        saveElbowPlot(file, processed_seurat(), suggested_dims())
      }
    )
    
    output$suggestedDims <- renderText({
      req(suggested_dims())
      paste("Suggested number of dimensions (based on elbow point):", suggested_dims(),
            ". Please adjust the number of PC for reduction if needed.")
    })
    
    # Render dimension controls
    output$dimControls <- renderUI({
      req(processed_seurat(), suggested_dims())
      createDimControls(ns, suggested_dims(), ncol(Embeddings(processed_seurat(), "pca")))
    })
    
    # Handle 2D UMAP generation
    observeEvent(input$run2DUMAP, {
      req(processed_seurat(), input$nDims)
      withProgress(message = 'Computing 2D UMAP...', {
        seurat <- runUMAP(processed_seurat(), input$nDims, n_components = 2)
        values$seurat_with_umap <- seurat
        values$umap_type <- "2D"
      })
    })
    
    # Handle 3D UMAP generation
    observeEvent(input$run3DUMAP, {
      req(processed_seurat(), input$nDims)
      withProgress(message = 'Computing 3D UMAP...', {
        seurat <- runUMAP(processed_seurat(), input$nDims, n_components = 3)
        values$seurat_with_umap <- seurat
        values$umap_type <- "3D"
      })
    })
    
    # Render clustering controls
    output$clusterControls <- renderUI({
      req(values$seurat_with_umap)
      createClusteringControls(ns)
    })
    
    # Run clustering
    observeEvent(input$runClustering, {
      req(values$seurat_with_umap, input$nDims, input$resolution)
      withProgress(message = 'Clustering...', {
        seurat <- runClustering(values$seurat_with_umap, input$nDims, input$resolution)
        values$clustered_seurat <- seurat
      })
    })
    
    # Render UMAP section based on type
    output$umapSection <- renderUI({
      req(values$seurat_with_umap)
      req("umap" %in% names(values$seurat_with_umap@reductions))
      
      if (values$umap_type == "2D") {
        create2DUmapUI(ns)
      } else {
        create3DUmapUI(ns)
      }
    })
    
    # 2D UMAP plot
    umap_plot <- reactive({
      req(values$seurat_with_umap)
      req("umap" %in% names(values$seurat_with_umap@reductions))
      req(values$umap_type == "2D")
      
      create2DUmapPlot(values$seurat_with_umap)
    })
    
    # 3D UMAP plot
    umap3d_plot <- reactive({
      req(values$seurat_with_umap)
      req("umap" %in% names(values$seurat_with_umap@reductions))
      req(values$umap_type == "3D")
      
      create3DUmapPlot(values$seurat_with_umap)
    })
    
    # Render UMAP plots based on type
    output$umapPlot <- renderPlot({
      req(values$umap_type == "2D")
      umap_plot()
    })
    
    output$umap3DPlot <- renderPlotly({
      req(values$umap_type == "3D")
      umap3d_plot()
    })
    
    # Download handlers for UMAP plots
    output$downloadUMAPPlot <- downloadHandler(
      filename = function() {
        paste("umap_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = umap_plot(), device = "png", width = 8, height = 8, dpi = 300)
      }
    )
    
    # Cluster visualization section
    output$clusterSection <- renderUI({
      req(values$clustered_seurat)
      req("umap" %in% names(values$clustered_seurat@reductions))
      req("seurat_clusters" %in% colnames(values$clustered_seurat@meta.data))
      
      if (values$umap_type == "2D") {
        createCluster2DUI(ns)
      } else {
        createCluster3DUI(ns)
      }
    })
    
    # 2D Cluster plot
    cluster_plot <- reactive({
      req(values$clustered_seurat)
      req("umap" %in% names(values$clustered_seurat@reductions))
      req(values$umap_type == "2D")
      
      createCluster2DPlot(values$clustered_seurat)
    })
    
    # 3D Cluster plot
    cluster3d_plot <- reactive({
      req(values$clustered_seurat)
      req("umap" %in% names(values$clustered_seurat@reductions))
      req(values$umap_type == "3D")
      req("seurat_clusters" %in% colnames(values$clustered_seurat@meta.data))
      
      umap_dims <- dim(values$clustered_seurat@reductions$umap@cell.embeddings)[2]
      req(umap_dims >= 3, "3D UMAP not available")
      
      createCluster3DPlot(values$clustered_seurat)
    })
    
    # Render cluster plots
    output$clusterPlot <- renderPlot({
      req(values$umap_type == "2D")
      cluster_plot()
    })
    
    output$cluster3DPlot <- renderPlotly({
      req(values$umap_type == "3D")
      cluster3d_plot()
    })
    
    # Download handlers for cluster plots
    output$downloadClusterPlot <- downloadHandler(
      filename = function() {
        paste("cluster_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = cluster_plot(), device = "png", width = 8, height = 8, dpi = 300)
      }
    )
    
    return(reactive({ values$clustered_seurat }))
  })
}

# Helper Functions

find_elbow <- function(x, y) {
  # Focus on the rate of change and look for where the rate of change stabilizes
  diffs <- diff(y) / diff(x)
  
  window_size <- 3 # Adjustable
  rolling_std <- sapply(1:(length(diffs) - window_size), function(i) {
    sd(diffs[i:(i + window_size)])
  })
  
  threshold <- mean(rolling_std) * 0.1
  elbow_idx <- which(rolling_std < threshold)[1]
  
  return(x[elbow_idx])
}

saveElbowPlot <- function(file, seurat_obj, suggested_dims) {
  png(file, width = 3000, height = 1800, res = 300)
  
  print(ElbowPlot(seurat_obj, 
                  ndims = ncol(Embeddings(seurat_obj, "pca"))) +
          geom_vline(xintercept = suggested_dims, 
                     color = "red", 
                     linetype = "dashed"))
  
  dev.off()
}

createDimControls <- function(ns, suggested_dims, max_dims) {
  tagList(
    tags$div(
      id = ns("dimension_controls"),
      numericInput(ns("nDims"), 
                   "Number of dimensions for UMAP:", 
                   value = suggested_dims,
                   min = 2,
                   max = max_dims),
      div(style = "display: flex; gap: 10px; margin-top: 10px;",
          actionButton(ns("run2DUMAP"), "Run 2D UMAP"),
          actionButton(ns("run3DUMAP"), "Run Interactive 3D UMAP")
      )
    )
  )
}

runUMAP <- function(seurat_obj, n_dims, n_components = 2) {
  RunUMAP(seurat_obj, dims = 1:n_dims, n.components = n_components)
}

createClusteringControls <- function(ns) {
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
}

runClustering <- function(seurat_obj, n_dims, resolution) {
  seurat <- FindNeighbors(seurat_obj, dims = 1:n_dims)
  seurat <- FindClusters(seurat, resolution = resolution)
  return(seurat)
}

create2DUmapUI <- function(ns) {
  div(
    div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
        h4(style = "margin: 0;", "2D UMAP Visualization"),
        downloadButton(ns("downloadUMAPPlot"), "Save Plot", 
                       class = "btn-sm btn-success")
    ),
    plotOutput(ns("umapPlot"), height = "800px")
  )
}

create3DUmapUI <- function(ns) {
  div(
    h4(style = "margin: 0;", "3D UMAP Visualization"),
    plotlyOutput(ns("umap3DPlot"), height = "800px")
  )
}

create2DUmapPlot <- function(seurat_obj) {
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    DimPlot(seurat_obj, reduction = "umap", group.by = "sample")
  } else {
    DimPlot(seurat_obj, reduction = "umap")
  }
}

create3DUmapPlot <- function(seurat_obj) {
  # Extract UMAP coordinates
  umap_data <- Embeddings(seurat_obj[["umap"]])
  
  # Extract sample or cluster information
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    color_var <- as.factor(seurat_obj$seurat_clusters)
    plot_title <- "3D UMAP - Colored by Clusters"
  } else if ("sample" %in% colnames(seurat_obj@meta.data)) {
    color_var <- as.factor(seurat_obj$sample)
    plot_title <- "3D UMAP - Colored by Sample"
  } else {
    color_var <- rep("Sample", ncol(seurat_obj))
    plot_title <- "3D UMAP"
  }
  
  # Create custom color palette
  if (length(unique(color_var)) > 1) {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(color_var)))
  } else {
    colors <- "blue"
  }
  
  # Create 3D plot
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
      d$showlegend <- TRUE
      d$marker$sizeref <- 0.2
    }
    return(d)
  })
  
  # Adjust legend position and size
  p <- p %>% layout(
    margin = list(r = 120),
    legend = list(
      font = list(size = 12),
      itemsizing = "constant",
      y = 0.5
    )
  )
  
  return(p)
}

createCluster2DUI <- function(ns) {
  div(
    div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
        h4(style = "margin: 0;", "Cluster Visualization (2D)"),
        downloadButton(ns("downloadClusterPlot"), "Save Plot", 
                       class = "btn-sm btn-success")
    ),
    plotOutput(ns("clusterPlot"), height = "800px")
  )
}

createCluster3DUI <- function(ns) {
  div(
    h4(style = "margin: 0;", "Cluster Visualization (3D)"),
    plotlyOutput(ns("cluster3DPlot"), height = "800px")
  )
}

createCluster2DPlot <- function(seurat_obj) {
  DimPlot(seurat_obj, reduction = "umap", label = TRUE)
}

createCluster3DPlot <- function(seurat_obj) {
  # Extract UMAP coordinates
  umap_data <- Embeddings(seurat_obj[["umap"]])
  if (ncol(umap_data) < 3) {
    return(plotly_empty() %>% 
             layout(title = "Error: 3D UMAP data not available"))
  }
  
  # Extract cluster information
  clusters <- seurat_obj$seurat_clusters
  clusters_factor <- as.factor(clusters)
  
  # Create custom color palette
  colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(levels(clusters_factor)))
  
  # Create 3D plot
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
      d$showlegend <- TRUE
      d$marker$sizeref <- 0.2
    }
    return(d)
  })
  
  # Adjust legend position and size
  p <- p %>% layout(
    margin = list(r = 120),
    legend = list(
      font = list(size = 12),
      itemsizing = "constant",
      y = 0.5
    )
  )
  
  return(p)
}