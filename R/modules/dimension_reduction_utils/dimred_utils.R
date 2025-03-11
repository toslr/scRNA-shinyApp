# R/modules/dimension_reduction_utils/dimred_utils.R

# ==================== Helper Functions for Dimension Reduction ====================

# Find elbow point for optimal dimensions
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

# Save elbow plot to file
saveElbowPlot <- function(file, seurat_obj, suggested_dims) {
  png(file, width = 3000, height = 1800, res = 300)
  
  print(ElbowPlot(seurat_obj, 
                  ndims = ncol(Embeddings(seurat_obj, "pca"))) +
          geom_vline(xintercept = suggested_dims, 
                     color = "red", 
                     linetype = "dashed"))
  
  dev.off()
}

# Create dimension controls UI
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

# Run UMAP with specified dimensions and components
runUMAP <- function(seurat_obj, n_dims, n_components = 2) {
  RunUMAP(seurat_obj, dims = 1:n_dims, n.components = n_components)
}

# Create clustering controls UI
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

# Run clustering with specified dimensions and resolution
runClustering <- function(seurat_obj, n_dims, resolution) {
  seurat <- FindNeighbors(seurat_obj, dims = 1:n_dims)
  seurat <- FindClusters(seurat, resolution = resolution)
  return(seurat)
}

# ==================== 2D UMAP Visualization Functions ====================

# Create 2D UMAP UI
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

# Create 2D UMAP plot
create2DUmapPlot <- function(seurat_obj) {
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    DimPlot(seurat_obj, reduction = "umap", group.by = "sample")
  } else {
    DimPlot(seurat_obj, reduction = "umap")
  }
}

# Create 2D cluster UI
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

# Create 2D cluster plot
createCluster2DPlot <- function(seurat_obj) {
  DimPlot(seurat_obj, reduction = "umap", label = TRUE)
}

# ==================== 3D UMAP Visualization Functions ====================

# Create 3D UMAP UI
create3DUmapUI <- function(ns) {
  div(
    h4(style = "margin: 0;", "3D UMAP Visualization"),
    plotlyOutput(ns("umap3DPlot"), height = "800px")
  )
}

# Create 3D UMAP plot
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

# Create 3D cluster UI
createCluster3DUI <- function(ns) {
  div(
    h4(style = "margin: 0;", "Cluster Visualization (3D)"),
    plotlyOutput(ns("cluster3DPlot"), height = "800px")
  )
}

# Create 3D cluster plot
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

# ==================== Gene UMAP Visualization Functions ====================

# Create Gene UMAP UI
createGeneUmapUI <- function(ns) {
  tagList(
    div(
      id = ns("gene_umap_section"),
      h4("Gene Expression UMAP"),
      p("Visualize expression of a specific gene on the UMAP."),
      fluidRow(
        column(6,
               div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
                   textInput(ns("geneQuery"), "Search for gene:", placeholder = "e.g. Sox10, Mbp, or ENSMUSG..."),
                   actionButton(ns("searchGene"), "Search", class = "btn-primary")
               )
        ),
        column(6,
               uiOutput(ns("geneSelectionUI"))
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] !== undefined && input['%s'] !== ''", ns("selectedGene"), ns("selectedGene")),
        div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
            h5(style = "margin: 0;", textOutput(ns("selectedGeneLabel")))
        ),
        uiOutput(ns("geneUmapVisSection"))
      )
    )
  )
}

# Create 2D Gene UMAP UI
create2DGeneUmapUI <- function(ns) {
  div(
    style = "margin-top: 15px;",
    div(style = "text-align: right; margin-bottom: 10px;",
        downloadButton(ns("downloadGeneUmap2DPlot"), "Save Plot", class = "btn-sm btn-success")
    ),
    plotOutput(ns("geneUmapPlot2D"), height = "800px")
  )
}

# Create 3D Gene UMAP UI
create3DGeneUmapUI <- function(ns) {
  div(
    style = "margin-top: 15px;",
    div(style = "text-align: right; margin-bottom: 10px;",
        downloadButton(ns("downloadGeneUmap3DPlot"), "Save Plot", class = "btn-sm btn-success")
    ),
    plotlyOutput(ns("geneUmapPlot3D"), height = "800px")
  )
}

# Create 2D Gene UMAP plot
createGeneUmapPlot <- function(seurat_obj, gene_id, type = "2D") {
  # Check if gene exists in the dataset
  if (!(gene_id %in% rownames(seurat_obj))) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("Gene", gene_id, "not found in dataset")) + 
             theme_void())
  }
  
  # Get gene symbol if available
  gene_symbol <- if (!is.null(seurat_obj@misc$gene_mapping) && 
                     gene_id %in% names(seurat_obj@misc$gene_mapping) &&
                     !is.na(seurat_obj@misc$gene_mapping[gene_id])) {
    seurat_obj@misc$gene_mapping[gene_id]
  } else {
    gene_id
  }
  
  # Create feature plot without cluster labels
  p <- FeaturePlot(seurat_obj, 
                   features = gene_id, 
                   reduction = "umap",
                   order = TRUE,  # Plot high-expression cells on top
                   label = FALSE,  # No cluster labels
                   cols = c("lightgrey", "red"),
                   pt.size = 1.2) 
  
  # Add custom title and theme elements
  p + theme_minimal() +
    labs(title = paste("Expression of", gene_symbol, "in UMAP")) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}

# Create 3D Gene UMAP plot
createGeneUmap3DPlot <- function(seurat_obj, gene_id) {
  # Check if gene exists in the dataset
  if (!(gene_id %in% rownames(seurat_obj))) {
    return(plotly_empty() %>% 
             layout(title = paste("Gene", gene_id, "not found in dataset")))
  }
  
  # Get gene symbol if available
  gene_symbol <- if (!is.null(seurat_obj@misc$gene_mapping) && 
                     gene_id %in% names(seurat_obj@misc$gene_mapping) &&
                     !is.na(seurat_obj@misc$gene_mapping[gene_id])) {
    seurat_obj@misc$gene_mapping[gene_id]
  } else {
    gene_id
  }
  
  # Get expression values
  expr_values <- GetAssayData(seurat_obj, slot = "data")[gene_id, ]
  
  # Extract UMAP coordinates
  umap_data <- Embeddings(seurat_obj[["umap"]])
  if (ncol(umap_data) < 3) {
    return(plotly_empty() %>% 
             layout(title = "Error: 3D UMAP data not available"))
  }
  
  # Extract cluster information for hover text only
  clusters <- seurat_obj$seurat_clusters
  
  # Determine color scale for expression values
  # Normalize expression to a 0-1 scale for better color gradients
  expr_min <- min(expr_values)
  expr_max <- max(expr_values)
  
  # Create normalized expression values for coloring (avoid division by zero)
  if (expr_min == expr_max) {
    normalized_expr <- rep(0.5, length(expr_values))
  } else {
    normalized_expr <- (expr_values - expr_min) / (expr_max - expr_min)
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(
    umap1 = umap_data[, 1],
    umap2 = umap_data[, 2],
    umap3 = umap_data[, 3],
    expr = expr_values,
    norm_expr = normalized_expr,
    cluster = as.factor(clusters)
  )
  
  # Create 3D plot without cluster labels
  p <- plot_ly(
    data = plot_data,
    x = ~umap1, 
    y = ~umap2, 
    z = ~umap3,
    color = ~expr,
    colors = colorRamp(c("lightgrey", "red")),
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 4,
      opacity = 0.7,
      colorbar = list(
        title = "Expression"
      )
    ),
    hoverinfo = "text",
    text = ~paste("Cluster:", cluster, "<br>Expression:", round(expr, 3))
  ) %>%
    layout(
      title = paste("Expression of", gene_symbol, "in UMAP"),
      scene = list(
        xaxis = list(title = "UMAP 1"),
        yaxis = list(title = "UMAP 2"),
        zaxis = list(title = "UMAP 3")
      )
    )
  
  return(p)
}