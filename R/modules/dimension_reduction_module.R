# R/modules/dimension_reduction_module.R

dimensionReductionUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Elbow plot section (unchanged)
    div(
      uiOutput(ns("elbowSaveButton")),
      plotOutput(ns("elbowPlot"), height = "400px"),
      textOutput(ns("suggestedDims")),
    ),
    
    # Dimension confirmation section
    div(
      style = "margin-top: 15px; margin-bottom: 20px;",
      fluidRow(
        column(6, 
               numericInput(ns("nDims"), 
                            "Number of dimensions for analysis:", 
                            value = 15,
                            min = 2,
                            max = 50)
        ),
        column(6,
               div(style = "margin-top: 25px;",
                   actionButton(ns("confirmDims"), "Confirm Dimensions", 
                                class = "btn-primary")
               )
        )
      )
    ),
    
    # Clustering controls (shown after dimensions are confirmed)
    uiOutput(ns("clusteringSection")),
    
    # UMAP visualization section (shown after clustering)
    uiOutput(ns("umapSection"))
  )
}

dimensionReductionServer <- function(id, processed_seurat, sample_management = NULL, condition_management = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # State management
    values <- reactiveValues(
      dims_confirmed = FALSE,
      clustering_done = FALSE,
      seurat_with_umap = NULL,
      clustered_seurat = NULL,
      left_umap_type = "2D",
      right_umap_type = "2D",
      left_color_by = "sample",
      right_color_by = "cluster"
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
    
    # Update suggested dimensions
    observe({
      req(suggested_dims())
      updateNumericInput(session, "nDims", value = suggested_dims())
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
    
    # Handle dimension confirmation
    observeEvent(input$confirmDims, {
      req(processed_seurat(), input$nDims)
      
      # Update state
      values$dims_confirmed <- TRUE
      
      # Pre-compute UMAPs for both 2D and 3D
      withProgress(message = 'Computing UMAPs...', {
        # Start with the original seurat object
        filtered_seurat <- processed_seurat()
        
        # Apply sample filtering if available
        if (!is.null(sample_management) && is.function(sample_management$getActiveSampleIds)) {
          active_samples <- sample_management$getActiveSampleIds()
          if (!is.null(active_samples) && length(active_samples) > 0) {
            filtered_seurat <- filterBySamples(filtered_seurat, active_samples)
          }
        }
        
        # Apply condition filtering if available
        if (!is.null(condition_management) && 
            is.function(condition_management$getConditionColumn) && 
            is.function(condition_management$getActiveConditions)) {
          
          condition_column <- condition_management$getConditionColumn()
          active_conditions <- condition_management$getActiveConditions()
          
          if (!is.null(condition_column) && !is.null(active_conditions) && length(active_conditions) > 0) {
            filtered_seurat <- filterByConditions(filtered_seurat, condition_column, active_conditions)
          }
        }
        
        # Run 2D UMAP
        incProgress(0.4, detail = "Computing 2D UMAP")
        seurat_2d <- runUMAP(filtered_seurat, input$nDims, n_components = 2, 
                             reduction.name = "umap2d", reduction.key = "UMAP2D_")
        
        # Run 3D UMAP
        incProgress(0.4, detail = "Computing 3D UMAP")
        seurat_3d <- runUMAP(seurat_2d, input$nDims, n_components = 3,
                             reduction.name = "umap3d", reduction.key = "UMAP3D_")
        
        # Store combined object
        values$seurat_with_umap <- seurat_3d
        
        incProgress(0.2, detail = "Done")
      })
    })
    
    # Clustering section UI
    output$clusteringSection <- renderUI({
      req(values$dims_confirmed)
      
      div(
        style = "margin-top: 20px; margin-bottom: 20px;",
        h3("Clustering"),
        fluidRow(
          column(6, 
                 numericInput(ns("resolution"), 
                              "Clustering resolution:", 
                              0.5, 
                              min = 0, 
                              max = 2, 
                              step = 0.1)
          ),
          column(6,
                 div(style = "margin-top: 25px;",
                     actionButton(ns("runClustering"), "Run Clustering", 
                                  class = "btn-primary")
                 )
          )
        ),
        # Show only after clustering is done
        if (values$clustering_done) {
          div(
            style = "margin-top: 15px;",
            downloadButton(ns("downloadClusterStats"), "Download Cluster Statistics", 
                           class = "btn-sm btn-success")
          )
        }
      )
    })
    
    # Run clustering
    observeEvent(input$runClustering, {
      req(values$seurat_with_umap, input$nDims, input$resolution)
      
      withProgress(message = 'Clustering...', {
        # Run clustering
        seurat <- FindNeighbors(values$seurat_with_umap, dims = 1:input$nDims)
        seurat <- FindClusters(seurat, resolution = input$resolution)
        
        # Store clustered object
        values$clustered_seurat <- seurat
        values$clustering_done <- TRUE
      })
    })
    
    # UMAP visualization section
    output$umapSection <- renderUI({
      req(values$clustering_done)
      
      div(
        style = "margin-top: 30px;",
        h3("UMAP Visualization"),
        
        # Side-by-side UMAP layout
        fluidRow(
          # Left UMAP
          column(6,
                 wellPanel(
                   div(
                     style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                     h4("UMAP 1"),
                     div(
                       style = "display: flex; gap: 10px;",
                       actionButton(ns("leftUMAP2D"), "2D", 
                                    class = if(values$left_umap_type == "2D") "btn-primary" else "btn-default"),
                       actionButton(ns("leftUMAP3D"), "3D", 
                                    class = if(values$left_umap_type == "3D") "btn-primary" else "btn-default")
                     )
                   ),
                   div(
                     style = "margin-bottom: 15px;",
                     selectInput(ns("leftColorBy"), "Color by:", 
                                 choices = getColoringOptions(values$clustered_seurat),
                                 selected = values$left_color_by)
                   ),
                   # Gene search for left UMAP (shown only when "gene" is selected)
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'gene'", ns("leftColorBy")),
                     div(
                       style = "margin-bottom: 15px;",
                       textInput(ns("leftGeneQuery"), "Search gene:", placeholder = "e.g. Sox10"),
                       actionButton(ns("searchLeftGene"), "Search", class = "btn-sm btn-primary"),
                       uiOutput(ns("leftGeneSelection"))
                     )
                   ),
                   # UMAP Visualization
                   div(
                     style = "text-align: right; margin-bottom: 5px;",
                     downloadButton(ns("downloadLeftUMAP"), "Save Plot", class = "btn-sm btn-success")
                   ),
                   # Plot output based on 2D/3D selection
                   if (values$left_umap_type == "2D") {
                     plotOutput(ns("leftUMAPPlot"), height = "500px")
                   } else {
                     plotlyOutput(ns("leftUMAPPlot3D"), height = "500px")
                   }
                 )
          ),
          
          # Right UMAP
          column(6,
                 wellPanel(
                   div(
                     style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                     h4("UMAP 2"),
                     div(
                       style = "display: flex; gap: 10px;",
                       actionButton(ns("rightUMAP2D"), "2D", 
                                    class = if(values$right_umap_type == "2D") "btn-primary" else "btn-default"),
                       actionButton(ns("rightUMAP3D"), "3D", 
                                    class = if(values$right_umap_type == "3D") "btn-primary" else "btn-default")
                     )
                   ),
                   div(
                     style = "margin-bottom: 15px;",
                     selectInput(ns("rightColorBy"), "Color by:", 
                                 choices = getColoringOptions(values$clustered_seurat),
                                 selected = values$right_color_by)
                   ),
                   # Gene search for right UMAP (shown only when "gene" is selected)
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'gene'", ns("rightColorBy")),
                     div(
                       style = "margin-bottom: 15px;",
                       textInput(ns("rightGeneQuery"), "Search gene:", placeholder = "e.g. Sox10"),
                       actionButton(ns("searchRightGene"), "Search", class = "btn-sm btn-primary"),
                       uiOutput(ns("rightGeneSelection"))
                     )
                   ),
                   # UMAP Visualization
                   div(
                     style = "text-align: right; margin-bottom: 5px;",
                     downloadButton(ns("downloadRightUMAP"), "Save Plot", class = "btn-sm btn-success")
                   ),
                   # Plot output based on 2D/3D selection
                   if (values$right_umap_type == "2D") {
                     plotOutput(ns("rightUMAPPlot"), height = "500px")
                   } else {
                     plotlyOutput(ns("rightUMAPPlot3D"), height = "500px")
                   }
                 )
          )
        )
      )
    })
    
    # Toggle UMAP types
    observeEvent(input$leftUMAP2D, {
      values$left_umap_type <- "2D"
    })
    
    observeEvent(input$leftUMAP3D, {
      values$left_umap_type <- "3D"
    })
    
    observeEvent(input$rightUMAP2D, {
      values$right_umap_type <- "2D"
    })
    
    observeEvent(input$rightUMAP3D, {
      values$right_umap_type <- "3D"
    })
    
    # Update color by selections
    observeEvent(input$leftColorBy, {
      values$left_color_by <- input$leftColorBy
    })
    
    observeEvent(input$rightColorBy, {
      values$right_color_by <- input$rightColorBy
    })
    
    # Gene search for left UMAP
    left_gene_search_results <- reactiveVal(NULL)
    left_selected_gene <- reactiveVal(NULL)
    
    observeEvent(input$searchLeftGene, {
      req(values$clustered_seurat, input$leftGeneQuery)
      left_gene_search_results(searchGenes(values$clustered_seurat, input$leftGeneQuery))
    })
    
    output$leftGeneSelection <- renderUI({
      results <- left_gene_search_results()
      if (is.null(results) || nrow(results) == 0) {
        div(p("No matching genes found"))
      } else {
        div(
          style = "margin-top: 10px;",
          selectInput(ns("leftSelectedGene"), 
                      paste("Select from", nrow(results), "genes:"),
                      choices = setNames(results$ensembl_id, results$gene_symbol),
                      selected = results$ensembl_id[1])
        )
      }
    })
    
    observeEvent(input$leftSelectedGene, {
      left_selected_gene(input$leftSelectedGene)
    })
    
    # Gene search for right UMAP
    right_gene_search_results <- reactiveVal(NULL)
    right_selected_gene <- reactiveVal(NULL)
    
    observeEvent(input$searchRightGene, {
      req(values$clustered_seurat, input$rightGeneQuery)
      right_gene_search_results(searchGenes(values$clustered_seurat, input$rightGeneQuery))
    })
    
    output$rightGeneSelection <- renderUI({
      results <- right_gene_search_results()
      if (is.null(results) || nrow(results) == 0) {
        div(p("No matching genes found"))
      } else {
        div(
          style = "margin-top: 10px;",
          selectInput(ns("rightSelectedGene"), 
                      paste("Select from", nrow(results), "genes:"),
                      choices = setNames(results$ensembl_id, results$gene_symbol),
                      selected = results$ensembl_id[1])
        )
      }
    })
    
    observeEvent(input$rightSelectedGene, {
      right_selected_gene(input$rightSelectedGene)
    })
    
    # Left UMAP plot functions
    output$leftUMAPPlot <- renderPlot({
      req(values$clustered_seurat, values$left_umap_type == "2D")
      
      # Generate the plot based on coloring option
      if (values$left_color_by == "sample") {
        DimPlot(values$clustered_seurat, reduction = "umap2d", group.by = "sample")
      } else if (values$left_color_by == "cluster") {
        DimPlot(values$clustered_seurat, reduction = "umap2d", label = TRUE)
      } else if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
        FeaturePlot(values$clustered_seurat, features = left_selected_gene(), reduction = "umap2d")
      } else {
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Please select a valid coloring option") + theme_void()
      }
    })
    
    output$leftUMAPPlot3D <- renderPlotly({
      req(values$clustered_seurat, values$left_umap_type == "3D")
      
      # Generate the plot based on coloring option
      if (values$left_color_by == "sample") {
        create3DUmapPlot(values$clustered_seurat, color_by = "sample", reduction = "umap3d")
      } else if (values$left_color_by == "cluster") {
        create3DUmapPlot(values$clustered_seurat, color_by = "cluster", reduction = "umap3d")
      } else if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
        create3DGeneUmapPlot(values$clustered_seurat, gene_id = left_selected_gene(), reduction = "umap3d")
      } else {
        plotly_empty() %>% layout(title = "Please select a valid coloring option")
      }
    })
    
    # Right UMAP plot functions
    output$rightUMAPPlot <- renderPlot({
      req(values$clustered_seurat, values$right_umap_type == "2D")
      
      # Generate the plot based on coloring option
      if (values$right_color_by == "sample") {
        DimPlot(values$clustered_seurat, reduction = "umap2d", group.by = "sample")
      } else if (values$right_color_by == "cluster") {
        DimPlot(values$clustered_seurat, reduction = "umap2d", label = TRUE)
      } else if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
        FeaturePlot(values$clustered_seurat, features = right_selected_gene(), reduction = "umap2d")
      } else {
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Please select a valid coloring option") + theme_void()
      }
    })
    
    output$rightUMAPPlot3D <- renderPlotly({
      req(values$clustered_seurat, values$right_umap_type == "3D")
      
      # Generate the plot based on coloring option
      if (values$right_color_by == "sample") {
        create3DUmapPlot(values$clustered_seurat, color_by = "sample", reduction = "umap3d")
      } else if (values$right_color_by == "cluster") {
        create3DUmapPlot(values$clustered_seurat, color_by = "cluster", reduction = "umap3d")
      } else if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
        create3DGeneUmapPlot(values$clustered_seurat, gene_id = right_selected_gene(), reduction = "umap3d")
      } else {
        plotly_empty() %>% layout(title = "Please select a valid coloring option")
      }
    })
    
    # Download handlers
    output$downloadLeftUMAP <- downloadHandler(
      filename = function() {
        suffix <- if (values$left_umap_type == "2D") "png" else "html"
        paste0("umap_left_", values$left_color_by, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", suffix)
      },
      content = function(file) {
        if (values$left_umap_type == "2D") {
          # Handle 2D plot
          p <- if (values$left_color_by == "sample") {
            DimPlot(values$clustered_seurat, reduction = "umap2d", group.by = "sample")
          } else if (values$left_color_by == "cluster") {
            DimPlot(values$clustered_seurat, reduction = "umap2d", label = TRUE)
          } else if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
            FeaturePlot(values$clustered_seurat, features = left_selected_gene(), reduction = "umap2d")
          }
          ggsave(file, plot = p, device = "png", width = 8, height = 8, dpi = 300)
        } else {
          # Handle 3D plot
          p <- if (values$left_color_by == "sample") {
            create3DUmapPlot(values$clustered_seurat, color_by = "sample", reduction = "umap3d")
          } else if (values$left_color_by == "cluster") {
            create3DUmapPlot(values$clustered_seurat, color_by = "cluster", reduction = "umap3d")
          } else if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
            create3DGeneUmapPlot(values$clustered_seurat, gene_id = left_selected_gene(), reduction = "umap3d")
          }
          saveWidget(p, file)
        }
      }
    )
    
    output$downloadRightUMAP <- downloadHandler(
      filename = function() {
        suffix <- if (values$right_umap_type == "2D") "png" else "html"
        paste0("umap_right_", values$right_color_by, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", suffix)
      },
      content = function(file) {
        if (values$right_umap_type == "2D") {
          # Handle 2D plot
          p <- if (values$right_color_by == "sample") {
            DimPlot(values$clustered_seurat, reduction = "umap2d", group.by = "sample")
          } else if (values$right_color_by == "cluster") {
            DimPlot(values$clustered_seurat, reduction = "umap2d", label = TRUE)
          } else if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
            FeaturePlot(values$clustered_seurat, features = right_selected_gene(), reduction = "umap2d")
          }
          ggsave(file, plot = p, device = "png", width = 8, height = 8, dpi = 300)
        } else {
          # Handle 3D plot
          p <- if (values$right_color_by == "sample") {
            create3DUmapPlot(values$clustered_seurat, color_by = "sample", reduction = "umap3d")
          } else if (values$right_color_by == "cluster") {
            create3DUmapPlot(values$clustered_seurat, color_by = "cluster", reduction = "umap3d")
          } else if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
            create3DGeneUmapPlot(values$clustered_seurat, gene_id = right_selected_gene(), reduction = "umap3d")
          }
          saveWidget(p, file)
        }
      }
    )
    
    # Download cluster statistics
    output$downloadClusterStats <- downloadHandler(
      filename = function() {
        paste0("cluster_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        req(values$clustered_seurat)
        
        # Get cluster statistics
        cluster_stats <- getClusterStats(values$clustered_seurat)
        write.csv(cluster_stats, file, row.names = FALSE)
      }
    )
    
    # Return the clustered Seurat object
    return(reactive({ values$clustered_seurat }))
  })
}

# Helper Functions

# Get available coloring options
getColoringOptions <- function(seurat_obj) {
  options <- c("cluster" = "cluster")
  
  if (!is.null(seurat_obj) && "sample" %in% colnames(seurat_obj@meta.data)) {
    options <- c("sample" = "sample", options)
  }
  
  # Always add gene option
  options <- c(options, "gene" = "gene")
  
  return(options)
}

# Search genes in Seurat object
searchGenes <- function(seurat_obj, query) {
  if (is.null(seurat_obj) || trimws(query) == "") {
    return(NULL)
  }
  
  # Convert query to lowercase for case-insensitive search
  query <- tolower(trimws(query))
  
  # Get gene mapping
  gene_mapping <- seurat_obj@misc$gene_mapping
  
  # Initialize empty results
  results <- NULL
  
  if (!is.null(gene_mapping)) {
    # Create a data frame of gene mappings for easier searching
    gene_df <- data.frame(
      ensembl_id = names(gene_mapping),
      gene_symbol = as.character(gene_mapping),
      stringsAsFactors = FALSE
    )
    
    # Remove any rows with NA values
    gene_df <- gene_df[!is.na(gene_df$ensembl_id) & !is.na(gene_df$gene_symbol), ]
    
    # Search by gene symbol (prioritize this)
    symbol_matches <- grep(query, tolower(gene_df$gene_symbol), value = FALSE)
    
    # Search by Ensembl ID if needed
    ensembl_matches <- grep(query, tolower(gene_df$ensembl_id), value = FALSE)
    
    # Combine matches (prioritize gene symbol matches)
    all_matches <- unique(c(symbol_matches, ensembl_matches))
    
    if (length(all_matches) > 0) {
      # Create results data frame
      results <- gene_df[all_matches, , drop = FALSE]
      
      # Only keep genes that are in the dataset
      results <- results[results$ensembl_id %in% rownames(seurat_obj), , drop = FALSE]
      
      # Sort by gene symbol
      if (nrow(results) > 0) {
        results <- results[order(results$gene_symbol), , drop = FALSE]
      } else {
        results <- NULL
      }
    }
  } else {
    # If no gene mapping is available, search directly in rownames
    matches <- grep(query, tolower(rownames(seurat_obj)), value = TRUE)
    
    if (length(matches) > 0) {
      results <- data.frame(
        ensembl_id = matches,
        gene_symbol = matches,
        stringsAsFactors = FALSE
      )
      results <- results[order(results$gene_symbol), , drop = FALSE]
    }
  }
  
  return(results)
}

# Run UMAP with specified dimensions and components
runUMAP <- function(seurat_obj, n_dims, n_components = 2, reduction.name = "umap", reduction.key = "UMAP_") {
  RunUMAP(seurat_obj, 
          dims = 1:n_dims, 
          n.components = n_components, 
          reduction.name = reduction.name,
          reduction.key = reduction.key)
}

# Create 3D UMAP plot
create3DUmapPlot <- function(seurat_obj, color_by = "cluster", reduction = "umap3d") {
  # Extract UMAP coordinates
  umap_data <- Embeddings(seurat_obj[[reduction]])
  
  # Extract coloring variable
  if (color_by == "cluster" && "seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    color_var <- as.factor(seurat_obj$seurat_clusters)
    plot_title <- "3D UMAP - Colored by Clusters"
  } else if (color_by == "sample" && "sample" %in% colnames(seurat_obj@meta.data)) {
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

# Create 3D Gene UMAP plot
create3DGeneUmapPlot <- function(seurat_obj, gene_id, reduction = "umap3d") {
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
umap_data <- Embeddings(seurat_obj[[reduction]])

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

# Get cluster statistics
getClusterStats <- function(seurat_obj) {
  # Check if clustering has been performed
  if (!("seurat_clusters" %in% colnames(seurat_obj@meta.data))) {
    return(data.frame(Error = "No clustering data available"))
  }
  
  # Get cluster identities
  clusters <- seurat_obj$seurat_clusters
  
  # Count cells per cluster
  cell_counts <- table(clusters)
  
  # Calculate percentage of cells per cluster
  cell_percentages <- prop.table(cell_counts) * 100
  
  # Create output data frame
  stats <- data.frame(
    Cluster = names(cell_counts),
    CellCount = as.numeric(cell_counts),
    Percentage = round(as.numeric(cell_percentages), 2)
  )
  
  # If sample information is available, calculate distribution by sample
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    # Create a contingency table of clusters by sample
    cluster_by_sample <- table(seurat_obj$seurat_clusters, seurat_obj$sample)
    
    # For each sample, add a column to the stats data frame
    for (sample_name in colnames(cluster_by_sample)) {
      sample_counts <- cluster_by_sample[, sample_name]
      stats[[paste0("Count_", sample_name)]] <- sample_counts
      
      # Calculate percentage within each cluster
      for (i in 1:nrow(stats)) {
        cluster_id <- stats$Cluster[i]
        cluster_total <- stats$CellCount[i]
        sample_count <- sample_counts[cluster_id]
        stats[[paste0("Pct_", sample_name)]][i] <- round((sample_count / cluster_total) * 100, 2)
      }
    }
  }
  
  return(stats)
}