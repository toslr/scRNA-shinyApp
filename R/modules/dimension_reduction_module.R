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
    uiOutput(ns("clusterSection")),
    uiOutput(ns("geneUmapSection"))
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
    
    # Gene expression UMAP section
    output$geneUmapSection <- renderUI({
      req(values$clustered_seurat)
      req("umap" %in% names(values$clustered_seurat@reductions))
      
      createGeneUmapUI(ns)
    })
    
    # Store gene search results
    gene_search_results <- reactiveVal(NULL)
    
    # Handle gene search
    observeEvent(input$searchGene, {
      req(values$clustered_seurat, input$geneQuery)
      
      # Get gene symbols from the Seurat object
      gene_mapping <- values$clustered_seurat@misc$gene_mapping
      
      # Make sure input is not empty
      if (trimws(input$geneQuery) == "") {
        gene_search_results(NULL)
        return()
      }
      
      # Convert query to lowercase for case-insensitive search
      query <- tolower(trimws(input$geneQuery))
      
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
          results <- results[results$ensembl_id %in% rownames(values$clustered_seurat), , drop = FALSE]
          
          # Sort by gene symbol
          if (nrow(results) > 0) {
            results <- results[order(results$gene_symbol), , drop = FALSE]
          } else {
            results <- NULL
          }
        }
      } else {
        # If no gene mapping is available, search directly in rownames
        matches <- grep(query, tolower(rownames(values$clustered_seurat)), value = TRUE)
        
        if (length(matches) > 0) {
          results <- data.frame(
            ensembl_id = matches,
            gene_symbol = matches,
            stringsAsFactors = FALSE
          )
          results <- results[order(results$gene_symbol), , drop = FALSE]
        }
      }
      
      # Update reactive value with results
      gene_search_results(results)
      
      # Log search results
      if (!is.null(results) && nrow(results) > 0) {
        print(paste("Found", nrow(results), "genes matching", input$geneQuery))
      } else {
        print(paste("No genes found matching", input$geneQuery))
      }
    })
    
    # Dynamic UI for gene selection
    output$geneSelectionUI <- renderUI({
      results <- gene_search_results()
      
      if (is.null(results) || nrow(results) == 0) {
        div(
          style = "margin-top: 32px;",
          p("No matching genes found. Try a different search term.")
        )
      } else {
        # Create selection with gene symbols as labels and Ensembl IDs as values
        # Prioritize showing gene symbols
        choices <- setNames(
          results$ensembl_id, 
          ifelse(
            !is.na(results$gene_symbol) & results$gene_symbol != results$ensembl_id,
            paste(results$gene_symbol, " (", results$ensembl_id, ")", sep = ""),
            results$ensembl_id
          )
        )
        
        selectInput(ns("selectedGene"), 
                    paste("Select from", nrow(results), "matching genes:"),
                    choices = choices,
                    selected = choices[1])
      }
    })
    
    output$selectedGeneLabel <- renderText({
      req(input$selectedGene)
      
      # Get gene symbol if we have gene mapping
      gene_id <- input$selectedGene
      
      if (!is.null(values$clustered_seurat@misc$gene_mapping) && 
          gene_id %in% names(values$clustered_seurat@misc$gene_mapping) &&
          !is.na(values$clustered_seurat@misc$gene_mapping[gene_id])) {
        gene_symbol <- values$clustered_seurat@misc$gene_mapping[gene_id]
        paste("Gene Expression UMAP:", gene_symbol, "(", gene_id, ")")
      } else {
        paste("Gene Expression UMAP:", gene_id)
      }
    })
    
    # Gene UMAP visualization section
    output$geneUmapVisSection <- renderUI({
      req(values$clustered_seurat, input$selectedGene)
      req("umap" %in% names(values$clustered_seurat@reductions))
      
      if (values$umap_type == "2D") {
        create2DGeneUmapUI(ns)
      } else {
        create3DGeneUmapUI(ns)
      }
    })
    
    # 2D Gene UMAP plot
    gene_umap_2d_plot <- reactive({
      req(values$clustered_seurat, input$selectedGene)
      req("umap" %in% names(values$clustered_seurat@reductions))
      req(values$umap_type == "2D")
      
      createGeneUmapPlot(values$clustered_seurat, input$selectedGene)
    })
    
    # 3D Gene UMAP plot
    gene_umap_3d_plot <- reactive({
      req(values$clustered_seurat, input$selectedGene)
      req("umap" %in% names(values$clustered_seurat@reductions))
      req(values$umap_type == "3D")
      
      createGeneUmap3DPlot(values$clustered_seurat, input$selectedGene)
    })
    
    # Render 2D gene UMAP plot
    output$geneUmapPlot2D <- renderPlot({
      req(input$selectedGene, values$umap_type == "2D")
      gene_umap_2d_plot()
    })
    
    # Render 3D gene UMAP plot
    output$geneUmapPlot3D <- renderPlotly({
      req(input$selectedGene, values$umap_type == "3D")
      gene_umap_3d_plot()
    })
    
    # Download handlers for gene UMAP plots
    output$downloadGeneUmap2DPlot <- downloadHandler(
      filename = function() {
        gene_name <- input$selectedGene
        paste("gene_umap_2d_", gene_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = gene_umap_2d_plot(), device = "png", width = 8, height = 8, dpi = 300)
      }
    )
    
    output$downloadGeneUmap3DPlot <- downloadHandler(
      filename = function() {
        gene_name <- input$selectedGene
        paste("gene_umap_3d_", gene_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html", sep = "")
      },
      content = function(file) {
        saveWidget(gene_umap_3d_plot(), file)
      }
    )
    
    # Return the clustered Seurat object
    return(reactive({ values$clustered_seurat }))
  })
}