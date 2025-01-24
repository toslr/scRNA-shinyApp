library(shiny)
library(patchwork)
library(shinyFiles)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scCustomize)


options(shiny.maxRequestSize = 500*1024^2)

ui <- fluidPage(
  titlePanel("Single-Cell RNA Analysis"),
  sidebarLayout(
    sidebarPanel(
      textOutput("intro"),
      shinyDirButton("dir", "Select Data Directory", "Choose directory"),
      verbatimTextOutput("dirpath"),
      actionButton("processData", "Read Data")
    ),
    mainPanel(
      textOutput("qcText"),
      plotOutput("qcPlot", height = "600px"),
      plotOutput("elbowPlot", height = "400px"),
      plotOutput("umapPlot", height = "400px"),
      plotOutput("clusterPlot", height = "400px"),
      #selectInput("cluster1", "Select first cluster:", choices = NULL),
      #selectInput("cluster2", "Select second cluster:", choices = NULL),
      selectInput("targetCluster", "Select cluster to compare against all others:", choices = NULL),
      actionButton("runDE", "Run Differential Expression"),
      plotOutput("volcanoPlot", height = "400px"),
      DT::dataTableOutput("deTable")
    )
  )
)

server <- function(input, output, session) {
  output$intro <- renderText({"This tool is able to read .txt.gz files, with ENSEMBL gene names as rows and cell names as columns. Select the directory containing those files and start the analysis."})

  volumes <- c(Home = '~/Desktop/Stanford/RA')
  shinyDirChoose(input, 'dir', roots = volumes, session = session)
  
  selected_dir <- reactive({
    req(input$dir)
    parseDirPath(volumes, input$dir)
  })
  
  output$dirpath <- renderText({
    selected_dir()
  })
  
  seurat_data <- reactiveVal()
  
  observeEvent(input$processData, {
    req(selected_dir())
    withProgress(message = 'Reading data...', value = 0, {
      # Read gene conversion file first
      incProgress(0.1, detail = "Reading gene name conversion")
      gene_conversion <- read.csv(file.path(selected_dir(), "gene_conversion_results.csv"))
      gene_mapping <- setNames(gene_conversion$external_gene_name, 
                               gene_conversion$ensembl_gene_id)
      
      files <- list.files(selected_dir(), pattern="1.txt.gz$")
      n_files <- length(files)
      
      incProgress(0.2, detail = "Reading expression data")
      GEO_data <- Read_GEO_Delim(data_dir = selected_dir(), 
                                 file_suffix = '1.txt.gz')
      
      #for (i in seq_along(GEO_data)) {
      #  ensembl_names <- rownames(GEO_data[[i]])
      #  converted_names <- gene_mapping[ensembl_names]
        
      #  # Replace NA and empty values with original names
      #  invalid_names <- is.na(converted_names) | 
      #    converted_names == "" | 
      #    grepl("^\\s*$", converted_names)
      #  converted_names[invalid_names] <- ensembl_names[invalid_names]
      #  # Update matrix with clean row names
      #  GEO_data[[i]] <- GEO_data[[i]][ensembl_names,]
      #  rownames(GEO_data[[i]]) <- converted_names
      #}
      
      incProgress(0.4, detail = "Creating Seurat object")
      seurat <- CreateSeuratObject(counts = GEO_data, project = "DS1")
      seurat@misc$gene_mapping <- gene_mapping
      incProgress(0.6, detail = "Converting gene names")
      gene_names <- rownames(GEO_data[[1]])
      #converted_names <- gene_mapping[gene_names]
      # Keep ENSEMBL ID if no match found
      #converted_names[is.na(converted_names)] <- gene_names[is.na(converted_names)]
      
      count_matrices <- names(seurat@assays$RNA@layers)
      
      for (i in 1:length(count_matrices)) {
        matrix_name <- count_matrices[i]
        cell_names <- colnames(GEO_data[[i]])
        rownames(seurat@assays$RNA@layers[[matrix_name]]) <- gene_names #converted_names
        colnames(seurat@assays$RNA@layers[[matrix_name]]) <- cell_names
      #  incProgress(0.2/length(count_matrices), 
      #              detail = sprintf("Processing matrix %d of %d", i, length(count_matrices)))
      }
      
      sample_names <- gsub("counts\\.GSM[0-9]+_(.*)", "\\1", count_matrices)
      cell_sample_info <- character(0)
      for (i in 1:length(count_matrices)) {
        matrix_name <- count_matrices[i]
        n_cells <- ncol(seurat@assays$RNA@layers[[matrix_name]])
        cell_sample_info <- c(cell_sample_info, 
                              rep(sample_names[i], n_cells))
      }
      seurat$sample <- cell_sample_info
      
      seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, 
                                                     pattern = "^ENSMUSG00000064")
      
      seurat_data(seurat)
    })
  })
  
  output$qcText <- renderText({req(input$processData)
    "Quality control"})
  output$qcPlot <- renderPlot({
    req(seurat_data())
    VlnPlot(seurat_data(), 
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
            group.by = "sample",
            ncol = 3)
  })
  
  # Add new UI elements after initial plot is rendered
  observeEvent(seurat_data(), {
    insertUI(
      selector = "#qcPlot",
      where = "afterEnd",
      ui = {
        tagList(
          renderPrint("Please adjust filtering parameters:"),
          numericInput("minFeature", "Minimum Features:", 500),
          numericInput("maxFeature", "Maximum Features:", 5000),
          numericInput("maxMT", "Maximum MT %:", 5),
          actionButton("processSeurat", "Filter and run PCA")
        )
      }
    )
  })
  
  # New reactive value for processed data
  processed_seurat <- reactiveVal()
  
  # Process Seurat object when button is clicked
  observeEvent(input$processSeurat, {
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
      processed_seurat(seurat)
    })
  })
  
  observeEvent(processed_seurat(), {
    removeUI(selector = "#dimension_controls")
    insertUI(
      selector = "#elbowPlot",
      where = "afterEnd",
      ui = tags$div(
        id = "dimension_controls",
        renderPrint("Please adjust the number of PC for reduction"),
        numericInput("nDims", "Number of dimensions for UMAP:", 10),
        actionButton("runUMAP", "Run UMAP")
      )
    )
  })
  
  # Only show clustering controls after UMAP is run
  observeEvent(input$runUMAP, {
    req(processed_seurat(), input$nDims)
    withProgress(message = 'Computing UMAP...', {
      seurat <- processed_seurat()
      seurat <- RunUMAP(seurat, dims = 1:input$nDims)
      processed_seurat(seurat)
    })
    removeUI(selector = "#clustering_controls")
    insertUI(
      selector = "#umapPlot",
      where = "afterEnd",
      ui = tags$div(
        id = "clustering_controls",
        renderPrint("Please adjust clustering resolution:"),
        numericInput("resolution", "Clustering Resolution:", 0.5, min = 0, max = 2, step = 0.01),
        actionButton("runClustering", "Run Clustering")
      )
    )
  })
  
  observeEvent(input$runClustering, {
    req(processed_seurat(), input$nDims, input$resolution)
    seurat <- processed_seurat()
    seurat <- FindNeighbors(seurat, dims = 1:input$nDims)
    seurat <- FindClusters(seurat, resolution = input$resolution)
    processed_seurat(seurat)
  })
  
  # Add elbow plot output
  output$elbowPlot <- renderPlot({
    req(processed_seurat())
    ElbowPlot(processed_seurat(), 
              ndims = ncol(Embeddings(processed_seurat(), "pca")))
  })
  
  # Add UMAP plot output
  output$umapPlot <- renderPlot({
    req(processed_seurat())
    req("umap" %in% names(processed_seurat()@reductions))
    isolate({
      DimPlot(processed_seurat(), reduction = "umap")
    })
  })
  
  # Add clustering plot output
  output$clusterPlot <- renderPlot({
    req(processed_seurat())
    req("umap" %in% names(processed_seurat()@reductions))
    DimPlot(processed_seurat(), reduction = "umap", label = TRUE)
  })
  
  # Update cluster choices after clustering
  observe({
    req(processed_seurat())
    if ("seurat_clusters" %in% colnames(processed_seurat()@meta.data)) {
      clusters <- sort(unique(processed_seurat()$seurat_clusters))
      updateSelectInput(session, "targetCluster", choices = clusters)
      #updateSelectInput(session, "cluster1", choices = clusters)
      #updateSelectInput(session, "cluster2", choices = clusters)
    }
  })
  
  # Run DE analysis
  observeEvent(input$runDE, {
    req(processed_seurat(), input$targetCluster)
    #req(processed_seurat(), input$cluster1, input$cluster2)
    withProgress(message = 'Computing differential expression...', {
      de_results <- FindMarkers(processed_seurat(),
                                ident.1 = input$targetCluster,
                                ident.2 = NULL,
                                #ident.1 = input$cluster1,
                                #ident.2 = input$cluster2,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
      
      # Add gene names as column
      de_results$gene <- processed_seurat()@misc$gene_mapping[rownames(de_results)]
      
      # Store results
      de_genes(de_results)
    })
  })
  
  # Create reactive value for DE results
  de_genes <- reactiveVal()
  
  # Render volcano plot
  output$volcanoPlot <- renderPlot({
    req(de_genes())
    ggplot(de_genes(), aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point() +
      theme_classic() +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed")
  })
  
  # Render results table
  output$deTable <- DT::renderDataTable({
    req(de_genes())
    DT::datatable(de_genes(),
                  options = list(pageLength = 10),
                  rownames = FALSE)
  })
  
}

shinyApp(ui = ui, server = server)