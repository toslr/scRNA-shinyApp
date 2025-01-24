dataInputUI <- function(id) {
  ns <- NS(id)
  tagList(
    textOutput(ns("intro")),
    shinyDirButton(ns("dir"), "Select Data Directory", "Choose directory"),
    verbatimTextOutput(ns("dirpath")),
    actionButton(ns("processData"), "Read Data")
  )
}

dataInputServer <- function(id, volumes = c(Home = '~/Desktop/Stanford/RA')) {
  moduleServer(id, function(input, output, session) {
    output$intro <- renderText({
      "This tool is able to read .txt.gz files, with ENSEMBL gene names as rows and cell names as columns. Select the directory containing those files and start the analysis."
    })
    
    shinyDirChoose(input, 'dir', roots = volumes, session = session)
    
    selected_dir <- reactive({
      req(input$dir)
      parseDirPath(volumes, input$dir)
    })
    
    output$dirpath <- renderText({
      selected_dir()
    })
    
    seurat_data <- eventReactive(input$processData, {
      req(selected_dir())
      withProgress(message = 'Reading data...', value = 0, {
        # Data processing code from original server function
        incProgress(0.1, detail = "Reading gene name conversion")
        gene_conversion <- read.csv(file.path(selected_dir(), "gene_conversion_results.csv"))
        gene_mapping <- setNames(gene_conversion$external_gene_name, 
                                 gene_conversion$ensembl_gene_id)
        
        files <- list.files(selected_dir(), pattern="1.txt.gz$")
        
        incProgress(0.2, detail = "Reading expression data")
        GEO_data <- Read_GEO_Delim(data_dir = selected_dir(), 
                                   file_suffix = '1.txt.gz')
        
        incProgress(0.4, detail = "Creating Seurat object")
        seurat <- CreateSeuratObject(counts = GEO_data, project = "DS1")
        seurat@misc$gene_mapping <- gene_mapping
        
        incProgress(0.6, detail = "Converting gene names")
        gene_names <- rownames(GEO_data[[1]])
        count_matrices <- names(seurat@assays$RNA@layers)

        for (i in 1:length(count_matrices)) {
          matrix_name <- count_matrices[i]
          cell_names <- colnames(GEO_data[[i]])
          rownames(seurat@assays$RNA@layers[[matrix_name]]) <- gene_names
          colnames(seurat@assays$RNA@layers[[matrix_name]]) <- cell_names
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
        seurat
        #seurat_data(seurat)
      })
    })
    
    return(seurat_data)
  })
}