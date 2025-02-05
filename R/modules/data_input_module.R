# R/modules/data_input_module.R

dataInputUI <- function(id) {
  ns <- NS(id)
  tagList(
    textOutput(ns("intro")),
    # Add GEO input field
    fluidRow(
      column(12,
             strong("Please input your GEO Series ID and browse to find your expression files."),
             h2(),
             textInput(ns("geoID"), "GEO Series ID (e.g., GSE182846)", 
                       placeholder = "GSExxxxxx"),
             actionButton(ns("fetchGEO"), "Fetch GEO Metadata"),
             verbatimTextOutput(ns("geoStatus")),
             h2(),
             strong("Select data directory"),
             shinyDirButton(ns("dir"), "Select Data Directory", "Choose directory"),
             verbatimTextOutput(ns("dirpath")),
             actionButton(ns("processData"), "Read Data")
      )
    )
  )
}

dataInputServer <- function(id, volumes = c(Home = '~/Desktop/Stanford/RA')) {
  moduleServer(id, function(input, output, session) {
    # Create reactive values
    geo_metadata <- reactiveVal(NULL)
    seurat_obj <- reactiveVal(NULL)
    selected_samples <- reactiveVal(NULL)
    
    # GEO metadata fetching
    observeEvent(input$fetchGEO, {
      req(input$geoID)
      tryCatch({
        withProgress(message = 'Fetching GEO metadata...', value = 0, {
          # Fetch GEO data
          gset <- getGEO(input$geoID, GSEMatrix = TRUE)
          if (is(gset, "list")) gset <- gset[[1]]
          
          # Extract phenoData
          pheno_data <- pData(phenoData(gset))
          
          # Create basic metadata dataframe with just geo_accession and title
          metadata <- data.frame(
            geo_accession = rownames(pheno_data),
            title = pheno_data$title,
            stringsAsFactors = FALSE
          )
          
          # Get all characteristics columns (ch1, ch1.1, ch2, etc)
          char_columns <- grep("characteristics_ch", colnames(pheno_data), value = TRUE)
          
          # Add each characteristics column as is
          for (col in char_columns) {
            metadata[[col]] <- pheno_data[[col]]
          }
          
          # Store metadata and initialize all samples as selected
          geo_metadata(metadata)
          selected_samples(metadata$geo_accession)
        })
        
        output$geoStatus <- renderText({
          "GEO metadata successfully fetched!"
        })
        
      }, error = function(e) {
        output$geoStatus <- renderText({
          paste("Error fetching GEO data:", e$message)
        })
      })
    })
    
    # Update selected_samples when checkboxes change
    observeEvent(input$selectedSamples, {
      selected_samples(input$selectedSamples)
    })
    
    # Directory selection
    shinyDirChoose(input, 'dir', roots = volumes, session = session)
    
    selected_dir <- reactive({
      req(input$dir)
      parseDirPath(volumes, input$dir)
    })
    
    output$dirpath <- renderText({
      selected_dir()
    })
    
    # Data processing with sample selection
    observeEvent(input$processData, {
      req(selected_dir(), selected_samples())
      req(length(selected_samples()) > 0)
      
      print("Starting data processing...")
      
      withProgress(message = 'Reading data...', value = 0, {
        tryCatch({
          incProgress(0.1, detail = "Reading gene name conversion")
          gene_conversion <- read.csv(file.path(selected_dir(), "gene_conversion_results.csv"))
          
          gene_mapping <- setNames(gene_conversion$external_gene_name, 
                                   gene_conversion$ensembl_gene_id)
          
          # Filter files based on selected samples
          all_files <- list.files(selected_dir(), pattern="1.txt.gz$")
          selected_gsm_pattern <- paste(selected_samples(), collapse="|")
          files <- grep(selected_gsm_pattern, all_files, value=TRUE)
          
          print(paste("Processing selected files:", paste(files, collapse=", ")))
          
          incProgress(0.2, detail = "Reading expression data")
          GEO_data <- Read_GEO_Delim(data_dir = selected_dir(), 
                                     file_pattern = paste0("(", selected_gsm_pattern, ").*1.txt.gz$"))
          
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
          
          if (!is.null(geo_metadata())) {
            metadata <- geo_metadata()
            gsm_numbers <- unique(gsub("counts\\.(GSM[0-9]+)_.*", "\\1", count_matrices))
            sample_metadata <- metadata[metadata$geo_accession %in% gsm_numbers, ]
            
            for (col in colnames(sample_metadata)) {
              if (col != "geo_accession") {
                meta_values <- setNames(
                  sample_metadata[[col]],
                  gsub("counts\\.(GSM[0-9]+)_.*", "\\1", count_matrices)
                )
                
                cell_meta <- sapply(
                  gsub("counts\\.(GSM[0-9]+)_.*", "\\1", count_matrices),
                  function(x) meta_values[x]
                )
                
                seurat[[col]] <- unname(cell_meta[match(seurat$sample, names(cell_meta))])
              }
            }
          }
          
          seurat[["percent.mt"]] <- PercentageFeatureSet(seurat,
                                                         pattern = "^ENSMUSG00000064")
          seurat_obj(seurat)
          
        }, error = function(e) {
          print(paste("Error in data processing:", e$message))
        })
      })
    })
    
    # Return both reactive values
    return(reactive({
      list(
        seurat = seurat_obj(),
        metadata = geo_metadata()
      )
    }))
  })
}