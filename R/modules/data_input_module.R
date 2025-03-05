# R/modules/data_input_module.R

dataInputUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             h2(),
             strong("Select data directory"),
             shinyDirButton(ns("dir"), "Select Data Directory", "Choose directory"),
             verbatimTextOutput(ns("dirpath")),
             actionButton(ns("processData"), "Read Data")
      )
    )
  )
}

dataInputServer <- function(id, volumes = c(Home = '~/Desktop/Stanford/RA'), metadata_module) {
  moduleServer(id, function(input, output, session) {
    # Create reactive values
    seurat_obj <- reactiveVal(NULL)
    
    # Directory selection
    shinyDirChoose(input, 'dir', roots = volumes, session = session)
    
    selected_dir <- reactive({
      req(input$dir)
      parseDirPath(volumes, input$dir)
    })
    
    output$dirpath <- renderText({
      selected_dir()
    })
    
    # Data processing
    observeEvent(input$processData, {
      req(selected_dir())
      req(metadata_module$selectedSamples())
      
      selected_samples <- metadata_module$selectedSamples()
      print(paste("Processing samples:", paste(selected_samples, collapse=", ")))
      
      withProgress(message = 'Reading data...', value = 0, {
        tryCatch({
          # Get ENSEMBL to name mapping
          incProgress(0.1, detail = "Reading gene name conversion")
          gene_conversion <- read.csv("gene_conversion_results.csv")
          
          gene_mapping <- setNames(gene_conversion$external_gene_name, 
                                   gene_conversion$ensembl_gene_id)
          
          # Get all available files
          all_files <- list.files(selected_dir(), pattern="\\.txt\\.gz$")
          print(paste("Available files:", paste(all_files, collapse=", ")))
          
          # Only process selected GSMs
          selected_files <- character(0)
          file_to_gsm <- list()
          
          for(gsm in selected_samples) {
            pattern <- paste0("^", gsm)
            matches <- grep(pattern, all_files, value=TRUE)
            print(paste("Looking for files matching GSM:", gsm))
            print(paste("Found matches:", paste(matches, collapse=", ")))
            selected_files <- c(selected_files, matches)
            file_to_gsm[matches] <- gsm
          }
          
          if(length(selected_files) == 0) {
            stop("No matching files found for selected GSMs")
          }
          
          print(paste("Processing files:", paste(selected_files, collapse=", ")))
          
          # Initialize list for all Seurat objects
          seurat_objects <- list()
          count = 0
          
          # Process each file
          for(file in selected_files) {
            count = count + 1
            incProgress(0.9/(length(selected_files)+1), 
                        detail = paste("Reading file", count, "of", length(selected_files)))
            
            # Read the file
            data <- Read_GEO_Delim(data_dir = selected_dir(), 
                                   file_suffix = file)
            
            if(length(data) > 0 && !is.null(data[[1]])) {
              # Create Seurat object for this sample
              gsm <- file_to_gsm[[file]]
              print(paste("Creating Seurat object for GSM:", gsm))
              
              seurat <- CreateSeuratObject(counts = data[[1]], project = gsm)
              print("Created Seurat object")
              seurat$sample <- gsm
              
              # Add GEO metadata if available
              meta_func <- metadata_module$getMetadata
              
              if (!is.null(meta_func)) {
                metadata_data <- meta_func()
                sample_meta <- metadata_data[metadata_data$geo_accession == gsm, ]
                
                for(col in setdiff(colnames(sample_meta), "geo_accession")) {
                  seurat[[col]] <- sample_meta[[col]][1]
                }
              }
              
              # Add percent.mt
              seurat[["percent.mt"]] <- PercentageFeatureSet(seurat,
                                                             pattern = "^ENSMUSG00000064")
              
              # Store gene mapping
              seurat@misc$gene_mapping <- gene_mapping
              
              # Store in list
              seurat_objects[[gsm]] <- seurat
            }
          }
          
          # If multiple objects, merge them
          if(length(seurat_objects) > 1) {
            print("Merging Seurat objects...")
            gene_mapping_to_preserve <- seurat_objects[[1]]@misc$gene_mapping
            final_seurat <- merge(seurat_objects[[1]], 
                                  y = seurat_objects[2:length(seurat_objects)],
                                  add.cell.ids = names(seurat_objects))
            final_seurat@misc$gene_mapping <- gene_mapping_to_preserve
          } else {
            final_seurat <- seurat_objects[[1]]
          }
          
          # Store the object
          print("Storing final Seurat object...")
          seurat_obj(final_seurat)
          
        }, error = function(e) {
          print(paste("Error in data processing:", e$message))
          print(traceback())
        })
      })
    })
    
    # Trigger cleanup of downstream analyses
    #invalidateLater(100, session)
    
    # Return reactive value
    return(seurat_obj)
  })
}