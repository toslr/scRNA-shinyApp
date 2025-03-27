#' @title Data Input Module UI
#' @description Creates the UI for data input which allows users to select and
#'   process scRNA-seq data directories.
#' @param id The module ID
#' @return A Shiny UI element containing the data input interface
#' @export
dataInputUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             strong("Select data directory"),
             shinyDirButton(ns("dir"), "Select Data Directory", "Choose directory"),
             verbatimTextOutput(ns("dirpath")),
             actionButton(ns("processData"), "Read Data")
      )
    )
  )
}

#' @title Data Input Module Server
#' @description Server logic for the data input module which loads and processes 
#'   scRNA-seq data files into Seurat objects.
#' @param id The module ID
#' @param volumes Volumes configuration for directory selection, defaults to Home desktop
#' @param metadata_module Metadata module instance to access selected samples and annotations
#' @return A reactive expression containing the processed Seurat object
#' @export
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
      
      withProgress(message = 'Reading data...', value = 0, {
        tryCatch({
          # Get ENSEMBL to name mapping
          incProgress(0.1, detail = "Reading gene name conversion")
          gene_conversion <- read.csv("gene_conversion_results.csv")
          
          gene_mapping <- setNames(gene_conversion$external_gene_name, 
                                   gene_conversion$ensembl_gene_id)
          
          # Get all available files
          all_files <- list.files(selected_dir(), pattern="\\.txt\\.gz$")
          
          # Only process selected GSMs
          selected_files <- character(0)
          file_to_gsm <- list()
          for(gsm in selected_samples) {
            pattern <- paste0("^", gsm)
            matches <- grep(pattern, all_files, value=TRUE)
            selected_files <- c(selected_files, matches)
            for(match in matches) {
              file_to_gsm[[match]] <- gsm
            }
          }
          
          if(length(selected_files) == 0) {
            stop("No matching files found for selected GSMs")
          }
          
          # Initialize list for all Seurat objects
          seurat_objects <- list()
          count = 0
          
          # Process each file
          for(file in selected_files) {
            count = count + 1
            incProgress(0.9/(length(selected_files)+1), 
                        detail = paste("Reading file", count, "of", length(selected_files)))
            
            data <- Read_GEO_Delim(data_dir = selected_dir(), 
                                   file_suffix = file)
            
            if(length(data) > 0 && !is.null(data[[1]])) {
              gsm <- file_to_gsm[[file]]
              seurat <- CreateSeuratObject(counts = data[[1]], project = gsm)
              seurat$sample <- gsm
              
              # Add GEO metadata if available
              meta_func <- metadata_module$getMetadata
              if (!is.null(meta_func)) {
                metadata_data <- meta_func()
                if (!is.null(metadata_data)) {
                  sample_meta <- metadata_data[metadata_data$geo_accession == gsm, ]
                  for(col in setdiff(colnames(sample_meta), "geo_accession")) {
                    seurat[[col]] <- sample_meta[[col]][1]
                  }
                }
              }
              
              seurat[["percent.mt"]] <- PercentageFeatureSet(seurat,
                                                             pattern = "^ENSMUSG00000064")
              seurat@misc$gene_mapping <- gene_mapping
              seurat_objects[[gsm]] <- seurat
            }
          }
          
          # If multiple objects, merge them
          if(length(seurat_objects) > 1) {
            gene_mapping_to_preserve <- seurat_objects[[1]]@misc$gene_mapping
            final_seurat <- merge(seurat_objects[[1]], 
                                  y = seurat_objects[2:length(seurat_objects)],
                                  add.cell.ids = names(seurat_objects))
            final_seurat@misc$gene_mapping <- gene_mapping_to_preserve
          } else if(length(seurat_objects) == 1) {
            final_seurat <- seurat_objects[[1]]
          } else {
            stop("No Seurat objects were created")
          }
          
          # Store the object
          seurat_obj(final_seurat)
          
        }, error = function(e) {
          showNotification(paste("Error in data processing:", e$message), type = "error")
        })
      })
    })
    
    # Return reactive value
    return(seurat_obj)
  })
}