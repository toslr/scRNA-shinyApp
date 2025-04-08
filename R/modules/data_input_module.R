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
             div(style = "margin-top: 10px;",
                 uiOutput(ns("speciesUI"))
             ),
             actionButton(ns("processData"), "Read Data"),
             div(style = "margin-top: 10px; color: #666;",
                 "Supported formats: .txt.gz, .h5, 10X MTX")
      )
    ),
    fluidRow(
      column(12,
             div(style = "margin-top: 10px; border: 1px solid #ddd; padding: 10px; border-radius: 4px; background-color: #f9f9f9;",
                 strong("Processing Status:"),
                 verbatimTextOutput(ns("statusOutput"))
             )
      )
    )
  )
}

dataInputServer <- function(id, volumes = c(Home = '~/Desktop/Stanford/RA'), metadata_module) {
  moduleServer(id, function(input, output, session) {
    # Create reactive values
    seurat_obj <- reactiveVal(NULL)
    processing_status <- reactiveVal("")
    species <- reactiveVal("auto")  # Default to auto detection
    
    # Add UI for species selection
    output$speciesUI <- renderUI({
      ns <- session$ns
      selectInput(ns("species"), "Species:",
                  choices = c("Auto-detect" = "auto",
                              "Human (Homo sapiens)" = "human",
                              "Mouse (Mus musculus)" = "mouse",
                              "Rat (Rattus norvegicus)" = "rat",
                              "Zebrafish (Danio rerio)" = "zebrafish",
                              "Fruit fly (D. melanogaster)" = "fly",
                              "C. elegans" = "worm"),
                  selected = "auto")
    })
    
    # Update species when input changes
    observeEvent(input$species, {
      species(input$species)
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
    
    # Status output
    output$statusOutput <- renderText({
      processing_status()
    })
    
    # Data processing
    observeEvent(input$processData, {
      req(selected_dir())
      req(metadata_module$selectedSamples())
      
      selected_samples <- metadata_module$selectedSamples()
      selected_species <- species()
      
      withProgress(message = 'Reading data...', value = 0, {
        tryCatch({
          # Get or create gene name mapping based on selected species
          incProgress(0.1, detail = "Setting up gene name conversion")
          processing_status("Setting up gene name conversion")
          
          # Use the gene mapping utility function
          if (selected_species == "auto" || selected_species == "") {
            # First try to use the existing gene mapping file
            if (file.exists("gene_conversion_results.csv")) {
              gene_conversion <- read.csv("gene_conversion_results.csv")
              gene_mapping <- setNames(gene_conversion$external_gene_name, 
                                       gene_conversion$ensembl_gene_id)
              processing_status("Using existing gene mapping file")
            } else {
              # Default to mouse if no specific species
              gene_mapping <- get_gene_mapping("mouse")
              processing_status("Using default mouse gene mapping")
            }
          } else {
            # Use species-specific mapping
            gene_mapping <- get_gene_mapping(selected_species)
            processing_status(paste("Using", selected_species, "gene mapping"))
          }
          
          # Get all available files
          processing_status("Looking for files in the selected directory")
          all_files <- list.files(selected_dir())
          
          # Only process selected GSMs
          selected_files <- character(0)
          file_to_gsm <- list()
          file_formats <- list()
          
          for(gsm in selected_samples) {
            processing_status(paste("Detecting file format for", gsm))
            
            # Determine file format type
            format_type <- Detect_File_Format_Type(selected_dir(), gsm)
            
            if (format_type == "unknown") {
              processing_status(paste("Warning: No supported files found for", gsm))
              next
            }
            
            processing_status(paste("Found", format_type, "format for", gsm))
            
            if (format_type == "txt.gz") {
              # For text files, find matching files
              pattern <- paste0("^", gsm)
              matches <- grep(pattern, all_files, value=TRUE)
              
              for(match in matches) {
                selected_files <- c(selected_files, match)
                file_to_gsm[[match]] <- gsm
                file_formats[[match]] <- format_type
              }
            } else {
              # For H5 or MTX formats, just use the GSM ID
              selected_files <- c(selected_files, gsm)
              file_to_gsm[[gsm]] <- gsm
              file_formats[[gsm]] <- format_type
            }
          }
          
          processing_status(paste("Found", length(selected_files), "files for selected GSMs"))
          
          if(length(selected_files) == 0) {
            processing_status("No matching files found for selected GSMs")
            stop("No matching files found for selected GSMs")
          }
          
          # Initialize list for all Seurat objects
          seurat_objects <- list()
          count = 0
          detected_species <- NULL
          
          # Process each file
          for(file in selected_files) {
            count = count + 1
            format_type <- file_formats[[file]]
            processing_status(paste("Reading file", count, "of", length(selected_files), ":", file, "(", format_type, "format)"))
            incProgress(0.9/(length(selected_files)+1), 
                        detail = paste("Reading file", count, "of", length(selected_files)))
            
            data <- Read_Data_File(selected_dir(), file)
            
            if(length(data) > 0 && !is.null(data[[1]])) {
              gsm <- file_to_gsm[[file]]
              processing_status(paste("Creating Seurat object for sample", gsm))
              seurat <- CreateSeuratObject(counts = data[[1]], project = gsm)
              seurat$sample <- gsm
              
              # Detect species if auto mode
              if (selected_species == "auto" && is.null(detected_species)) {
                processing_status("Auto-detecting species from gene IDs...")
                detected_species <- detect_species(seurat)
                processing_status(paste("Detected species:", detected_species))
                # Update selected species
                selected_species <- detected_species
              }
              
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
              
              # Add species information
              seurat$species <- selected_species
              
              # Calculate mitochondrial percentage with species info
              seurat <- Calculate_MT_Percent(seurat, species = selected_species)
              
              seurat@misc$gene_mapping <- gene_mapping
              seurat_objects[[gsm]] <- seurat
            }
          }
          
          processing_status(paste("Created", length(seurat_objects), "Seurat objects"))
          
          # If multiple objects, merge them
          if(length(seurat_objects) > 1) {
            processing_status("Merging multiple Seurat objects")
            gene_mapping_to_preserve <- seurat_objects[[1]]@misc$gene_mapping
            species_to_preserve <- seurat_objects[[1]]$species[1]
            
            final_seurat <- merge(seurat_objects[[1]], 
                                  y = seurat_objects[2:length(seurat_objects)],
                                  add.cell.ids = names(seurat_objects))
            
            final_seurat@misc$gene_mapping <- gene_mapping_to_preserve
            final_seurat$species <- species_to_preserve
          } else if(length(seurat_objects) == 1) {
            processing_status("Using single Seurat object")
            final_seurat <- seurat_objects[[1]]
          } else {
            processing_status("No Seurat objects were created")
            stop("No Seurat objects were created")
          }
          
          # Store the object
          processing_status("Data processing completed successfully")
          seurat_obj(final_seurat)
          
        }, error = function(e) {
          processing_status(paste("Error in data processing:", e$message))
          showNotification(paste("Error in data processing:", e$message), type = "error")
        })
      })
    })
    
    # Return reactive value and processing status
    return(list(
      data = seurat_obj,
      status = processing_status,
      species = species
    ))
  })
}