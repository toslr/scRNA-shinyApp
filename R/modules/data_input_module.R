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

#' @title Data Input Module Server
#' @description Server logic for the data input module that handles directory selection 
#'   and data loading for scRNA-seq analysis.
#' @param id The module ID
#' @param volumes Named vector of root directories to use for file selection (default uses system defaults)
#' @param metadata_module Metadata module instance for accessing metadata information
#' @return A list of reactive expressions for accessing the loaded data
#' @export
dataInputServer <- function(id, volumes = NULL, metadata_module) {
  moduleServer(id, function(input, output, session) {
    # Create reactive values
    seurat_obj <- reactiveVal(NULL)
    processing_status <- reactiveVal("")
    species <- reactiveVal("auto")  # Default to auto detection
    
    # Define default root volumes if none are provided
    if (is.null(volumes)) {
      # Create a list of root locations that are likely to exist across systems
      volumes <- c(
        Home = Sys.getenv("HOME"),
        Desktop = file.path(Sys.getenv("HOME"), "Desktop"),
        Documents = file.path(Sys.getenv("HOME"), "Documents"),
        Downloads = file.path(Sys.getenv("HOME"), "Downloads")
      )
      
      # On Windows, also add common drive letters
      if (.Platform$OS.type == "windows") {
        for (drive in c("C:", "D:", "E:")) {
          if (dir.exists(drive)) {
            volumes[drive] <- drive
          }
        }
      }
      
      # Filter out paths that don't exist to avoid errors
      volumes <- volumes[sapply(volumes, dir.exists)]
      
      # If no valid paths were found, default to current working directory
      if (length(volumes) == 0) {
        volumes <- c("Current" = getwd())
      }
    }
    
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
    
    # Add a function to provide available root directories for UI
    getRootDirectories <- reactive({
      volumes
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
          
          # Only process selected samples - enhanced for custom metadata
          selected_files <- character(0)
          file_to_gsm <- list()
          file_formats <- list()
          
          for(sample_id in selected_samples) {
            processing_status(paste("Detecting file format for", sample_id))
            
            # Determine file format type
            format_type <- Detect_File_Format_Type(selected_dir(), sample_id)
            
            if (format_type == "unknown") {
              # Try checking for files that might contain the sample ID
              matching_files <- grep(sample_id, all_files, value=TRUE)
              
              if (length(matching_files) > 0) {
                # Check if any of the matching files have a known format
                for (match_file in matching_files) {
                  match_format <- Detect_File_Format(match_file)
                  if (match_format != "unknown") {
                    format_type <- match_format
                    processing_status(paste("Found", format_type, "format for", sample_id))
                    
                    selected_files <- c(selected_files, match_file)
                    file_to_gsm[[match_file]] <- sample_id
                    file_formats[[match_file]] <- format_type
                    break
                  }
                }
              }
              
              if (format_type == "unknown") {
                processing_status(paste("Warning: No supported files found for", sample_id))
                next
              }
            } else {
              processing_status(paste("Found", format_type, "format for", sample_id))
              
              if (format_type == "txt.gz") {
                # For text files, find matching files
                pattern <- paste0("^", sample_id)
                matches <- grep(pattern, all_files, value=TRUE)
                
                for(match in matches) {
                  selected_files <- c(selected_files, match)
                  file_to_gsm[[match]] <- sample_id
                  file_formats[[match]] <- format_type
                }
              } else {
                # For H5 or MTX formats, just use the sample ID
                selected_files <- c(selected_files, sample_id)
                file_to_gsm[[sample_id]] <- sample_id
                file_formats[[sample_id]] <- format_type
              }
            }
          }
          
          processing_status(paste("Found", length(selected_files), "files for selected samples"))
          
          if(length(selected_files) == 0) {
            processing_status("No matching files found for selected samples")
            stop("No matching files found for selected samples")
          }
          
          # Initialize list for all Seurat objects
          seurat_objects <- list()
          count = 0
          detected_species <- NULL
          
          # Process each file using createSeuratObjects
          seurat_objects <- createSeuratObjects(
            selected_dir(), 
            selected_files, 
            file_to_gsm, 
            file_formats, 
            gene_mapping, 
            metadata_module, 
            selected_species, 
            processing_status, 
            incProgress
          )
          
          processing_status(paste("Created", length(seurat_objects), "Seurat objects"))
          
          # If multiple objects, merge them
          if(length(seurat_objects) > 1) {
            processing_status("Merging multiple Seurat objects")
            
            # Preserve important information from first object
            gene_mapping_to_preserve <- seurat_objects[[1]]@misc$gene_mapping
            species_to_preserve <- seurat_objects[[1]]$species[1]
            
            # Note: cell names already have GSM prefixes, so we don't need add.cell.ids
            # Using merge with proper checking to avoid cell name conflicts
            final_seurat <- merge(
              x = seurat_objects[[1]],
              y = seurat_objects[2:length(seurat_objects)],
              add.cell.ids = NULL, # Cell IDs already have GSM prefixes
              project = "merged_samples"
            )
            
            # Restore preserved information
            final_seurat@misc$gene_mapping <- gene_mapping_to_preserve
            final_seurat$species <- species_to_preserve
            
            processing_status(paste("Successfully merged", length(seurat_objects), 
                                    "samples into one Seurat object with", 
                                    ncol(final_seurat), "cells and",
                                    nrow(final_seurat), "genes"))
            
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
          # Print complete error message for debugging
          print(e)
          processing_status(paste("Error in data processing:", e$message))
          showNotification(paste("Error in data processing:", e$message), type = "error", duration = NULL)
        })
      })
    })
    
    # Return reactive value and processing status
    return(list(
      data = seurat_obj,
      status = processing_status,
      species = species,
      getRootDirectories = getRootDirectories
    ))
  })
}

#' @title Create Seurat Objects from Files
#' @description Creates Seurat objects from data files, adding metadata and gene mappings.
#' @param dir_path Path to the directory containing data files
#' @param selected_files Character vector of selected file names to process
#' @param file_to_gsm Mapping of file names to sample IDs
#' @param file_formats Mapping of file names to format types
#' @param gene_mapping Named vector for mapping Ensembl IDs to gene names
#' @param metadata_module Metadata module instance for accessing metadata information
#' @param selected_species Selected species for mitochondrial gene detection
#' @param processing_status Function to update processing status
#' @param incProgress Function to update progress
#' @return List of Seurat objects, one for each sample
#' @export
createSeuratObjects <- function(dir_path, selected_files, file_to_gsm, file_formats, 
                                gene_mapping, metadata_module, selected_species,
                                processing_status, incProgress) {
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
    
    # Read the data file
    data <- Read_Data_File(dir_path, file)
    
    if(length(data) > 0 && !is.null(data[[1]])) {
      sample_id <- file_to_gsm[[file]]
      processing_status(paste("Creating Seurat object for sample", sample_id))
      
      # Create a Seurat object with proper cell naming to avoid conflicts
      # Add sample prefix to cell names to ensure uniqueness across samples
      cell_names <- colnames(data[[1]])
      new_cell_names <- paste0(sample_id, "_", cell_names)
      colnames(data[[1]]) <- new_cell_names
      
      # Create Seurat object with the modified cell names
      seurat <- CreateSeuratObject(counts = data[[1]], project = sample_id)
      
      # Add sample identifier
      seurat$sample <- sample_id
      
      # Detect species if auto mode
      if (selected_species == "auto" && is.null(detected_species)) {
        processing_status("Auto-detecting species from gene IDs...")
        detected_species <- detect_species(seurat)
        processing_status(paste("Detected species:", detected_species))
        # Update selected species
        selected_species <- detected_species
      }
      
      # Add metadata if available
      meta_func <- metadata_module$getMetadata
      if (!is.null(meta_func)) {
        metadata_data <- meta_func()
        if (!is.null(metadata_data)) {
          # Try to find metadata row matching this sample
          sample_meta <- NULL
          
          # First check geo_accession column
          if ("geo_accession" %in% colnames(metadata_data)) {
            sample_meta <- metadata_data[metadata_data$geo_accession == sample_id, ]
          }
          
          # If not found, try sample_id column
          if ((is.null(sample_meta) || nrow(sample_meta) == 0) && 
              "sample_id" %in% colnames(metadata_data)) {
            sample_meta <- metadata_data[metadata_data$sample_id == sample_id, ]
          }
          
          # If we found metadata, add all columns to the Seurat object
          if (!is.null(sample_meta) && nrow(sample_meta) > 0) {
            for(col in colnames(sample_meta)) {
              # Skip special columns that would conflict with existing data
              if (!col %in% c("sample")) {
                seurat[[col]] <- sample_meta[[col]][1]
              }
            }
          }
        }
      }
      
      # Add species information
      seurat$species <- selected_species
      
      # Calculate mitochondrial percentage with species info
      seurat <- Calculate_MT_Percent(seurat, species = selected_species)
      
      seurat@misc$gene_mapping <- gene_mapping
      seurat_objects[[sample_id]] <- seurat
      
      processing_status(paste("Successfully created Seurat object for", sample_id,
                              "with", ncol(seurat), "cells and", nrow(seurat), "genes"))
    } else {
      processing_status(paste("Warning: No data found for", file))
    }
  }
  
  return(seurat_objects)
}