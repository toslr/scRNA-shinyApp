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
    
    # Modified data processing part to handle multiple files without joining layers
    observeEvent(input$processData, {
      req(selected_dir(), input$selectedSamples)
      req(length(input$selectedSamples) > 0)
      
      print("Starting data processing...")
      
      withProgress(message = 'Reading data...', value = 0, {
        tryCatch({
          incProgress(0.1, detail = "Reading gene name conversion")
          gene_conversion <- read.csv(file.path(selected_dir(), "gene_conversion_results.csv"))
          
          gene_mapping <- setNames(gene_conversion$external_gene_name, 
                                   gene_conversion$ensembl_gene_id)
          
          # First get all available files
          all_files <- list.files(selected_dir(), pattern="\\.txt\\.gz$")
          
          # Only process actually selected GSMs
          selected_files <- character(0)
          file_to_gsm <- list()  # Keep track of which file belongs to which GSM
          
          for(gsm in input$selectedSamples) {
            pattern <- paste0("^", gsm, "_")
            matches <- grep(pattern, all_files, value=TRUE)
            selected_files <- c(selected_files, matches)
            file_to_gsm[matches] <- gsm
          }
          
          print(paste("Selected files:", paste(selected_files, collapse=", ")))
          
          if(length(selected_files) == 0) {
            stop("No matching files found for selected GSMs")
          }
          
          # Initialize list for all Seurat objects
          seurat_objects <- list()
          count = 0
          # Process each file
          for(file in selected_files) {
            count = count + 1
            incProgress(0.9/(length(selected_files)+1), detail = paste("Reading file",count,"of",length(selected_files)))
            # Read the file
            data <- Read_GEO_Delim(data_dir = selected_dir(), 
                                   file_suffix = file)
            
            if(length(data) > 0 && !is.null(data[[1]])) {
              # Create Seurat object for this sample
              gsm <- file_to_gsm[[file]]
              
              seurat <- CreateSeuratObject(counts = data[[1]], project = gsm)
              seurat$sample <- gsm
              
              # Add GEO metadata if available
              if (!is.null(geo_metadata())) {
                metadata <- geo_metadata()
                sample_meta <- metadata[metadata$geo_accession == gsm, ]
                
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
            final_seurat <- merge(seurat_objects[[1]], 
                                  y = seurat_objects[2:length(seurat_objects)],
                                  add.cell.ids = names(seurat_objects))
          } else {
            final_seurat <- seurat_objects[[1]]
          }
          
          # Store the object
          seurat_obj(final_seurat)
          
        }, error = function(e) {
          print(paste("Error in data processing:", e$message))
          print(traceback())
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