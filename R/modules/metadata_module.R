#' @title Metadata Module UI
#' @description Creates the UI for the metadata module which allows users to input GEO
#'   Series IDs and fetch metadata for scRNA-seq analysis.
#' @param id The module ID
#' @return A Shiny UI element containing the metadata input interface
#' @export
metadataUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Intro text
    p("Welcome to SC Explorer! Please start by selecting your GEO dataset."),

    tabsetPanel(
      id = ns("metadataSource"),
      tabPanel(
        "GEO Metadata",
        div(style = "margin-top: 15px;",
            # GEO input section
            strong("Input GEO Series ID:"),
            div(style = "margin-top: 10px;",
                textInput(ns("geoID"), 
                          label = NULL,
                          placeholder = "GSExxxxxx")
            ),
            div(style = "margin-top: 5px;",
                actionButton(ns("fetchGEO"), "Fetch GEO Metadata")
            )
        )
      ),
      tabPanel(
        "Custom Metadata",
        div(style = "margin-top: 15px;",
            p("Upload a CSV/TSV file with your metadata. The file should contain:"),
            tags$ul(
              tags$li("One row per sample"),
              tags$li("A column named 'geo_accession' or 'sample_id' to match with data files"),
              tags$li("Additional columns with sample information")
            ),
            fileInput(ns("metadataFile"), "Choose Metadata File",
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        "text/tab-separated-values",
                        ".csv", ".tsv"
                      )),
            checkboxInput(ns("headerRow"), "File contains header row", TRUE),
            selectInput(ns("separator"), "Column separator:",
                        choices = c(Comma = ",", Tab = "\t", Semicolon = ";"),
                        selected = ","),
            actionButton(ns("parseMetadata"), "Parse Metadata File")
        )
      )
    ),
    
    # Status and preview area
    div(style = "margin-top: 10px;",
        verbatimTextOutput(ns("metadataStatus")),
        uiOutput(ns("metadataPreview"))
    )
  )
}

#' @title Metadata Module Server
#' @description Server logic for the metadata module which fetches and processes
#'   GEO metadata or parses custom metadata for scRNA-seq analysis.
#' @param id The module ID
#' @return A list of reactive expressions for accessing metadata and selected samples
#' @export
metadataServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize reactive values
    geo_metadata <- reactiveVal(NULL)
    selected_samples <- reactiveVal(NULL)
    metadata_status <- reactiveVal("")
    
    # GEO metadata fetching
    observeEvent(input$fetchGEO, {
      req(input$geoID)
      metadata_status("Fetching GEO metadata...")
      
      tryCatch({
        withProgress(message = 'Fetching GEO metadata...', value = 0, {
          # Fetch GEO data
          gset <- getGEO(input$geoID, GSEMatrix = TRUE)
          if (is(gset, "list")) gset <- gset[[1]]
          
          # Extract phenoData
          pheno_data <- pData(phenoData(gset))
          
          # Create basic metadata dataframe
          metadata <- data.frame(
            geo_accession = rownames(pheno_data),
            title = pheno_data$title,
            stringsAsFactors = FALSE
          )
          
          # Get all characteristics columns
          char_columns <- grep("characteristics_ch", colnames(pheno_data), value = TRUE)
          for (col in char_columns) {
            metadata[[col]] <- pheno_data[[col]]
          }
          
          geo_metadata(metadata)
          selected_samples(rownames(pheno_data))
          metadata_status("GEO metadata successfully fetched!")
        })
      }, error = function(e) {
        metadata_status(paste("Error fetching GEO data:", e$message))
      })
    })
    
    # Custom metadata parsing
    observeEvent(input$parseMetadata, {
      req(input$metadataFile)
      metadata_status("Parsing custom metadata file...")
      
      tryCatch({
        withProgress(message = 'Parsing metadata file...', value = 0, {
          # Read the file based on selected separator
          separator <- input$separator
          has_header <- input$headerRow
          
          file_path <- input$metadataFile$datapath
          
          # Determine if it's CSV or TSV based on extension
          if (endsWith(tolower(input$metadataFile$name), ".tsv")) {
            separator <- "\t"
          }
          
          # Read the metadata file
          metadata <- read.table(
            file_path,
            header = has_header,
            sep = separator,
            stringsAsFactors = FALSE,
            check.names = FALSE,
            quote = "\"",
            comment.char = "",
            fill = TRUE
          )
          
          # Check for required columns and rename if needed
          required_cols <- c("geo_accession", "sample_id")
          if (!any(required_cols %in% colnames(metadata))) {
            # If neither column exists, try to use the first column as sample_id
            if (ncol(metadata) > 0) {
              # Rename first column to sample_id
              names(metadata)[1] <- "sample_id"
              metadata_status("Warning: No 'geo_accession' or 'sample_id' column found. Using first column as sample IDs.")
            } else {
              stop("Metadata file must contain at least one column for sample IDs")
            }
          }
          
          # Standardize column name for sample ID
          if ("sample_id" %in% colnames(metadata) && !("geo_accession" %in% colnames(metadata))) {
            # Add geo_accession column as a copy of sample_id for compatibility
            metadata$geo_accession <- metadata$sample_id
          } else if ("geo_accession" %in% colnames(metadata) && !("sample_id" %in% colnames(metadata))) {
            # Add sample_id column as a copy of geo_accession
            metadata$sample_id <- metadata$geo_accession
          }
          
          # Add title column if missing (for consistency with GEO format)
          if (!("title" %in% colnames(metadata))) {
            metadata$title <- metadata$geo_accession
          }
          
          # Store the metadata and selected samples
          geo_metadata(metadata)
          selected_samples(metadata$geo_accession)
          metadata_status("Custom metadata successfully loaded!")
        })
      }, error = function(e) {
        metadata_status(paste("Error parsing metadata file:", e$message))
      })
    })
    
    # Status output
    output$metadataStatus <- renderText({
      metadata_status()
    })
    
    # Metadata preview
    output$metadataPreview <- renderUI({
      req(geo_metadata())
      
      # Create a preview of the metadata
      metadata <- geo_metadata()
      
      if (nrow(metadata) > 0) {
        # Show a small preview table
        tagList(
          h4("Metadata Preview:"),
          div(style = "max-height: 300px; overflow-y: auto;",
              renderDT(
                metadata[1:min(5, nrow(metadata)), 1:min(5, ncol(metadata))],
                options = list(dom = 't', scrollX = TRUE)
              )
          ),
          p(paste("Total samples:", nrow(metadata), "| Total columns:", ncol(metadata))),
          div(style = "margin-top: 10px;",
              actionButton(ns("showMetadataTable"), "View Full Metadata Table")
          )
        )
      } else {
        p("No metadata available to preview.")
      }
    })
    
    # Show full metadata table when button is clicked
    observeEvent(input$showMetadataTable, {
      req(geo_metadata())
      
      showModal(modalDialog(
        title = "Full Metadata Table",
        size = "l",
        easyClose = TRUE,
        
        DTOutput(ns("fullMetadataTable")),
        
        footer = modalButton("Close")
      ))
    })
    
    # Render full metadata table
    output$fullMetadataTable <- renderDT({
      req(geo_metadata())
      metadata <- geo_metadata()
      
      datatable(
        metadata,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = 'ltip',
          columnDefs = list(
            list(
              targets = 0,
              render = JS("
                function(data, type, row, meta) {
                  return '<input type=\"checkbox\" class=\"sample-select\" checked data-gsm=\"' + row[0] + '\">';
                }
              "),
              orderable = FALSE,
              width = '30px',
              className = 'dt-center',
              title = '<input type=\"checkbox\" id=\"select-all-samples\" checked>'
            )
          ),
          initComplete = JS("
            function(settings, json) {
              var table = this.api();
              
              // Handle select-all checkbox
              $('#select-all-samples').on('change', function() {
                var checked = this.checked;
                table.$('input.sample-select').prop('checked', checked);
                
                // Update selected samples
                var selectedGSMs = [];
                if (checked) {
                  table.$('input.sample-select').each(function() {
                    selectedGSMs.push($(this).data('gsm'));
                  });
                }
                Shiny.setInputValue('metadata-selected', selectedGSMs);
              });
              
              // Handle individual checkboxes
              table.on('change', 'input.sample-select', function() {
                var allChecked = table.$('input.sample-select:not(:checked)').length === 0;
                $('#select-all-samples').prop('checked', allChecked);
                
                // Update selected samples
                var selectedGSMs = [];
                table.$('input.sample-select:checked').each(function() {
                  selectedGSMs.push($(this).data('gsm'));
                });
                Shiny.setInputValue('metadata-selected', selectedGSMs);
              });
            }
          ")
        ),
        selection = 'none',
        escape = FALSE,
        class = 'cell-border stripe',
        rownames = FALSE
      )
    })
    
    # Update selected_samples when input changes
    observeEvent(input$selected, {
      selected_samples(input$selected)
    })
    
    # Set selected samples (for restoring state)
    setSelectedSamples <- function(samples) {
      if (!is.null(samples) && length(samples) > 0) {
        selected_samples(samples)
        
        # Update checkboxes if the DT is already rendered
        session$sendCustomMessage(
          type = "updateSampleCheckboxes",
          message = list(samples = samples)
        )
      }
    }
    
    # Return reactive expressions
    return(list(
      getMetadata = reactive({ geo_metadata() }),
      selectedSamples = reactive({ selected_samples() }),
      updateMetadata = function(new_metadata) {
        geo_metadata(new_metadata)
      },
      setSelectedSamples = setSelectedSamples
    ))
  })
}