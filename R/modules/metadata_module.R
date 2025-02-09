# R/modules/metadata_module.R

metadataUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Intro text
    p("Welcome to scRNA Analysis! Please start by selecting your GEO dataset."),
    # GEO input section
    strong("Input GEO Series ID:"),
    div(style = "margin-top: 10px;",
        textInput(ns("geoID"), 
                  label = NULL,
                  placeholder = "GSExxxxxx")
    ),
    div(style = "margin-top: 5px;",
        actionButton(ns("fetchGEO"), "Fetch GEO Metadata")
    ),
    div(style = "margin-top: 5px;",
        verbatimTextOutput(ns("geoStatus"))
    ),
  )
}

metadataServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    print("Initializing metadata server")
    ns <- session$ns
    
    # Initialize reactive values
    geo_metadata <- reactiveVal(NULL)
    selected_samples <- reactiveVal(NULL)
    
    # GEO metadata fetching
    observeEvent(input$fetchGEO, {
      print("Fetch GEO button clicked")
      req(input$geoID)
      tryCatch({
        withProgress(message = 'Fetching GEO metadata...', value = 0, {
          # Fetch GEO data
          gset <- getGEO(input$geoID, GSEMatrix = TRUE)
          if (is(gset, "list")) gset <- gset[[1]]
          print("GEO data fetched successfully")
          
          # Extract phenoData
          pheno_data <- pData(phenoData(gset))
          print("Phenotype data extracted")
          
          # Create basic metadata dataframe
          metadata <- data.frame(
            geo_accession = rownames(pheno_data),
            title = pheno_data$title,
            stringsAsFactors = FALSE
          )
          
          # Get all characteristics columns
          char_columns <- grep("characteristics_ch", colnames(pheno_data), value = TRUE)
          
          # Add each characteristics column
          for (col in char_columns) {
            metadata[[col]] <- pheno_data[[col]]
          }
          
          print("Metadata dataframe created with dimensions:")
          print(dim(metadata))
          
          # Store metadata
          geo_metadata(metadata)
          print("Metadata stored in reactive value")
        })
        
        output$geoStatus <- renderText({
          "GEO metadata successfully fetched!"
        })
        
      }, error = function(e) {
        print(paste("Error in GEO fetching:", e$message))
        output$geoStatus <- renderText({
          paste("Error fetching GEO data:", e$message)
        })
      })
    })
    
    # Render metadata table
    output$metadataTable <- renderDT({
      print("Attempting to render metadata table...")
      req(geo_metadata())
      md <- geo_metadata()
      print("Metadata available for table:")
      print(dim(md))
      
      # Add empty column for checkboxes
      md <- cbind(
        Select = rep("", nrow(md)),
        md
      )
      
      datatable(
        md,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = 'tlip',
          columnDefs = list(
            list(
              targets = 0,
              render = JS("
                function(data, type, row, meta) {
                  return '<input type=\"checkbox\" class=\"sample-select\" checked data-gsm=\"' + row[1] + '\">';
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
        colnames = c("", names(md[-1])),
        rownames = FALSE
      )
    })
    
    # Update selected_samples when input changes
    observeEvent(input$selected, {
      selected_samples(input$selected)
    })
    
    print("Metadata server initialization complete")
    
    # Return reactive expressions
    list(
      getMetadata = reactive({ 
        print("Getting metadata from module")
        geo_metadata() 
      }),
      selectedSamples = reactive({ selected_samples() })
    )
  })
}