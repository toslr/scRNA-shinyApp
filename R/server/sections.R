# R/server/sections.R

setupSections <- function(input, output, seurat_data, metadata, processed_seurat, 
                          clustered_seurat, session) {  # Add input and session parameters
  # Keep track of current checkbox states
  selected_states <- reactiveVal()
  
  # Metadata section
  output$metadataSection <- renderUI({
    req(metadata())
    div(id = "metadata-section",
        h3(class = "section-header", "Sample Metadata"),
        renderDT({
          # Add empty column for checkboxes at the beginning
          md <- cbind(
            Select = rep("", nrow(metadata())),
            metadata()
          )
          
          # Get current checkbox states or initialize if NULL
          current_states <- selected_states()
          if(is.null(current_states)) {
            current_states <- rep(TRUE, nrow(md))
            selected_states(current_states)
          }
          
          datatable(
            md,
            options = list(
              pageLength = 15,
              scrollX = TRUE,
              dom = 'tlip',
              columnDefs = list(
                list(
                  targets = 0,
                  render = JS(sprintf("
                    function(data, type, row, meta) {
                      return '<input type=\"checkbox\" class=\"sample-select\" ' + 
                             (meta.row in %s ? 'checked ' : '') +
                             'data-gsm=\"' + row[1] + '\">';
                    }",
                                      jsonlite::toJSON(which(current_states) - 1)
                  )),
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
                    var selectedGSMs = [];
                    
                    table.$('input.sample-select').each(function() {
                      $(this).prop('checked', checked);
                      if (checked) {
                        selectedGSMs.push($(this).data('gsm'));
                      }
                    });
                    
                    Shiny.setInputValue('dataInput-selectedSamples', selectedGSMs);
                  });
                  
                  // Handle individual checkboxes
                  table.on('change', 'input.sample-select', function() {
                    var selectedGSMs = [];
                    table.$('input.sample-select:checked').each(function() {
                      selectedGSMs.push($(this).data('gsm'));
                    });
                    Shiny.setInputValue('dataInput-selectedSamples', selectedGSMs);
                    
                    // Update header checkbox
                    var allChecked = table.$('input.sample-select:not(:checked)').length === 0;
                    $('#select-all-samples').prop('checked', allChecked);
                  });
                }
              ")
            ),
            selection = 'none',
            escape = FALSE,
            class = 'cell-border stripe',
            colnames = c("", names(metadata())),
            rownames = FALSE
          )
        })
    )
  })
  
  # Other sections...
  output$qcSection <- renderUI({
    req(seurat_data())
    seurat_obj <- seurat_data()
    req(seurat_obj)
    div(id = "qc-section",
        h3(class = "section-header", "Quality Control"),
        qcUI("qc")
    )
  })
  
  output$dimredSection <- renderUI({
    req(processed_seurat())
    seurat_obj <- processed_seurat()
    req(seurat_obj)
    div(id = "dimred-section",
        h3(class = "section-header", "Dimension Reduction"),
        dimensionReductionUI("dimRed")
    )
  })
  
  output$deSection <- renderUI({
    req(clustered_seurat())
    seurat_obj <- clustered_seurat()
    req(seurat_obj)
    req("seurat_clusters" %in% colnames(seurat_obj@meta.data))
    div(id = "de-section",
        h3(class = "section-header", "Differential Expression"),
        deAnalysisUI("de")
    )
  })
  
  # Update checkbox states when they change
  observeEvent(input$`dataInput-selectedSamples`, {
    if(!is.null(metadata())) {
      current_states <- metadata()$geo_accession %in% input$`dataInput-selectedSamples`
      selected_states(current_states)
    }
  })
}