# R/server/sections.R

setupSections <- function(output, seurat_data, metadata, processed_seurat, clustered_seurat, session) {
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
          
          datatable(
            md,
            options = list(
              pageLength = 15,
              scrollX = TRUE,
              dom = 'tlip',
              columnDefs = list(
                list(
                  targets = 0,
                  render = JS("function(data, type, row, meta) {
                    return '<input type=\"checkbox\" class=\"sample-select\" checked data-gsm=\"' + row[1] + '\">';
                  }"),
                  orderable = FALSE,
                  width = '30px',
                  className = 'dt-center',
                  title = '<input type=\"checkbox\" id=\"select-all-samples\" checked>'
                )
              ),
              initComplete = JS("
                function(settings, json) {
                  var table = this.api();
                  
                  // Initialize selected samples
                  var initialSelected = [];
                  table.$('input.sample-select:checked').each(function() {
                    initialSelected.push($(this).data('gsm'));
                  });
                  Shiny.setInputValue('dataInput-selectedSamples', initialSelected);
                  
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
  
  # Rest of the sections remain unchanged
  output$qcSection <- renderUI({
    req(seurat_data)
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
}