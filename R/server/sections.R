# R/server/sections.R

setupSections <- function(output, seurat_data, metadata, processed_seurat, clustered_seurat, session) {
  # Metadata section
  output$metadataSection <- renderUI({
    req(metadata())
    div(id = "metadata-section",
        h3(class = "section-header", "Sample Metadata"),
        renderDT({
          # Add a selection column to the metadata
          md <- metadata()
          md$selected <- TRUE  # Default all to selected
          
          datatable(
            md,
            options = list(
              pageLength = 15,
              scrollX = TRUE,
              dom = 'tlip',
              columnDefs = list(
                list(
                  targets = 0,  # First column (selection)
                  render = JS("function(data, type, row, meta) {
                    return '<input type=\"checkbox\" ' + (data ? 'checked' : '') + '>';
                  }"),
                  orderable = FALSE
                )
              ),
              # Custom initialization to set up checkbox handling
              initComplete = JS("
                function(settings, json) {
                  var table = this.api();
                  // Handle checkbox changes
                  table.on('change', 'input[type=\"checkbox\"]', function() {
                    var $row = $(this).closest('tr');
                    var data = table.row($row).data();
                    var selectedGSMs = [];
                    table.$('input[type=\"checkbox\"]:checked').each(function() {
                      var rowData = table.row($(this).closest('tr')).data();
                      selectedGSMs.push(rowData[1]); // geo_accession is in second column
                    });
                    Shiny.setInputValue('selectedSamples', selectedGSMs);
                  });
                }
              ")
            ),
            selection = 'none',
            escape = FALSE,
            class = 'cell-border stripe'
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