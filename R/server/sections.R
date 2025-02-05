# R/server/sections.R

setupSections <- function(output, seurat_data, metadata, processed_seurat, clustered_seurat, session) {
  # Metadata section with side-by-side layout
  output$metadataSection <- renderUI({
    req(metadata())
    ns <- session$ns  # Get namespace from session
    
    div(id = "metadata-section",
        h3(class = "section-header", "Sample Metadata"),
        fluidRow(
          # Left column for checkboxes
          column(2,
                 class = "sample-select-column",
                 div(class = "sample-select-buttons",
                     actionButton(ns("selectAll"), "Select All"),
                     actionButton(ns("deselectAll"), "Deselect All")
                 ),
                 div(class = "checkbox-container",
                     checkboxGroupInput(ns("selectedSamples"), 
                                        label = NULL,
                                        choices = setNames(metadata()$geo_accession, 
                                                           metadata()$geo_accession))
                 )
          ),
          # Right column for metadata table
          column(10,
                 renderDT({
                   datatable(metadata(),
                             options = list(
                               pageLength = 15,
                               scrollX = TRUE,
                               dom = 'tlip'
                             ),
                             rownames = FALSE,
                             selection = 'none',
                             class = 'cell-border stripe')
                 })
          )
        )
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