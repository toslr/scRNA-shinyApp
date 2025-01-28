# R/sections.R

setupSections <- function(output, seurat_data, metadata, processed_seurat, clustered_seurat) {
  # Conditional section renders
  
  output$metadataSection <- renderUI({
    req(metadata())
    div(id = "metadata-section",
        h3(class = "section-header", "Sample Metadata"),
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
  })
  
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