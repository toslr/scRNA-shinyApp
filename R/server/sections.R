# R/sections.R

setupSections <- function(output, seurat_data, processed_seurat, clustered_seurat) {
  # Conditional section renders
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