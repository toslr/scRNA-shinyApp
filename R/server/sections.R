# R/sections.R

setupSections <- function(output, seurat_data, processed_seurat, clustered_seurat) {
  # Conditional section renders
  output$qcSection <- renderUI({
    req(seurat_data())
    div(id = "qc-section",
        h3(class = "section-header", "Quality Control"),
        qcUI("qc")
    )
  })
  
  output$dimredSection <- renderUI({
    req(processed_seurat())
    div(id = "dimred-section",
        h3(class = "section-header", "Dimension Reduction"),
        dimensionReductionUI("dimRed")
    )
  })
  
  output$deSection <- renderUI({
    req(clustered_seurat())
    req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
    div(id = "de-section",
        h3(class = "section-header", "Differential Expression"),
        deAnalysisUI("de")
    )
  })
}