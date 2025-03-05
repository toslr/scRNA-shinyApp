# R/server/sections.R

setupSections <- function(input, output, seurat_data, metadata_handler, processed_seurat, 
                          clustered_seurat, session) {
  
  output$metadataSection <- renderUI({
    # Call the reactive function to get metadata module
    print("Checking metadata handler")
    req(metadata_handler)
    metadata_module <- metadata_handler()
    
    req(metadata_module$getMetadata())
    
    div(id = "metadata-section",
        h3(class = "section-header", "Sample Metadata"),
        DTOutput(NS("metadata", "metadataTable"))
    )
  })
  
  # QC section
  output$qcSection <- renderUI({
    req(seurat_data())
    div(id = "qc-section",
        h3(class = "section-header", "Quality Control"),
        div(class = "section-content",
            qcUI("qc")
        )
    )
  })
  
  # Dimension reduction section
  output$dimredSection <- renderUI({
    req(processed_seurat())
    div(id = "dimred-section",
        h3(class = "section-header", "Dimension Reduction"),
        div(class = "section-content",
            dimensionReductionUI("dimRed")
        )
    )
  })
  
  # DE section
  output$deSection <- renderUI({
    req(clustered_seurat())
    req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
    req("umap" %in% names(clustered_seurat()@reductions))
    div(id = "de-section",
        h3(class = "section-header", "Differential Expression"),
        div(class = "section-content",
            deAnalysisUI("de")
        )
    )
  })
}