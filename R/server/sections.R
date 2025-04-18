# R/server/sections.R

#' @title Setup Main Content Sections
#' @description Sets up the main content sections of the single-cell RNA-seq analysis application.
#'   This function configures the dynamic rendering of the metadata, quality control, dimension 
#'   reduction, and differential expression sections based on data availability and analysis state.
#' @param input Shiny input object
#' @param output Shiny output object for rendering UI elements
#' @param seurat_data Reactive expression containing the initial Seurat object
#' @param metadata_handler Reactive expression providing access to the metadata module
#' @param processed_seurat Reactive expression containing the QC-processed Seurat object
#' @param clustered_seurat Reactive expression containing the clustered Seurat object
#' @param session The current Shiny session
#' @return None (used for its side effects of setting up UI sections)
#' @export
setupSections <- function(input, output, seurat_data, metadata_handler, processed_seurat, 
                          clustered_seurat, session) {
  
  # Metadata section
  output$metadataSection <- renderUI({
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
    
    # Check if we have clustering information in the Seurat object
    valid_clustering <- FALSE
    valid_umap_exists <- FALSE
    
    tryCatch({
      # Check if seurat_clusters column exists in metadata
      valid_clustering <- "seurat_clusters" %in% colnames(clustered_seurat()@meta.data)
      
      # Check if umap reductions exist
      valid_umap_exists <- any(c("umap2d", "umap3d", "umap") %in% names(clustered_seurat()@reductions))
    }, error = function(e) {
      print(paste("Error checking clustered_seurat:", e$message))
    })
    
    req(valid_clustering)
    req(valid_umap_exists)
    
    div(id = "de-section",
        h3(class = "section-header", "Differential Expression"),
        div(class = "section-content",
            deAnalysisUI("de")
        )
    )
  })
}