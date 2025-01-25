# R/server.R

buildServer <- function() {
  function(input, output, session) {
    # Chain the reactive values through the modules
    seurat_data <- dataInputServer("dataInput")
    processed_seurat <- qcServer("qc", seurat_data)
    clustered_seurat <- dimensionReductionServer("dimRed", processed_seurat)
    de_module <- deAnalysisServer("de", clustered_seurat)
    
    # Create reactive values to track completion of each step
    steps_completed <- reactiveValues(
      data_input = FALSE,
      qc = FALSE,
      dimred = FALSE,
      clustering = FALSE,
      de = FALSE
    )
    
    # Setup different components
    setupObservers(steps_completed, seurat_data, processed_seurat, 
                   clustered_seurat, de_module)
    setupSections(output, seurat_data, processed_seurat, clustered_seurat)
    setupNavigation(output, steps_completed)
  }
}