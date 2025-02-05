# R/server.R

buildServer <- function() {
  function(input, output, session) {
    # Get data from input module
    data_input <- dataInputServer("dataInput")
    
    # Extract Seurat object with proper reactive chain
    seurat_data <- reactive({
      req(data_input())
      input_data <- data_input()
      req(input_data$seurat)
      input_data$seurat
    })
    
    # Extract metadata if needed
    metadata <- reactive({
      req(data_input())
      input_data <- data_input()
      req(input_data$metadata)
      input_data$metadata
    })
    
    # Chain the reactive values through the modules
    processed_seurat <- qcServer("qc", seurat_data)
    clustered_seurat <- dimensionReductionServer("dimRed", processed_seurat)
    de_module <- deAnalysisServer("de", clustered_seurat)
    
    # Create reactive values to track completion of each step
    steps_completed <- reactiveValues(
      data_input = FALSE,
      metadata = FALSE,
      qc = FALSE,
      dimred = FALSE,
      clustering = FALSE,
      de = FALSE
    )
    
    # Setup different components using the modularized functions
    setupObservers(steps_completed, seurat_data, metadata, processed_seurat, 
                   clustered_seurat, de_module)
    setupSections(input,output, seurat_data, metadata, processed_seurat, 
                  clustered_seurat, session)  # Pass session here
    setupNavigation(output, steps_completed)
  }
}