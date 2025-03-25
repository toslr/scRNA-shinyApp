# R/server.R

buildServer <- function() {
  function(input, output, session) {
    print("Starting server initialization")
    
    # Initialize metadata module first
    metadata_module <- metadataServer("metadata")
    
    metadata_handler <- reactive({
      metadata_module
    })
    
    # Initialize data input module
    data_input <- dataInputServer("dataInput", metadata_module = metadata_module)
    
    # Extract Seurat object with proper reactive chain
    seurat_data <- reactive({
      req(data_input())
      data_input()
    })
    
    # Chain the reactive values through the modules
    processed_seurat <- qcServer("qc", seurat_data)
    clustered_seurat <- dimensionReductionServer("dimRed", processed_seurat)
    
    cluster_management <- clusterManagementServer("clusterManagement", clustered_seurat)
    
    de_module <- deAnalysisServer("de", clustered_seurat, cluster_management)
    
    # Create reactive values to track completion of each step
    steps_completed <- reactiveValues(
      data_input = FALSE,
      metadata = FALSE,
      qc = FALSE,
      dimred = FALSE,
      clustering = FALSE,
      de = FALSE
    )
    
    # Initialize save/load module
    loaded_analysis <- saveLoadServer("saveLoad", 
                                      seurat_data, 
                                      metadata_module, 
                                      processed_seurat, 
                                      clustered_seurat, 
                                      de_module,
                                      steps_completed,
                                      session)
    
    observe({
      data <- loaded_analysis()
      if (is.null(data)) return()
      
      # Update your data modules with the loaded data
      if (!is.null(data$seurat_data)) {
        data_input(data$seurat_data)  # Here data_input is the reactiveVal that seurat_data depends on
      }
    })
    
    # Render cluster management UI in sidebar only when clustering is done
    output$clusterManagementSection <- renderUI({
      req(clustered_seurat())
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      
      clusterManagementUI("clusterManagement")
    })
    
    # Setup observers
    setupObservers(steps_completed, seurat_data, metadata_module, processed_seurat, 
                   clustered_seurat, de_module)
    
    # Setup sections
    setupSections(input, output, seurat_data, metadata_handler, processed_seurat, 
                  clustered_seurat, session)
    
    # Setup navigation
    setupNavigation(output, steps_completed)
    

    print("Server initialization complete")
  }
}