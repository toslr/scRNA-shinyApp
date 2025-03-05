# R/modules/save_load_module.R

saveLoadUI <- function(id) {
  ns <- NS(id)
  
  div(class = "save-load-buttons",
      style = "display: flex; gap: 10px;",
      actionButton(ns("showSave"), "Save Analysis", 
                   class = "btn-primary",
                   icon = icon("save")),
      actionButton(ns("showLoad"), "Load Analysis", 
                   class = "btn-primary",
                   icon = icon("folder-open"))
  )
}

saveLoadServer <- function(id, seurat_data, metadata_module, processed_seurat, 
                           clustered_seurat, de_module, steps_completed, session) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Create save directory if it doesn't exist
    save_dir <- "analysis_saves"
    if (!dir.exists(save_dir)) {
      dir.create(save_dir)
    }
    
    # Get list of saved analyses
    get_saved_analyses <- reactive({
      files <- list.files(save_dir, pattern = "\\.rds$", full.names = FALSE)
      files
    })
    
    # Show save dialog
    observeEvent(input$showSave, {
      showModal(modalDialog(
        title = "Save Analysis",
        textInput(ns("saveName"), "Analysis Name:", placeholder = "my_analysis"),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("saveAnalysis"), "Save", class = "btn-primary")
        )
      ))
    })
    
    # Show load dialog
    observeEvent(input$showLoad, {
      saved_files <- get_saved_analyses()
      
      showModal(modalDialog(
        title = "Load Analysis",
        if (length(saved_files) > 0) {
          selectInput(ns("loadFile"), "Select Analysis:", 
                      choices = setNames(saved_files, tools::file_path_sans_ext(saved_files)))
        } else {
          p("No saved analyses found")
        },
        footer = tagList(
          modalButton("Cancel"),
          if (length(saved_files) > 0) {
            actionButton(ns("loadAnalysis"), "Load", class = "btn-primary")
          }
        )
      ))
    })
    
    # Create a reactive to hold loaded data
    loaded_data <- reactiveVal(NULL)
    
    # Handle save
    observeEvent(input$saveAnalysis, {
      req(input$saveName)
      
      withProgress(message = 'Saving analysis...', value = 0, {
        tryCatch({
          print("Starting save process...")
          
          # Get data safely
          current_seurat <- NULL
          current_processed <- NULL
          current_clustered <- NULL
          
          # Use isolate to get values without creating dependencies
          isolate({
            # Get Seurat objects
            if (is.reactive(seurat_data)) {
              current_seurat <- seurat_data()
            }
            
            if (is.reactive(processed_seurat)) {
              current_processed <- processed_seurat()
            }
            
            if (is.reactive(clustered_seurat)) {
              current_clustered <- clustered_seurat()
            }
            
            # Get metadata
            current_metadata <- metadata_module$getMetadata()
            current_selected_samples <- metadata_module$selectedSamples()
            
            # Get DE data
            if (!is.null(de_module)) {
              current_de_results <- de_module$results()
              current_cluster_labels <- de_module$labels()
              current_active_clusters <- de_module$active()
            } else {
              current_de_results <- NULL
              current_cluster_labels <- NULL
              current_active_clusters <- NULL
            }
            
            # Get step completion status
            current_steps <- reactiveValuesToList(steps_completed)
          })
          
          # Create save object
          save_object <- list(
            timestamp = Sys.time(),
            seurat_data = current_seurat,
            processed_seurat = current_processed,
            clustered_seurat = current_clustered,
            metadata = current_metadata,
            selected_samples = current_selected_samples,
            de_results = current_de_results,
            cluster_labels = current_cluster_labels,
            active_clusters = current_active_clusters,
            steps_completed = current_steps
          )
          
          # Create save directory if it doesn't exist
          if (!dir.exists(save_dir)) {
            dir.create(save_dir, recursive = TRUE)
          }
          
          # Save to file
          filename <- file.path(save_dir, paste0(input$saveName, ".rds"))
          saveRDS(save_object, filename)
          
          removeModal()
          showNotification("Analysis saved successfully!", type = "message")
          
        }, error = function(e) {
          print(paste("Error in save process:", e$message))
          print(traceback())
          showNotification(paste("Error saving analysis:", e$message), type = "error")
        })
      })
    })
    
    # Handle load
    observeEvent(input$loadAnalysis, {
      req(input$loadFile)
      
      withProgress(message = 'Loading analysis...', value = 0, {
        tryCatch({
          # Load the file
          filename <- file.path(save_dir, input$loadFile)
          data <- readRDS(filename)
          
          # Store the data
          loaded_data(data)
          
          # Update completion steps first
          for (step in names(data$steps_completed)) {
            steps_completed[[step]] <- data$steps_completed[[step]]
          }
          
          # Update metadata
          metadata_module$updateMetadata(data$metadata)
          
          # Update DE module if it exists
          if (!is.null(de_module)) {
            if (!is.null(data$de_results)) {
              de_module$results(data$de_results)
            }
            if (!is.null(data$cluster_labels)) {
              de_module$labels(data$cluster_labels)
            }
            if (!is.null(data$active_clusters)) {
              de_module$active(data$active_clusters)
            }
          }
          
          removeModal()
          showNotification("Analysis loaded successfully!", type = "message")
          
        }, error = function(e) {
          print(paste("Error loading analysis:", e$message))
          print(traceback())
          showNotification(paste("Error loading analysis:", e$message), type = "error")
        })
      })
    })
    
    # Return the loaded data
    return(loaded_data)
  })
}