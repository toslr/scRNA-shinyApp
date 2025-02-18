# R/modules/save_load_module.R

#' Helper function to safely convert reactive values to list
safe_reactive_to_list <- function(rv) {
  if (is.null(rv)) return(NULL)
  if (is.function(rv)) return(NULL)  # Skip functions
  if (inherits(rv, "reactivevalues")) {
    return(reactiveValuesToList(rv))
  }
  if (is.reactive(rv)) {
    tryCatch({
      value <- rv()
      if (is.function(value)) return(NULL)
      return(value)
    }, error = function(e) {
      return(NULL)
    })
  }
  return(rv)
}

saveLoadUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(class = "save-load-controls",
        wellPanel(
          h4("Save Analysis"),
          textInput(ns("saveName"), "Analysis Name:", placeholder = "my_analysis"),
          selectInput(ns("saveStage"), "Save Stage:", 
                      choices = c(
                        "Post-QC" = "qc",
                        "Post-Clustering" = "clustering",
                        "Full Analysis" = "complete"
                      )),
          actionButton(ns("saveAnalysis"), "Save Current Analysis", 
                       class = "btn-primary")
        ),
        wellPanel(
          h4("Load Analysis"),
          uiOutput(ns("loadControls"))
        )
    )
  )
}

saveLoadServer <- function(id, seurat_data, metadata_module, processed_seurat, 
                           clustered_seurat, de_module, steps_completed) {
  moduleServer(id, function(input, output, session) {
    ns <- NS(id)
    
    # Create save directory if it doesn't exist
    save_dir <- "analysis_saves"
    if (!dir.exists(save_dir)) {
      dir.create(save_dir)
    }
    
    # Function to safely save an analysis stage
    save_stage <- function(name, stage) {
      tryCatch({
        # Create stage-specific directory
        stage_dir <- file.path(save_dir, paste0(name, "_", stage))
        if (!dir.exists(stage_dir)) {
          dir.create(stage_dir, recursive = TRUE)
        }
        
        # Get current values safely
        current_seurat <- safe_reactive_to_list(seurat_data)
        current_metadata <- safe_reactive_to_list(metadata_module$getMetadata)
        current_processed <- safe_reactive_to_list(processed_seurat)
        current_clustered <- safe_reactive_to_list(clustered_seurat)
        current_steps <- safe_reactive_to_list(steps_completed)
        
        # Create appropriate save object based on stage
        save_object <- list(
          stage = stage,
          timestamp = Sys.time(),
          r_version = R.version.string,
          seurat_version = packageVersion("Seurat")
        )
        
        # Add stage-specific data
        if (stage == "qc") {
          save_object$data <- list(
            seurat = current_seurat,
            metadata = current_metadata,
            processed = current_processed,
            steps = current_steps[c("data_input", "metadata", "qc")]
          )
        } else if (stage == "clustering") {
          save_object$data <- list(
            seurat = current_seurat,
            metadata = current_metadata,
            processed = current_processed,
            clustered = current_clustered,
            steps = current_steps[c("data_input", "metadata", "qc", "dimred", "clustering")]
          )
        } else if (stage == "complete") {
          save_object$data <- list(
            seurat = current_seurat,
            metadata = current_metadata,
            processed = current_processed,
            clustered = current_clustered,
            de_results = safe_reactive_to_list(de_module$results),
            cluster_labels = safe_reactive_to_list(de_module$labels),
            active_clusters = safe_reactive_to_list(de_module$active),
            steps = current_steps
          )
        }
        
        # Save the object
        saveRDS(save_object, file.path(stage_dir, "analysis_data.rds"))
        TRUE
        
      }, error = function(e) {
        warning(paste("Error saving analysis:", e$message))
        FALSE
      })
    }
    
    # Function to load an analysis
    load_stage <- function(file_path) {
      tryCatch({
        # Load the save object
        save_object <- readRDS(file_path)
        
        # Version check
        if (save_object$seurat_version != packageVersion("Seurat")) {
          warning("Warning: Different Seurat version detected")
        }
        
        # Restore data based on stage
        if (!is.null(save_object$data$seurat)) {
          seurat_data(save_object$data$seurat)
        }
        if (!is.null(save_object$data$metadata)) {
          metadata_module$updateMetadata(save_object$data$metadata)
        }
        if (!is.null(save_object$data$processed)) {
          processed_seurat(save_object$data$processed)
        }
        
        if (save_object$stage %in% c("clustering", "complete")) {
          if (!is.null(save_object$data$clustered)) {
            clustered_seurat(save_object$data$clustered)
          }
        }
        
        if (save_object$stage == "complete") {
          if (!is.null(save_object$data$de_results)) {
            de_module$results(save_object$data$de_results)
          }
          if (!is.null(save_object$data$cluster_labels)) {
            de_module$labels(save_object$data$cluster_labels)
          }
          if (!is.null(save_object$data$active_clusters)) {
            de_module$active(save_object$data$active_clusters)
          }
        }
        
        # Restore steps
        if (!is.null(save_object$data$steps)) {
          for (step in names(save_object$data$steps)) {
            steps_completed[[step]] <- save_object$data$steps[[step]]
          }
        }
        
        TRUE
        
      }, error = function(e) {
        warning(paste("Error loading analysis:", e$message))
        FALSE
      })
    }
    
    # Get list of saved analyses
    get_saved_analyses <- reactive({
      all_saves <- list.dirs(save_dir, full.names = FALSE, recursive = FALSE)
      saves_with_data <- all_saves[sapply(all_saves, function(x) {
        file.exists(file.path(save_dir, x, "analysis_data.rds"))
      })]
      
      # Extract metadata for each save
      lapply(saves_with_data, function(save_name) {
        save_data <- readRDS(file.path(save_dir, save_name, "analysis_data.rds"))
        list(
          name = sub("_[^_]+$", "", save_name),  # Remove stage suffix
          stage = save_data$stage,
          path = file.path(save_dir, save_name, "analysis_data.rds"),
          timestamp = save_data$timestamp
        )
      })
    })
    
    # Render load controls
    output$loadControls <- renderUI({
      saved_analyses <- get_saved_analyses()
      if (length(saved_analyses) == 0) {
        return(p("No saved analyses found"))
      }
      
      # Create choices with informative labels
      choices <- sapply(saved_analyses, function(x) {
        sprintf("%s (%s - %s)", 
                x$name, 
                x$stage, 
                format(x$timestamp, "%Y-%m-%d %H:%M"))
      })
      names(choices) <- sapply(saved_analyses, function(x) x$path)
      
      tagList(
        selectInput(ns("loadFile"), "Select Analysis:", choices = choices),
        actionButton(ns("loadAnalysis"), "Load Selected Analysis",
                     class = "btn-primary")
      )
    })
    
    # Save analysis
    observeEvent(input$saveAnalysis, {
      req(input$saveName, input$saveStage)
      
      withProgress(message = 'Saving analysis...', value = 0, {
        if (save_stage(input$saveName, input$saveStage)) {
          showNotification("Analysis saved successfully!", type = "message")
        } else {
          showNotification("Error saving analysis", type = "error")
        }
      })
    })
    
    # Load analysis
    observeEvent(input$loadAnalysis, {
      req(input$loadFile)
      
      withProgress(message = 'Loading analysis...', value = 0, {
        if (load_stage(input$loadFile)) {
          showNotification("Analysis loaded successfully!", type = "message")
        } else {
          showNotification("Error loading analysis", type = "error")
        }
      })
    })
  })
}