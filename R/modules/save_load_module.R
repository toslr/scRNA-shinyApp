# R/modules/save_load_module.R

saveLoadUI <- function(id) {
  ns <- NS(id)
  
  div(class = "save-load-buttons",
      style = "display: flex; flex-direction: column; gap: 10px;",
      div(style = "display: flex; gap: 10px;",
          actionButton(ns("showSave"), "Save Analysis", 
                       class = "btn-primary",
                       icon = icon("save")),
          actionButton(ns("showLoad"), "Load Analysis", 
                       class = "btn-primary",
                       icon = icon("folder-open"))
      ),
      # Add a checkbox to toggle between default and custom location
      checkboxInput(ns("useCustomLocation"), "Use custom location for saving/loading analyses", value = FALSE),
      # Container for directory path display (initially hidden)
      conditionalPanel(
        condition = sprintf("input['%s'] == true", ns("useCustomLocation")),
        div(style = "margin-top: 5px;",
            shinyDirButton(ns("selectDir"), "Select Directory", "Choose directory"),
            verbatimTextOutput(ns("selectedDirPath"))
        )
      )
  )
}

saveLoadServer <- function(id, seurat_data, metadata_module, processed_seurat, 
                           clustered_seurat, de_module, steps_completed, session,
                           sample_management = NULL, condition_management = NULL, cluster_management = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Define reactive values for directory paths
    values <- reactiveValues(
      save_dir = "analysis_saves",
      custom_save_dir = NULL
    )
    
    # Initialize default directory - use observe() to properly access reactive values
    observe({
      default_dir <- values$save_dir
      if (!dir.exists(default_dir)) {
        dir.create(default_dir)
      }
    })
    
    # Define roots for directory selection
    volumes <- c(
      Home = Sys.getenv("HOME"),
      Desktop = file.path(Sys.getenv("HOME"), "Desktop"),
      Documents = file.path(Sys.getenv("HOME"), "Documents")
    )
    
    # Add platform-specific drives for Windows
    if (.Platform$OS.type == "windows") {
      for (drive in c("C:", "D:", "E:")) {
        if (dir.exists(drive)) {
          volumes[drive] <- drive
        }
      }
    }
    
    # Setup directory selection
    shinyDirChoose(input, 'selectDir', roots = volumes, session = session)
    
    # Update selected directory path
    observe({
      req(input$selectDir)
      values$custom_save_dir <- parseDirPath(volumes, input$selectDir)
    })
    
    # Display selected directory path
    output$selectedDirPath <- renderText({
      if (!is.null(values$custom_save_dir)) {
        return(values$custom_save_dir)
      } else {
        return("No directory selected")
      }
    })
    
    # Get the active save directory
    get_active_save_dir <- reactive({
      if (input$useCustomLocation && !is.null(values$custom_save_dir)) {
        return(values$custom_save_dir)
      } else {
        return(values$save_dir)
      }
    })
    
    # Get list of saved analyses
    get_saved_analyses <- reactive({
      save_dir <- get_active_save_dir()
      if (!dir.exists(save_dir)) {
        return(character(0))
      }
      files <- list.files(save_dir, pattern = "\\.rds$", full.names = FALSE)
      files
    })
    
    # Show save dialog
    observeEvent(input$showSave, {
      save_dir <- get_active_save_dir()
      
      # Create directory if it doesn't exist
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
      }
      
      showModal(modalDialog(
        title = "Save Analysis",
        textInput(ns("saveName"), "Analysis Name:", placeholder = "my_analysis"),
        if (input$useCustomLocation) {
          tags$div(
            tags$p(paste("Analysis will be saved to:", save_dir)),
            tags$p("You can change this location using the 'Select Directory' button.")
          )
        },
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("saveAnalysis"), "Save", class = "btn-primary")
        )
      ))
    })
    
    # Show load dialog
    observeEvent(input$showLoad, {
      save_dir <- get_active_save_dir()
      saved_files <- get_saved_analyses()
      
      showModal(modalDialog(
        title = "Load Analysis",
        if (length(saved_files) > 0) {
          selectInput(ns("loadFile"), "Select Analysis:", 
                      choices = setNames(saved_files, tools::file_path_sans_ext(saved_files)))
        } else {
          if (dir.exists(save_dir)) {
            p(paste0("No saved analyses found in ", save_dir))
          } else {
            p(paste0("Directory does not exist: ", save_dir))
          }
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
      save_dir <- get_active_save_dir()
      
      withProgress(message = 'Saving analysis...', value = 0, {
        tryCatch({
          print(paste("Starting save process to directory:", save_dir))
          
          # Get data safely
          current_seurat <- NULL
          current_processed <- NULL
          current_clustered <- NULL
          
          # Use isolate to get values without creating dependencies
          isolate({
            # Get Seurat objects
            current_seurat <- NULL
            tryCatch({
              if (is.reactive(seurat_data)) {
                current_seurat <- seurat_data()
              }
            }, error = function(e) {
              print(paste("Error getting seurat_data:", e$message))
            })
            
            # Process data with error handling
            tryCatch({
              if (is.list(processed_seurat) && is.function(processed_seurat$data)) {
                current_processed <- processed_seurat$data()
              } else if (is.reactive(processed_seurat)) {
                current_processed <- processed_seurat()
              }
            }, error = function(e) {
              print(paste("Error getting processed_seurat:", e$message))
            })
            
            # Get clustered data with error handling
            tryCatch({
              if (is.list(clustered_seurat) && is.function(clustered_seurat$data)) {
                current_clustered <- clustered_seurat$data()
              } else if (is.reactive(clustered_seurat)) {
                current_clustered <- clustered_seurat()
              }
            }, error = function(e) {
              print(paste("Error getting clustered_seurat:", e$message))
            })
            
            # Get metadata
            current_metadata <- NULL
            current_selected_samples <- NULL
            
            tryCatch({
              if (is.list(metadata_module)) {
                if (is.function(metadata_module$getMetadata)) {
                  current_metadata <- metadata_module$getMetadata()
                }
                if (is.function(metadata_module$selectedSamples)) {
                  current_selected_samples <- metadata_module$selectedSamples()
                }
              }
            }, error = function(e) {
              print(paste("Error getting metadata:", e$message))
            })
            
            # Get DE data
            current_de_results <- NULL
            current_de_analysis_type <- NULL
            current_heatmap_data <- NULL
            current_general_heatmap_genes <- NULL
            
            tryCatch({
              if (!is.null(de_module) && is.list(de_module)) {
                if (is.function(de_module$results)) {
                  current_de_results <- de_module$results()
                }
                # Safely check if function exists before calling
                if ("getAnalysisState" %in% names(de_module) && is.function(de_module$getAnalysisState)) {
                  current_de_analysis_type <- de_module$getAnalysisState()
                }
                if ("getHeatmapData" %in% names(de_module) && is.function(de_module$getHeatmapData)) {
                  current_heatmap_data <- de_module$getHeatmapData()
                }
                if ("getGeneralHeatmapGenes" %in% names(de_module) && is.function(de_module$getGeneralHeatmapGenes)) {
                  current_general_heatmap_genes <- de_module$getGeneralHeatmapGenes()
                }
              }
            }, error = function(e) {
              print(paste("Error getting DE data:", e$message))
            })
            
            # Get management module states
            sample_management_state <- NULL
            condition_management_state <- NULL
            cluster_management_state <- NULL
            
            tryCatch({
              if (!is.null(sample_management) && is.list(sample_management) && 
                  "getFullState" %in% names(sample_management) && is.function(sample_management$getFullState)) {
                sample_management_state <- sample_management$getFullState()
                print(paste("Sample management state successfully captured with", 
                            length(sample_management_state$all_samples), "samples"))
              }
            }, error = function(e) {
              print(paste("Error getting sample management state:", e$message))
            })
            
            tryCatch({
              if (!is.null(condition_management) && is.list(condition_management) && 
                  "getFullState" %in% names(condition_management) && is.function(condition_management$getFullState)) {
                condition_management_state <- condition_management$getFullState()
                print(paste("Condition management state successfully captured with column:", 
                            condition_management_state$condition_column))
              }
            }, error = function(e) {
              print(paste("Error getting condition management state:", e$message))
            })
            
            tryCatch({
              if (!is.null(cluster_management) && is.list(cluster_management) && 
                  "getFullState" %in% names(cluster_management) && is.function(cluster_management$getFullState)) {
                cluster_management_state <- cluster_management$getFullState()
                print(paste("Cluster management state successfully captured with", 
                            length(cluster_management_state$all_clusters), "clusters"))
                
                # Log active clusters
                active_clusters <- names(cluster_management_state$active_clusters[
                  unlist(cluster_management_state$active_clusters) == TRUE])
                print(paste("Active clusters:", paste(active_clusters, collapse = ", ")))
              }
            }, error = function(e) {
              print(paste("Error getting cluster management state:", e$message))
            })
            
            # Get step completion status
            current_steps <- reactiveValuesToList(steps_completed)
            
            # Capture QC parameter values from the global input
            qc_params <- NULL
            tryCatch({
              qc_params <- list(
                minFeature = input$`qc-minFeature`,
                maxFeature = input$`qc-maxFeature`,
                maxMT = input$`qc-maxMT`,
                qc_processed = !is.null(current_processed)
              )
            }, error = function(e) {
              print(paste("Error capturing QC parameters:", e$message))
            })
            
            # Capture PCA/UMAP parameters
            pca_params <- NULL
            tryCatch({
              dimred_complete <- FALSE
              if (!is.null(current_processed) && "pca" %in% names(current_processed@reductions)) {
                dimred_complete <- TRUE
              }
              
              pca_params <- list(
                nDims = input$`dimRed-nDims`,
                dims_confirmed = dimred_complete
              )
            }, error = function(e) {
              print(paste("Error capturing PCA parameters:", e$message))
            })
            
            # Capture clustering parameters
            clustering_params <- NULL
            tryCatch({
              clustering_params <- list(
                resolution = input$`dimRed-resolution`,
                clustering_done = !is.null(current_clustered) && 
                  "seurat_clusters" %in% colnames(current_clustered@meta.data)
              )
            }, error = function(e) {
              print(paste("Error capturing clustering parameters:", e$message))
            })
            
            # Capture DE parameters
            de_params <- NULL
            tryCatch({
              de_analysis_done <- !is.null(current_de_results) && length(current_de_results) > 0
              # For dataframes, check rows
              if (is.data.frame(current_de_results)) {
                de_analysis_done <- nrow(current_de_results) > 0
              }
              
              de_params <- list(
                target_cluster_all = input$`de-targetClusterAll`,
                target_cluster1 = input$`de-targetCluster1`,
                target_cluster2 = input$`de-targetCluster2`,
                genes_per_cluster = input$`de-genesPerCluster`,
                de_analysis_done = de_analysis_done,
                analysis_type = current_de_analysis_type
              )
            }, error = function(e) {
              print(paste("Error capturing DE parameters:", e$message))
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
              de_analysis_type = current_de_analysis_type,
              de_heatmap_data = current_heatmap_data,
              de_general_heatmap_genes = current_general_heatmap_genes,
              sample_management_state = sample_management_state,
              condition_management_state = condition_management_state,
              cluster_management_state = cluster_management_state,
              steps_completed = current_steps,
              ui_state = list(
                qc_params = qc_params,
                pca_params = pca_params,
                clustering_params = clustering_params,
                de_params = de_params
              )
            )
            
            # Create save directory if it doesn't exist
            if (!dir.exists(save_dir)) {
              dir.create(save_dir, recursive = TRUE)
            }
            
            # Save to file
            filename <- file.path(save_dir, paste0(input$saveName, ".rds"))
            saveRDS(save_object, filename)
            
            removeModal()
            showNotification(paste0("Analysis saved successfully to ", filename), type = "message")
          })
          
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
      save_dir <- get_active_save_dir()
      
      withProgress(message = 'Loading analysis...', value = 0, {
        tryCatch({
          # Load the file
          filename <- file.path(save_dir, input$loadFile)
          print(paste("Loading analysis from:", filename))
          
          if (!file.exists(filename)) {
            showNotification(paste("File not found:", filename), type = "error")
            return(NULL)
          }
          
          data <- readRDS(filename)
          
          # Store the data for other modules to access
          loaded_data(data)
          
          # Update completion steps first
          incProgress(0.1, detail = "Restoring analysis state")
          for (step in names(data$steps_completed)) {
            steps_completed[[step]] <- data$steps_completed[[step]]
          }
          
          # Update metadata
          incProgress(0.2, detail = "Restoring metadata")
          tryCatch({
            if (is.list(metadata_module) && is.function(metadata_module$updateMetadata)) {
              metadata_module$updateMetadata(data$metadata)
            }
          }, error = function(e) {
            print(paste("Error updating metadata:", e$message))
          })
          
          # Update DE module if it exists
          incProgress(0.3, detail = "Restoring analysis results")
          tryCatch({
            if (!is.null(de_module) && is.list(de_module)) {
              if (!is.null(data$de_results) && is.function(de_module$results)) {
                de_module$results(data$de_results)
              }
              if (!is.null(data$cluster_labels) && is.function(de_module$labels)) {
                de_module$labels(data$cluster_labels)
              }
              if (!is.null(data$active_clusters) && is.function(de_module$active)) {
                de_module$active(data$active_clusters)
              }
            }
          }, error = function(e) {
            print(paste("Error updating DE module:", e$message))
          })
          
          removeModal()
          showNotification(paste0("Analysis loaded successfully from ", filename), type = "message")
          
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