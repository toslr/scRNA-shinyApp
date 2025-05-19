#' @title Differential Expression Analysis Module UI
#' @description Creates the UI for the differential expression analysis module, allowing users to 
#'   compare gene expression between clusters and visualize results.
#' @param id The module ID
#' @return A Shiny UI element containing the DE analysis interface
#' @export
deAnalysisUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Analysis UI
    uiOutput(ns("analysisUI"))
  )
}

#' @title Differential Expression Analysis Module Server
#' @description Server logic for the DE analysis module that performs differential expression 
#'   analysis between clusters and generates visualizations.
#' @param id The module ID
#' @param clustered_seurat Reactive expression containing the clustered Seurat object
#' @param cluster_management Cluster management module instance for accessing active clusters
#' @param sample_management Optional sample management module for filtering by samples
#' @param condition_management Optional condition management module for filtering by conditions
#' @return A list of reactive expressions containing DE results and status
#' @export
deAnalysisServer <- function(id, clustered_seurat, cluster_management, sample_management = NULL, condition_management = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize state management
    state <- setupReactiveState()
    
    # Track changes to cluster labels 
    observe({
      # This will only trigger when labels are actually saved via the button
      req(cluster_management)
      
      # Use the last_update as a trigger
      if (!is.null(cluster_management$last_update)) {
        trigger_time <- cluster_management$last_update()
        if (!is.null(trigger_time)) {
          # When cluster labels change, force update of results if we have them
          if (!is.null(state$de_genes()) && nrow(state$de_genes()) > 0) {
            # Get current DE results
            results <- state$de_genes()
            
            # Check if we have a comparison column to update
            if ("comparison" %in% colnames(results)) {
              # Get new cluster labels
              new_labels <- cluster_management$getClusterLabels()
              
              # Get all cluster numbers in the data
              all_clusters <- as.character(sort(unique(as.numeric(names(new_labels)))))
              
              # Generate a mapping of old labels to new labels
              label_mapping <- lapply(all_clusters, function(cluster_key) {
                if (cluster_key %in% names(new_labels)) {
                  return(list(
                    old = paste("Cluster", cluster_key),
                    new = new_labels[[cluster_key]]
                  ))
                }
                return(NULL)
              })
              
              # Filter out NULLs
              label_mapping <- label_mapping[!sapply(label_mapping, is.null)]
              
              # Apply the mapping to comparison text if we have mappings
              if (length(label_mapping) > 0) {
                # Get the current comparison text
                current_comparison <- results$comparison[1]
                
                # Apply each mapping
                for (mapping in label_mapping) {
                  current_comparison <- gsub(
                    mapping$old, 
                    mapping$new, 
                    current_comparison, 
                    fixed = TRUE
                  )
                }
                
                # Update all rows
                results$comparison <- current_comparison
                
                # Update the stored results
                state$de_genes(results)
              }
            }
          }
        }
      }
    })
    
    # Track sample and condition filtering
    state$filtered_samples <- reactiveVal(NULL)
    state$filtered_conditions <- reactiveVal(NULL)
    state$condition_column <- reactiveVal(NULL)
    state$filter_update_trigger <- reactiveVal(runif(1))
    
    # Watch for changes in sample management active status
    observe({
      req(sample_management)
      
      # Get active samples
      active_samples <- sample_management$getActiveSampleIds()
      
      # Update filtered samples
      if (!identical(state$filtered_samples(), active_samples)) {
        state$filtered_samples(active_samples)
        
        # Force UI refresh when filters change
        state$filter_update_trigger(runif(1))
      }
    })
    
    # Watch for changes in condition management active status
    observe({
      req(condition_management)
      
      # Get active conditions and the condition column
      active_conditions <- condition_management$getActiveConditions()
      condition_column <- condition_management$getConditionColumn()
      
      # Update filtered conditions and condition column
      if (!identical(state$filtered_conditions(), active_conditions) || 
          !identical(state$condition_column(), condition_column)) {
        
        state$filtered_conditions(active_conditions)
        state$condition_column(condition_column)
        
        # Force UI refresh when filters change
        state$filter_update_trigger(runif(1))
      }
    })
    
    # Track cluster labels from cluster management module
    observe({
      req(cluster_management)
      # This will create a reactive dependency on the cluster labels
      state$cluster_labels <- cluster_management$getClusterLabels()
    })
    
    # Force UI refresh when cluster labels change
    observeEvent(state$cluster_labels, {
      output$analysisUI <- renderUI({
        setupAnalysisUI(ns, clustered_seurat, cluster_management, state, sample_management, condition_management)
      })
    })
    
    # Render Analysis UI
    output$analysisUI <- renderUI({
      setupAnalysisUI(ns, clustered_seurat, cluster_management, state, sample_management, condition_management)
    })
    
    # Setup DE analysis handlers
    setupAnalysisHandlers(input, output, clustered_seurat, cluster_management, state, session, sample_management, condition_management)
    
    setResults <- function(new_results, analysis_type = NULL, heatmap_data = NULL, general_heatmap_genes = NULL) {
      if (!is.null(new_results)) {
        state$de_genes(new_results)
        state$de_status("completed")
        
        # Set the analysis state based on input or determine from results
        if (!is.null(analysis_type)) {
          state$analysis_state(analysis_type)
        } else if ("comparison" %in% colnames(new_results)) {
          # Try to determine if it's one-vs-all or pairwise from results
          if (grepl("vs All", new_results$comparison[1])) {
            state$analysis_state("one_vs_all")
          } else {
            state$analysis_state("pairwise")
          }
        }
        
        # Set heatmap data if provided
        if (!is.null(heatmap_data)) {
          state$heatmap_data(heatmap_data)
          state$heatmap_type("specific")
        }
        
        # Set general heatmap data if provided
        if (!is.null(general_heatmap_genes)) {
          state$general_heatmap_genes(general_heatmap_genes)
          if (state$analysis_state() == "general_heatmap") {
            state$heatmap_type("general")
          }
        }
      }
    }
    
    # Return results, status, and labels
    return(list(
      results = state$de_genes,
      status = state$de_status,
      cluster_labels = reactive({ state$cluster_labels }),
      getAnalysisState = function() { state$analysis_state() },
      getHeatmapData = function() { state$heatmap_data() },
      getGeneralHeatmapGenes = function() { state$general_heatmap_genes() },
      setResults = setResults
    ))
  })
}

#' @title Setup Module Reactive State
#' @description Initializes the reactive state for the DE analysis module.
#' @return A list of reactive values to store module state
#' @keywords internal
setupReactiveState <- function() {
  # Initialize state
  list(
    # Heatmap state
    heatmap_data = reactiveVal(NULL),
    heatmap_type = reactiveVal(NULL),
    heatmap_status_message = NULL,
    
    # DE analysis state
    de_genes = reactiveVal(NULL),
    de_status = reactiveVal(NULL),
    general_heatmap_genes = reactiveVal(NULL),
    general_heatmap_clusters = reactiveVal(NULL),
    analysis_state = reactiveVal("none"),
    cluster_labels = NULL
  )
}

#' @title Setup Analysis UI
#' @description Creates the UI for differential expression analysis based on available clusters.
#' @param ns Namespace function
#' @param clustered_seurat Reactive expression containing the clustered Seurat object
#' @param cluster_management Cluster management module instance
#' @return UI elements for DE analysis
#' @keywords internal
setupAnalysisUI <- function(ns, clustered_seurat, cluster_management, state, sample_management = NULL, condition_management = NULL) {
  # Get cluster choices (existing code)
  active_cluster_ids <- cluster_management$getActiveClusterIds()
  current_labels <- if (!is.null(state$cluster_labels)) {
    state$cluster_labels
  } else {
    cluster_management$getClusterLabels()
  }
  
  cluster_choices <- setNames(
    active_cluster_ids,
    vapply(active_cluster_ids, function(x) {
      if (as.character(x) %in% names(current_labels)) {
        current_labels[[as.character(x)]]
      } else {
        paste("Cluster", x)
      }
    }, character(1))
  )
  
  # Get sample choices
  sample_choices <- NULL
  if (!is.null(sample_management)) {
    active_sample_ids <- sample_management$getActiveSampleIds()
    if (!is.null(active_sample_ids) && length(active_sample_ids) > 0) {
      sample_labels <- sample_management$getSampleLabels()
      
      sample_choices <- setNames(
        active_sample_ids,
        vapply(active_sample_ids, function(x) {
          if (!is.null(sample_labels) && x %in% names(sample_labels)) {
            sample_labels[[x]]
          } else {
            x  # Use the sample ID as its own label if no custom label
          }
        }, character(1))
      )
    }
  }
  
  # Get condition choices
  condition_choices <- NULL
  condition_column <- NULL
  if (!is.null(condition_management)) {
    active_conditions <- condition_management$getActiveConditions()
    condition_column <- condition_management$getConditionColumn()
    
    if (!is.null(active_conditions) && length(active_conditions) > 0 && !is.null(condition_column)) {
      condition_labels <- condition_management$getConditionLabels()
      
      condition_choices <- setNames(
        active_conditions,
        vapply(active_conditions, function(x) {
          if (!is.null(condition_labels) && x %in% names(condition_labels)) {
            condition_labels[[x]]
          } else {
            x  # Use the condition value as its own label if no custom label
          }
        }, character(1))
      )
    }
  }

  # Create analysis UI with all choices
  createDEAnalysisUI(ns, cluster_choices, sample_choices, condition_choices, condition_column)
}

#' @title Setup Analysis Handlers
#' @description Sets up event handlers and renders outputs for DE analysis.
#' @param input Shiny input object
#' @param output Shiny output object
#' @param clustered_seurat Reactive expression containing the clustered Seurat object
#' @param cluster_management Cluster management module instance
#' @param state Reactive state list
#' @param session Shiny session object
#' @keywords internal
setupAnalysisHandlers <- function(input, output, clustered_seurat, cluster_management, state, session, sample_management = NULL, condition_management = NULL) {
  ns <- session$ns
  
  observe({
    req(cluster_management)
    state$cluster_labels <- cluster_management$getClusterLabels()
  })
  
  # Handle unified Cluster comparison
  observeEvent(input$runDECluster, {
    req(clustered_seurat(), input$targetClusters1)
    
    # Get the second group (check if it's "All" or specific clusters)
    clusters2 <- NULL
    if (length(input$targetClusters2) == 1 && input$targetClusters2[1] == "All") {
      clusters2 <- "All"  # Special case
    } else {
      req(input$targetClusters2)
      clusters2 <- input$targetClusters2
    }
    
    # Determine analysis type based on second cluster selection
    if (length(clusters2) == 1 && clusters2[1] == "All") {
      clear_state("one_vs_all")
    } else {
      clear_state("pairwise")
    }
    
    # Run cluster analysis with the unified function that considers all management modules
    runClusterAnalysis(
      clustered_seurat(), 
      input$targetClusters1, 
      clusters2, 
      state,
      cluster_management,
      sample_management,
      condition_management
    )
  })
  
  # Handle Sample comparison
  observeEvent(input$runDESample, {
    req(clustered_seurat(), input$targetSamples1)
    
    # Get the second group (check if it's "All" or specific samples)
    samples2 <- NULL
    if (length(input$targetSamples2) == 1 && input$targetSamples2[1] == "All") {
      samples2 <- "All"  # Special case
    } else {
      req(input$targetSamples2)
      samples2 <- input$targetSamples2
    }
    
    # Check for overlap in selections when not using "All"
    if (samples2 != "All" && any(input$targetSamples1 %in% samples2)) {
      showNotification("Please select different samples for each group.", type = "warning")
      return(NULL)
    }
    
    clear_state("sample_comparison")
    
    # Run sample comparison with integrated filtering
    runSampleAnalysis(
      clustered_seurat(), 
      input$targetSamples1, 
      samples2, 
      state,
      sample_management,
      cluster_management,
      condition_management
    )
  })
  
  # Handle Condition comparison
  observeEvent(input$runDECondition, {
    req(clustered_seurat(), input$targetConditions1, condition_management)
    
    # Get the second group (check if it's "All" or specific conditions)
    conditions2 <- NULL
    if (length(input$targetConditions2) == 1 && input$targetConditions2[1] == "All") {
      conditions2 <- "All"  # Special case
    } else {
      req(input$targetConditions2)
      conditions2 <- input$targetConditions2
    }
    
    # Check for overlap in selections when not using "All"
    if (conditions2 != "All" && any(input$targetConditions1 %in% conditions2)) {
      showNotification("Please select different conditions for each group.", type = "warning")
      return(NULL)
    }
    
    condition_column <- condition_management$getConditionColumn()
    if (is.null(condition_column)) {
      showNotification("No condition column selected.", type = "error")
      return(NULL)
    }
    
    clear_state("condition_comparison")
    
    # Run condition comparison with integrated filtering
    runConditionAnalysis(
      clustered_seurat(), 
      input$targetConditions1, 
      conditions2, 
      condition_column,
      state,
      condition_management,
      cluster_management,
      sample_management
    )
  })
  
  observeEvent(input$runSampleHeatmap, {
    req(clustered_seurat(), input$genesPerSample)
    
    state$heatmap_type("general")
    
    # Run general heatmap analysis but for samples instead of clusters
    # You'll need to implement this function
    runSampleHeatmapAnalysis(
      clustered_seurat(), 
      sample_management$getActiveSampleIds(), 
      input$genesPerSample, 
      state
    )
  })
  
  observeEvent(input$runConditionHeatmap, {
    req(clustered_seurat(), input$genesPerCondition, condition_management)
    
    state$heatmap_type("general")
    
    # Run general heatmap analysis but for conditions instead of clusters
    # You'll need to implement this function
    runConditionHeatmapAnalysis(
      clustered_seurat(), 
      condition_management$getActiveConditions(), 
      condition_management$getConditionColumn(),
      input$genesPerCondition, 
      state
    )
  })
  
  
  # Helper to clear state when starting a new analysis
  clear_state <- function(new_state) {
    state$analysis_state(new_state)
    
    # List of all DE analysis states
    de_analysis_states <- c("one_vs_all", "pairwise", "sample_comparison", "condition_comparison")
    
    # List of all heatmap states
    heatmap_states <- c("general_heatmap", "sample_heatmap", "condition_heatmap")
    
    # Clear visualizations that aren't needed for the new state
    if (new_state %in% de_analysis_states) {
      # Reset general heatmap when switching to DE analysis
      state$general_heatmap_genes(NULL)
      state$general_heatmap_clusters(NULL)
    } else if (new_state %in% heatmap_states) {
      # Clear specific heatmap data when switching to general heatmaps
      state$heatmap_data(NULL)
    }
  }
  
  # Get currently active clusters
  active_cluster_list <- reactive({
    cluster_management$getActiveClusterList()
  })
  
  # Handle General heatmap analysis
  observeEvent(input$runGeneralHeatmap, {
    req(clustered_seurat(), input$genesPerCluster)
    current_active_list <- active_cluster_list()
    
    if (length(current_active_list) == 0) {
      showNotification("No active clusters selected.", type = "warning")
      # Clear any previous heatmap data to prevent errors
      state$general_heatmap_genes(NULL)
      state$general_heatmap_clusters(NULL)
      return(NULL)
    }
    
    clear_state("general_heatmap")
    
    # Run general heatmap analysis
    runGeneralHeatmapAnalysis(
      clustered_seurat(), 
      current_active_list, 
      input$genesPerCluster, 
      state
    )
  })
  
  # Handle heatmap generation button
  observeEvent(input$generate_heatmap, {
    req(clustered_seurat(), input$top_n_genes, state$de_genes())
    
    generateDEHeatmap(
      clustered_seurat(), 
      input$top_n_genes, 
      state$de_genes(), 
      state
    )
  })
  
  # Render DE results UI
  output$deResultsUI <- renderUI({
    condition_column <- NULL
    if (!is.null(condition_management) && is.function(condition_management$getConditionColumn)) {
      condition_column <- condition_management$getConditionColumn()
    }
    renderDEResultsUI(ns, state, clustered_seurat(), condition_column)
  })
  
  # Render general heatmap UI
  output$generalHeatmapUI <- renderUI({
    renderGeneralHeatmapUI(ns, state, active_cluster_list)
  })
  
  # Single volcano plot
  volcano_plot <- reactive({
    req(state$de_genes())
    results <- state$de_genes()
    if (nrow(results) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No significant DE genes found") + 
               theme_void())
    }
    
    createVolcanoPlot(results)
  })
  
  # Render volcano plot
  output$volcanoPlot <- renderPlot({
    volcano_plot()
  })
  
  # Download handler for volcano plot
  output$downloadVolcanoPlot <- downloadHandler(
    filename = function() {
      paste("volcano_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      # Create high-resolution PNG device directly
      png(file, width = 3000, height = 2400, res = 300)
      
      # Generate the plot directly to the device
      results <- state$de_genes()
      print(createVolcanoPlot(results))
      
      # Close the device to save the file
      dev.off()
    }
  )
  
  # DE results table
  output$deTable <- DT::renderDataTable({
    req(state$de_genes())
    results <- state$de_genes()
    
    # Sort by p-value adjusted and then fold change
    results <- results[order(results$p_val_adj, -abs(results$avg_log2FC)), ]
    
    # Remove the filter_context column from display (but keep it in the data)
    display_results <- results
    if ("filter_context" %in% colnames(display_results)) {
      display_results <- display_results[, !colnames(display_results) %in% "filter_context", drop = FALSE]
    }
    
    DT::datatable(display_results,
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE,
                    dom = 'lBfrtip',
                    buttons = c('csv')
                  ),
                  extensions = 'Buttons',
                  rownames = FALSE) %>%
      DT::formatSignif(columns = c("p_val", "p_val_adj"), digits = 3) %>%
      DT::formatRound(columns = c("avg_log2FC"), digits = 2)
  })
  
  # Specific heatmap plot with comprehensive error handling
  specific_heatmap_plot <- reactive({
    req(state$heatmap_data())
    req(state$heatmap_type() == "specific")
    req(clustered_seurat())
    
    # Get the latest cluster labels from state first, then from cluster_management
    cluster_labels <- if (!is.null(state$cluster_labels)) {
      state$cluster_labels
    } else if (!is.null(cluster_management)) {
      tryCatch({
        cluster_management$getClusterLabels()
      }, error = function(e) {
        NULL
      })
    } else {
      NULL
    }
    
    # Final fallback to default labels if still null
    if (is.null(cluster_labels)) {
      unique_clusters <- sort(unique(clustered_seurat()$seurat_clusters))
      cluster_labels <- setNames(
        paste("Cluster", unique_clusters),
        as.character(unique_clusters)
      )
    }
    
    # Get active clusters list
    active_clusters <- tryCatch({
      if (is.function(active_cluster_list)) {
        as.numeric(active_cluster_list())
      } else {
        NULL
      }
    }, error = function(e) {
      NULL
    })
    
    # Determine what grouping to use based on the analysis type
    group_by <- "seurat_clusters"  # Default
    
    # See if we have a stored group_by in state
    if (!is.null(state$group_by)) {
      group_by <- state$group_by
    } else {
      # Try to determine from the heatmap data or analysis state
      current_state <- state$analysis_state()
      if (current_state == "sample_comparison") {
        group_by <- "sample"
      } else if (current_state == "condition_comparison" && !is.null(condition_management)) {
        # Get condition column from condition management
        condition_col <- condition_management$getConditionColumn()
        if (!is.null(condition_col)) {
          group_by <- condition_col
        }
      }
    }
    
    # For non-cluster groupings, get the appropriate active items and labels
    if (group_by == "sample" && !is.null(sample_management)) {
      # Use sample management for active samples and labels
      active_items <- sample_management$getActiveSampleIds()
      sample_labels <- sample_management$getSampleLabels()
      
      # Now call the function with robust error handling
      tryCatch({
        createExpressionHeatmap(
          clustered_seurat(),
          state$heatmap_data(),
          sample_labels,  # Use sample labels
          active_items,   # Use active samples
          group_by = "sample"  # Group by sample
        )
      }, error = function(e) {
        # Return a simple error message plot
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = paste("Error creating heatmap:", e$message)) + 
          theme_void()
      })
    } else if (group_by != "seurat_clusters" && group_by != "sample" && 
               !is.null(condition_management) && 
               group_by == condition_management$getConditionColumn()) {
      # Use condition management for active conditions and labels
      active_items <- condition_management$getActiveConditions()
      condition_labels <- condition_management$getConditionLabels()
      
      # Now call the function with robust error handling
      tryCatch({
        createExpressionHeatmap(
          clustered_seurat(),
          state$heatmap_data(),
          condition_labels,  # Use condition labels
          active_items,      # Use active conditions
          group_by = group_by  # Group by condition column
        )
      }, error = function(e) {
        # Return a simple error message plot
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = paste("Error creating heatmap:", e$message)) + 
          theme_void()
      })
    } else {
      # Default cluster-based heatmap
      tryCatch({
        createExpressionHeatmap(
          clustered_seurat(),
          state$heatmap_data(),
          cluster_labels,
          active_clusters,
          group_by = "seurat_clusters"  # Group by clusters
        )
      }, error = function(e) {
        # Return a simple error message plot
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = paste("Error creating heatmap:", e$message)) + 
          theme_void()
      })
    }
  })
  
  general_heatmap_plot <- reactive({
    req(state$general_heatmap_genes())
    req(clustered_seurat())
    
    # Get the latest cluster labels from state first, then from cluster_management
    cluster_labels <- if (!is.null(state$cluster_labels)) {
      state$cluster_labels
    } else if (!is.null(cluster_management)) {
      tryCatch({
        cluster_management$getClusterLabels()
      }, error = function(e) {
        NULL
      })
    } else {
      NULL
    }
    
    # Final fallback to default labels if still null
    if (is.null(cluster_labels)) {
      unique_clusters <- sort(unique(clustered_seurat()$seurat_clusters))
      cluster_labels <- setNames(
        paste("Cluster", unique_clusters),
        as.character(unique_clusters)
      )
    }
    
    # Get active items based on the analysis type
    active_items <- NULL
    current_state <- state$analysis_state()
    
    if (current_state == "general_heatmap") {
      # Get active clusters for cluster-based heatmap
      active_items <- tryCatch({
        if (is.function(active_cluster_list)) {
          as.numeric(active_cluster_list())
        } else {
          NULL
        }
      }, error = function(e) {
        NULL
      })
      
      if (is.null(active_items) || length(active_items) == 0) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "No active clusters selected") + 
                 theme_void())
      }
      
      # Call createGeneralHeatmap with robust error handling
      tryCatch({
        createGeneralHeatmap(
          clustered_seurat(),
          state$general_heatmap_genes(),
          cluster_labels,
          as.numeric(active_items),
          state$general_heatmap_clusters()
        )
      }, error = function(e) {
        # Return a simple error message plot
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = paste("Error creating general heatmap:", e$message)) + 
          theme_void()
      })
    } else if (current_state == "sample_heatmap") {
      # For sample heatmap, use the sample labels from sample_management if available
      sample_labels <- NULL
      if (!is.null(sample_management) && is.function(sample_management$getSampleLabels)) {
        sample_labels <- sample_management$getSampleLabels()
      }
      
      # Call createSampleHeatmap with sample labels
      tryCatch({
        # Make a copy of the Seurat object to modify for this plot
        seurat_obj <- clustered_seurat()
        # Set sample as the identity for heatmap grouping
        Idents(seurat_obj) <- "sample"
        
        # Get active samples for filtering
        active_samples <- NULL
        if (!is.null(sample_management) && is.function(sample_management$getActiveSampleIds)) {
          active_samples <- sample_management$getActiveSampleIds()
        }
        
        # Create the general heatmap, using samples as the groups
        createGeneralHeatmap(
          seurat_obj,
          state$general_heatmap_genes(),
          sample_labels,
          NULL,  # No cluster filtering needed
          NULL   # No cluster ordering needed
        )
      }, error = function(e) {
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = paste("Error creating sample heatmap:", e$message)) + 
          theme_void()
      })
    } else if (current_state == "condition_heatmap") {
      # For condition heatmap, get condition column and labels
      condition_column <- NULL
      condition_labels <- NULL
      
      if (!is.null(condition_management)) {
        if (is.function(condition_management$getConditionColumn)) {
          condition_column <- condition_management$getConditionColumn()
        }
        if (is.function(condition_management$getConditionLabels)) {
          condition_labels <- condition_management$getConditionLabels()
        }
      }
      
      # Skip if no condition column is set
      if (is.null(condition_column)) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "No condition column selected") + 
                 theme_void())
      }
      
      # Call createConditionHeatmap with condition labels
      tryCatch({
        # Make a copy of the Seurat object to modify for this plot
        seurat_obj <- clustered_seurat()
        # Set the condition column as the identity for heatmap grouping
        Idents(seurat_obj) <- condition_column
        
        # Get active conditions for filtering
        active_conditions <- NULL
        if (!is.null(condition_management) && is.function(condition_management$getActiveConditions)) {
          active_conditions <- condition_management$getActiveConditions()
        }
        
        # Create the general heatmap, using conditions as the groups
        createGeneralHeatmap(
          seurat_obj,
          state$general_heatmap_genes(),
          condition_labels,
          NULL,  # No cluster filtering needed
          NULL   # No cluster ordering needed
        )
      }, error = function(e) {
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = paste("Error creating condition heatmap:", e$message)) + 
          theme_void()
      })
    } else {
      # Default error message if unknown state
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = "Unknown heatmap type requested") + 
        theme_void()
    }
  })
  
  # Render heatmap plots
  output$heatmapPlot <- renderPlot({
    specific_heatmap_plot()
  })
  
  output$generalHeatmapPlot <- renderPlot({
    general_heatmap_plot()
  })
  
  # Download handlers for heatmaps
  setupHeatmapDownloadHandlers(output, state, clustered_seurat, active_cluster_list)
  
  setupGeneSearch(input, output, session, clustered_seurat)
  
  # Add this download handler for the boxplot
  output$download_boxplot <- downloadHandler(
    filename = function() {
      gene_info <- if (!is.null(input$boxplot_gene) && input$boxplot_gene != "") {
        input$boxplot_gene
      } else {
        "gene"
      }
      paste0("gene_expression_", gene_info, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      # Create high-resolution PNG
      png(file, width = 3000, height = 2400, res = 300)
      
      # Get the gene ID
      gene_id <- input$boxplot_gene
      
      # Get group variable
      group_var <- input$boxplot_group
      
      # Get active groups for filtering
      active_groups <- NULL
      
      if (group_var == "seurat_clusters") {
        # Get active clusters from cluster management module
        if (!is.null(cluster_management) && is.function(cluster_management$getActiveClusterIds)) {
          active_groups <- cluster_management$getActiveClusterIds()
        } else {
          # Fallback - use all clusters
          active_groups <- unique(clustered_seurat()$seurat_clusters)
        }
      } else if (group_var == "sample") {
        # Get active samples from sample management module
        if (!is.null(sample_management) && is.function(sample_management$getActiveSampleIds)) {
          active_groups <- sample_management$getActiveSampleIds()
        } else {
          # Fallback - use all samples
          active_groups <- unique(clustered_seurat()$sample)
        }
      }
      
      # Get labels if available
      group_labels <- NULL
      if (group_var == "seurat_clusters" && 
          !is.null(cluster_management) && 
          is.function(cluster_management$getClusterLabels)) {
        group_labels <- cluster_management$getClusterLabels()
      }
      
      # Determine if we should show statistics
      show_stats <- input$boxplot_stats
      
      # Create comparisons list for statistics if needed
      comparisons <- NULL
      if (show_stats && group_var == "seurat_clusters" && !is.null(active_groups) && length(active_groups) > 1) {
        # If we have many groups, just compare adjacent ones to avoid cluttering
        if (length(active_groups) > 3) {
          comparisons <- lapply(1:(length(active_groups)-1), function(i) {
            c(as.character(active_groups[i]), as.character(active_groups[i+1]))
          })
        } else {
          # For few groups, compare all pairs
          comparisons <- combn(as.character(active_groups), 2, simplify = FALSE)
        }
      }
      
      # Generate and save the plot
      p <- createExpressionBoxplot(
        clustered_seurat(),
        gene_id = gene_id,
        group_by = group_var,
        comparisons = comparisons,
        active_groups = active_groups,
        cluster_labels = group_labels
      )
      
      print(p)
      dev.off()
    }
  )
  
  # Download handler for boxplot data
  output$download_boxplot_data <- downloadHandler(
    filename = function() {
      gene_info <- if (!is.null(input$boxplot_gene) && input$boxplot_gene != "") {
        input$boxplot_gene
      } else {
        "gene"
      }
      group_by <- input$boxplot_group
      paste0("gene_expression_", gene_info, "_by_", group_by, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(clustered_seurat())
      
      # Get the gene ID
      gene_id <- input$boxplot_gene
      
      # Check if gene exists
      if (is.null(gene_id) || gene_id == "" || !(gene_id %in% rownames(clustered_seurat()))) {
        write.csv(data.frame(error = "Gene not found in dataset"), file, row.names = FALSE)
        return()
      }
      
      # Get group variable
      group_var <- input$boxplot_group
      
      # Get active groups for filtering
      active_groups <- NULL
      
      if (group_var == "seurat_clusters") {
        # Get active clusters from cluster management module
        if (!is.null(cluster_management) && is.function(cluster_management$getActiveClusterIds)) {
          active_groups <- cluster_management$getActiveClusterIds()
        } else {
          # Fallback - use all clusters
          active_groups <- unique(clustered_seurat()$seurat_clusters)
        }
      } else if (group_var == "sample") {
        # Get active samples from sample management module
        if (!is.null(sample_management) && is.function(sample_management$getActiveSampleIds)) {
          active_groups <- sample_management$getActiveSampleIds()
        } else {
          # Fallback - use all samples
          active_groups <- unique(clustered_seurat()$sample)
        }
      } else if (!is.null(condition_management) && 
                 group_var == condition_management$getConditionColumn()) {
        # Get active conditions from condition management module
        active_groups <- condition_management$getActiveConditions()
      }
      
      # Extract boxplot data
      boxplot_data <- extractBoxplotData(
        clustered_seurat(),
        gene_id = gene_id,
        group_by = group_var,
        active_groups = active_groups
      )
      
      if (!is.null(boxplot_data) && nrow(boxplot_data) > 0) {
        # Add metadata header with attributes
        meta_info <- c(
          paste0("# Gene Expression Data Export"),
          paste0("# Date: ", Sys.time()),
          paste0("# Gene ID: ", attr(boxplot_data, "gene_id")),
          paste0("# Gene Symbol: ", attr(boxplot_data, "gene_symbol")),
          paste0("# Grouped by: ", attr(boxplot_data, "group_by"))
        )
        
        # Write metadata and data
        con <- file(file, "w")
        writeLines(meta_info, con)
        close(con)
        
        # Append the data
        write.table(boxplot_data, file, sep = ",", row.names = FALSE, col.names = TRUE, 
                    append = TRUE, quote = TRUE)
      } else {
        # Create empty file with message
        write.csv(data.frame(error = "No expression data available for this gene"), 
                  file, row.names = FALSE)
      }
    }
  )
  
  # Generate the boxplot with improved error handling
  output$gene_boxplot <- renderPlot({
    req(clustered_seurat())
    
    # Get the gene ID
    gene_id <- input$boxplot_gene
    
    # Check if gene_id is empty
    if (is.null(gene_id) || gene_id == "") {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "Please select a gene from the dropdown or type to search") + 
               theme_void())
    }
    
    # Check if gene exists in the dataset
    if (!(gene_id %in% rownames(clustered_seurat()))) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = paste("Gene", gene_id, "not found in dataset")) + 
               theme_void())
    }
    
    # Get group variable
    group_var <- input$boxplot_group
    
    # Get active groups for filtering
    active_groups <- NULL
    
    if (group_var == "seurat_clusters") {
      # Get active clusters
      tryCatch({
        if (!is.null(cluster_management) && is.function(cluster_management$getActiveClusterIds)) {
          active_groups <- cluster_management$getActiveClusterIds()
        }
      }, error = function(e) {
        print(paste("Error getting active clusters:", e$message))
      })
    } else if (group_var == "sample") {
      # Get active samples
      tryCatch({
        if (!is.null(sample_management) && is.function(sample_management$getActiveSampleIds)) {
          active_groups <- sample_management$getActiveSampleIds()
        }
      }, error = function(e) {
        print(paste("Error getting active samples:", e$message))
      })
    } else if (!is.null(condition_management) && 
               group_var == condition_management$getConditionColumn()) {
      # Get active conditions
      tryCatch({
        active_groups <- condition_management$getActiveConditions()
      }, error = function(e) {
        print(paste("Error getting active conditions:", e$message))
      })
    }
    
    # Get cluster labels if available
    group_labels <- NULL
    if (group_var == "seurat_clusters") {
      tryCatch({
        if (!is.null(cluster_management) && is.function(cluster_management$getClusterLabels)) {
          group_labels <- cluster_management$getClusterLabels()
        }
      }, error = function(e) {
        print(paste("Error getting cluster labels:", e$message))
      })
    }
    
    # Determine if we should show statistics
    show_stats <- input$boxplot_stats
    
    # Create the boxplot with robust error handling
    tryCatch({
      createExpressionBoxplot(
        clustered_seurat(),
        gene_id = gene_id,
        group_by = group_var,
        comparisons = if(show_stats) TRUE else NULL,
        active_groups = active_groups,
        cluster_labels = group_labels
      )
    }, error = function(e) {
      print(paste("Error in boxplot generation:", e$message))
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Error generating boxplot:", e$message)) + 
        theme_void()
    })
  })
}

#' @title Apply Global Filtering to Seurat Object
#' @description Applies filtering based on active clusters, samples, and conditions
#' @param seurat_obj Seurat object to filter
#' @param cluster_management Cluster management module
#' @param sample_management Sample management module 
#' @param condition_management Condition management module
#' @return Filtered Seurat object
#' @keywords internal
applyGlobalFiltering <- function(seurat_obj, cluster_management, sample_management, condition_management) {
  filtered_obj <- seurat_obj
  cells_to_keep <- rep(TRUE, ncol(filtered_obj))
  filter_applied <- FALSE
  
  # Apply cluster filtering if available
  if (!is.null(cluster_management)) {
    active_clusters <- cluster_management$getActiveClusterIds()
    if (!is.null(active_clusters) && length(active_clusters) > 0) {
      cluster_cells <- filtered_obj$seurat_clusters %in% active_clusters
      cells_to_keep <- cells_to_keep & cluster_cells
      filter_applied <- TRUE
    }
  }
  
  # Apply sample filtering if available
  if (!is.null(sample_management)) {
    active_samples <- sample_management$getActiveSampleIds()
    if (!is.null(active_samples) && length(active_samples) > 0) {
      sample_cells <- filtered_obj$sample %in% active_samples
      cells_to_keep <- cells_to_keep & sample_cells
      filter_applied <- TRUE
    }
  }
  
  # Apply condition filtering if available
  if (!is.null(condition_management)) {
    condition_column <- condition_management$getConditionColumn()
    active_conditions <- condition_management$getActiveConditions()
    
    if (!is.null(condition_column) && !is.null(active_conditions) && 
        length(active_conditions) > 0 && 
        condition_column %in% colnames(filtered_obj@meta.data)) {
      
      condition_cells <- filtered_obj@meta.data[[condition_column]] %in% active_conditions
      cells_to_keep <- cells_to_keep & condition_cells
      filter_applied <- TRUE
    }
  }
  
  # Only subset if filtering was applied and there are cells to keep
  if (filter_applied && sum(cells_to_keep) > 0) {
    filtered_obj <- subset(filtered_obj, cells = colnames(filtered_obj)[cells_to_keep])
  } else if (filter_applied && sum(cells_to_keep) == 0) {
    # Return NULL if no cells match all filters
    return(NULL)
  }
  
  return(filtered_obj)
}

#' @title Run Cluster DE Analysis
#' @description Performs differential expression analysis comparing groups of clusters.
#' @param seurat_obj Seurat object containing the data
#' @param clusters1 First group of clusters for comparison (vector)
#' @param clusters2 Second group of clusters for comparison (can include "All")
#' @param state Reactive state list
#' @param cluster_management Cluster management module instance
#' @param sample_management Sample management module instance
#' @param condition_management Condition management module instance
#' @keywords internal
runClusterAnalysis <- function(seurat_obj, clusters1, clusters2, state, 
                               cluster_management, sample_management = NULL, condition_management = NULL) {
  withProgress(message = 'Computing cluster differential expression...', {
    # First apply global filtering based on active items in all management modules
    filtered_obj <- applyGlobalFiltering(seurat_obj, cluster_management, sample_management, condition_management)
    
    # Check if there are cells left after filtering
    if (is.null(filtered_obj) || ncol(filtered_obj) == 0) {
      showNotification("No cells match the active filtering criteria.", type = "error")
      return(NULL)
    }
    
    # Get active clusters for limiting the scope
    active_clusters <- cluster_management$getActiveClusterIds()
    if (is.null(active_clusters) || length(active_clusters) == 0) {
      showNotification("No active clusters selected.", type = "warning")
      return(NULL)
    }
    
    # Convert cluster inputs to numeric vectors
    clusters1 <- as.numeric(clusters1)
    
    # Handle "All" case in clusters2
    is_one_vs_all <- FALSE
    if (length(clusters2) == 1 && clusters2[1] == "All") {
      is_one_vs_all <- TRUE
      # Use all active clusters except those in clusters1
      clusters2 <- setdiff(active_clusters, clusters1)
      
      if (length(clusters2) == 0) {
        showNotification("No other active clusters available for comparison.", type = "warning")
        return(NULL)
      }
    } else {
      # Convert to numeric for standard comparison
      clusters2 <- as.numeric(clusters2)
    }
    
    # Check for empty selections
    if (length(clusters1) == 0) {
      showNotification("Please select at least one cluster for the first group.", type = "warning")
      return(NULL)
    }
    if (length(clusters2) == 0) {
      showNotification("Please select at least one cluster for the second group.", type = "warning")
      return(NULL)
    }
    
    # Check for overlap between groups
    if (any(clusters1 %in% clusters2)) {
      showNotification("The two cluster groups should not overlap. Please select different clusters for each group.", type = "warning")
      return(NULL)
    }
    
    # Check if all selected clusters exist in the filtered data
    all_clusters <- unique(filtered_obj$seurat_clusters)
    if (!all(clusters1 %in% all_clusters) || !all(clusters2 %in% all_clusters)) {
      showNotification("One or more selected clusters are not present after filtering.", type = "error")
      return(NULL)
    }
    
    # Create a temporary grouping factor for FindMarkers
    filtered_obj$temp_group <- ifelse(
      filtered_obj$seurat_clusters %in% clusters1, 
      "group1", 
      ifelse(filtered_obj$seurat_clusters %in% clusters2, "group2", "other")
    )
    
    # Keep only cells from the selected groups
    cells_to_keep <- filtered_obj$temp_group %in% c("group1", "group2")
    if (sum(cells_to_keep) == 0) {
      showNotification("No cells match the selected cluster groups.", type = "error")
      return(NULL)
    }
    
    # Subset the data
    seurat_subset <- subset(filtered_obj, cells = colnames(filtered_obj)[cells_to_keep])
    
    # Set the identity to our temporary grouping
    Idents(seurat_subset) <- "temp_group"
    
    # Run DE analysis
    de_results <- FindMarkers(
      seurat_subset,
      ident.1 = "group1",
      ident.2 = "group2",
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    
    # Check results
    if (is.null(de_results) || nrow(de_results) == 0) {
      showNotification("No differential expression results found.", type = "warning")
      return(NULL)
    }
    
    # Add gene names
    de_results <- addGeneNames(de_results, seurat_obj)
    
    # Get cluster labels
    current_labels <- if (!is.null(state$cluster_labels)) {
      state$cluster_labels
    } else if (!is.null(cluster_management)) {
      cluster_management$getClusterLabels()
    } else {
      NULL
    }
    
    # Create group labels using cluster names
    group1_labels <- sapply(clusters1, function(cl) {
      cluster_str <- as.character(cl)
      if (!is.null(current_labels) && cluster_str %in% names(current_labels)) {
        current_labels[[cluster_str]]
      } else {
        paste("Cluster", cl)
      }
    })
    
    if (is_one_vs_all) {
      comparison_text <- paste0(
        "Clusters (", paste(group1_labels, collapse = ", "), ") vs All Others"
      )
    } else {
      group2_labels <- sapply(clusters2, function(cl) {
        cluster_str <- as.character(cl)
        if (!is.null(current_labels) && cluster_str %in% names(current_labels)) {
          current_labels[[cluster_str]]
        } else {
          paste("Cluster", cl)
        }
      })
      
      comparison_text <- paste0(
        "Clusters (", paste(group1_labels, collapse = ", "), 
        ") vs Clusters (", paste(group2_labels, collapse = ", "), ")"
      )
    }
    
    # Add comparison label
    de_results$comparison <- comparison_text
    
    # Add additional context about current filtering
    de_results$filter_context <- buildFilterContext(sample_management, condition_management, cluster_management)
    
    # Update state
    state$de_genes(de_results)
    state$de_status("completed")
    state$heatmap_data(NULL)
    state$heatmap_type(NULL)
  })
}

#' @title Run Sample DE Analysis
#' @description Performs differential expression analysis comparing groups of samples.
#' @param seurat_obj Seurat object containing the data
#' @param samples1 First group of samples for comparison (vector)
#' @param samples2 Second group of samples for comparison (can include "All")
#' @param state Reactive state list
#' @param sample_management Sample management module instance
#' @param cluster_management Cluster management module instance
#' @param condition_management Condition management module instance
#' @keywords internal
runSampleAnalysis <- function(seurat_obj, samples1, samples2, state, 
                              sample_management, cluster_management = NULL, condition_management = NULL) {
  withProgress(message = 'Computing sample comparison...', {
    # First apply global filtering based on active items in all management modules
    filtered_obj <- applyGlobalFiltering(seurat_obj, cluster_management, sample_management, condition_management)
    
    # Check if there are cells left after filtering
    if (is.null(filtered_obj) || ncol(filtered_obj) == 0) {
      showNotification("No cells match the active filtering criteria.", type = "error")
      return(NULL)
    }
    
    # Get all active samples
    all_active_samples <- sample_management$getActiveSampleIds()
    if (is.null(all_active_samples) || length(all_active_samples) == 0) {
      showNotification("No active samples available.", type = "error")
      return(NULL)
    }
    
    # Check for empty selections
    if (length(samples1) == 0) {
      showNotification("Please select at least one sample for the first group.", type = "warning")
      return(NULL)
    }
    
    # Handle "All" case
    is_one_vs_all <- FALSE
    if (length(samples2) == 1 && samples2[1] == "All") {
      is_one_vs_all <- TRUE
      # Use all active samples except those in samples1
      samples2 <- setdiff(all_active_samples, samples1)
      
      if (length(samples2) == 0) {
        showNotification("No other active samples to compare against.", type = "error")
        return(NULL)
      }
    } else if (length(samples2) == 0) {
      showNotification("Please select at least one sample for the second group.", type = "warning")
      return(NULL)
    }
    
    # Check for overlap between groups
    if (any(samples1 %in% samples2)) {
      showNotification("The two sample groups should not overlap. Please select different samples for each group.", type = "warning")
      return(NULL)
    }
    
    # Create a temporary grouping column for the analysis
    filtered_obj$temp_sample_group <- ifelse(
      filtered_obj$sample %in% samples1, 
      "group1", 
      ifelse(filtered_obj$sample %in% samples2, "group2", "other")
    )
    
    # Keep only cells from the selected groups
    cells_to_keep <- filtered_obj$temp_sample_group %in% c("group1", "group2")
    
    # Check if any cells match the selection
    if (sum(cells_to_keep) == 0) {
      showNotification("No cells found for the selected samples after filtering.", type = "error")
      return(NULL)
    }
    
    # Subset data to include only the selected groups
    seurat_subset <- subset(filtered_obj, cells = colnames(filtered_obj)[cells_to_keep])
    
    # Set identity to the temporary grouping
    Idents(seurat_subset) <- "temp_sample_group"
    
    # Run DE analysis
    de_results <- FindMarkers(
      seurat_subset,
      ident.1 = "group1",
      ident.2 = "group2",
      min.pct = 0.25,
      logfc.threshold = 0
    )
    
    # Check results
    if (is.null(de_results) || nrow(de_results) == 0) {
      showNotification("No differential expression results found.", type = "warning")
      return(NULL)
    }
    
    # Add gene names
    de_results <- addGeneNames(de_results, seurat_obj)
    
    # Get sample labels
    sample_labels <- sample_management$getSampleLabels()
    
    # Create group labels using sample names/labels
    group1_labels <- sapply(samples1, function(sample) {
      if (!is.null(sample_labels) && sample %in% names(sample_labels)) {
        sample_labels[[sample]]
      } else {
        sample
      }
    })
    
    if (is_one_vs_all) {
      comparison_text <- paste0(
        "Samples (", paste(group1_labels, collapse = ", "), ") vs All Others"
      )
    } else {
      group2_labels <- sapply(samples2, function(sample) {
        if (!is.null(sample_labels) && sample %in% names(sample_labels)) {
          sample_labels[[sample]]
        } else {
          sample
        }
      })
      
      comparison_text <- paste0(
        "Samples (", paste(group1_labels, collapse = ", "), 
        ") vs Samples (", paste(group2_labels, collapse = ", "), ")"
      )
    }
    
    # Add comparison label
    de_results$comparison <- comparison_text
    
    # Add additional context about current filtering
    de_results$filter_context <- buildFilterContext(sample_management, condition_management, cluster_management)
    
    # Update state
    state$de_genes(de_results)
    state$de_status("completed")
    state$heatmap_data(NULL)
    state$heatmap_type(NULL)
  })
}

#' @title Run Condition DE Analysis
#' @description Performs differential expression analysis comparing groups of conditions.
#' @param seurat_obj Seurat object containing the data
#' @param conditions1 First group of conditions for comparison (vector)
#' @param conditions2 Second group of conditions for comparison (can include "All")
#' @param condition_column Column name containing condition information
#' @param state Reactive state list
#' @param condition_management Condition management module instance
#' @param cluster_management Cluster management module instance
#' @param sample_management Sample management module instance
#' @keywords internal
runConditionAnalysis <- function(seurat_obj, conditions1, conditions2, condition_column, state, 
                                 condition_management, cluster_management = NULL, sample_management = NULL) {
  withProgress(message = 'Computing condition comparison...', {
    # First apply global filtering based on active items in all management modules
    filtered_obj <- applyGlobalFiltering(seurat_obj, cluster_management, sample_management, condition_management)
    
    # Check if there are cells left after filtering
    if (is.null(filtered_obj) || ncol(filtered_obj) == 0) {
      showNotification("No cells match the active filtering criteria.", type = "error")
      return(NULL)
    }
    
    # Get all active conditions
    all_active_conditions <- condition_management$getActiveConditions()
    if (is.null(all_active_conditions) || length(all_active_conditions) == 0) {
      showNotification("No active conditions available.", type = "error")
      return(NULL)
    }
    
    # Check for empty selections
    if (length(conditions1) == 0) {
      showNotification("Please select at least one condition for the first group.", type = "warning")
      return(NULL)
    }
    
    # Handle "All" case
    is_one_vs_all <- FALSE
    if (length(conditions2) == 1 && conditions2[1] == "All") {
      is_one_vs_all <- TRUE
      # Use all active conditions except those in conditions1
      conditions2 <- setdiff(all_active_conditions, conditions1)
      
      if (length(conditions2) == 0) {
        showNotification("No other active conditions to compare against.", type = "error")
        return(NULL)
      }
    } else if (length(conditions2) == 0) {
      showNotification("Please select at least one condition for the second group.", type = "warning")
      return(NULL)
    }
    
    # Check for overlap between groups
    if (any(conditions1 %in% conditions2)) {
      showNotification("The two condition groups should not overlap. Please select different conditions for each group.", type = "warning")
      return(NULL)
    }
    
    # Create a temporary grouping column for the analysis
    filtered_obj$temp_condition_group <- ifelse(
      filtered_obj@meta.data[[condition_column]] %in% conditions1, 
      "group1", 
      ifelse(filtered_obj@meta.data[[condition_column]] %in% conditions2, "group2", "other")
    )
    
    # Keep only cells from the selected groups
    cells_to_keep <- filtered_obj$temp_condition_group %in% c("group1", "group2")
    
    # Skip if no cells match
    if (sum(cells_to_keep) == 0) {
      showNotification("No cells found for the selected conditions after filtering.", type = "error")
      return(NULL)
    }
    
    # Subset data
    seurat_subset <- subset(filtered_obj, cells = colnames(filtered_obj)[cells_to_keep])
    
    # Set identity to the temporary grouping
    Idents(seurat_subset) <- "temp_condition_group"
    
    # Run DE analysis
    de_results <- FindMarkers(
      seurat_subset,
      ident.1 = "group1",
      ident.2 = "group2",
      min.pct = 0.25,
      logfc.threshold = 0
    )
    
    # Check results
    if (is.null(de_results) || nrow(de_results) == 0) {
      showNotification("No differential expression results found.", type = "warning")
      return(NULL)
    }
    
    # Add gene names
    de_results <- addGeneNames(de_results, seurat_obj)
    
    # Get condition labels
    condition_labels <- condition_management$getConditionLabels()
    
    # Create group labels using condition values/labels
    group1_labels <- sapply(conditions1, function(condition) {
      if (!is.null(condition_labels) && condition %in% names(condition_labels)) {
        condition_labels[[condition]]
      } else {
        condition
      }
    })
    
    if (is_one_vs_all) {
      comparison_text <- paste0(
        "Conditions (", paste(group1_labels, collapse = ", "), ") vs All Others"
      )
    } else {
      group2_labels <- sapply(conditions2, function(condition) {
        if (!is.null(condition_labels) && condition %in% names(condition_labels)) {
          condition_labels[[condition]]
        } else {
          condition
        }
      })
      
      comparison_text <- paste0(
        "Conditions (", paste(group1_labels, collapse = ", "), 
        ") vs Conditions (", paste(group2_labels, collapse = ", "), ")"
      )
    }
    
    # Add comparison label
    de_results$comparison <- comparison_text
    
    # Add additional context about current filtering
    de_results$filter_context <- buildFilterContext(sample_management, condition_management, cluster_management)
    
    # Update state
    state$de_genes(de_results)
    state$de_status("completed")
    state$heatmap_data(NULL)
    state$heatmap_type(NULL)
  })
}

#' @title Build Filtering Context String
#' @description Creates a descriptive string of the current filtering state
#' @param sample_management Sample management module
#' @param condition_management Condition management module
#' @param cluster_management Cluster management module
#' @return String describing active filters
#' @keywords internal
buildFilterContext <- function(sample_management, condition_management, cluster_management) {
  context_parts <- c()
  
  # Add cluster context
  if (!is.null(cluster_management)) {
    active_clusters <- cluster_management$getActiveClusterIds()
    if (!is.null(active_clusters) && length(active_clusters) > 0) {
      if (length(active_clusters) == 1) {
        cluster_context <- paste("Restricted to cluster", active_clusters)
      } else {
        cluster_context <- paste("Restricted to", length(active_clusters), "clusters")
      }
      context_parts <- c(context_parts, cluster_context)
    }
  }
  
  # Add sample context
  if (!is.null(sample_management)) {
    active_samples <- sample_management$getActiveSampleIds()
    if (!is.null(active_samples) && length(active_samples) > 0) {
      if (length(active_samples) == 1) {
        sample_context <- paste("Restricted to sample", active_samples)
      } else {
        sample_context <- paste("Restricted to", length(active_samples), "samples")
      }
      context_parts <- c(context_parts, sample_context)
    }
  }
  
  # Add condition context
  if (!is.null(condition_management)) {
    condition_column <- condition_management$getConditionColumn()
    active_conditions <- condition_management$getActiveConditions()
    if (!is.null(condition_column) && !is.null(active_conditions) && length(active_conditions) > 0) {
      if (length(active_conditions) == 1) {
        condition_context <- paste("Restricted to condition", active_conditions)
      } else {
        condition_context <- paste("Restricted to", length(active_conditions), "conditions")
      }
      context_parts <- c(context_parts, condition_context)
    }
  }
  
  # Combine all parts
  if (length(context_parts) > 0) {
    return(paste(context_parts, collapse = "; "))
  } else {
    return("No filtering applied")
  }
}

#' @title Run General Heatmap Analysis
#' @description Generates a heatmap of top marker genes for each cluster.
#' @param seurat_obj Seurat object containing the clustered data
#' @param active_clusters List of active clusters to include in analysis
#' @param genes_per_cluster Number of top genes to include per cluster
#' @param state Reactive state list
#' @keywords internal
runGeneralHeatmapAnalysis <- function(seurat_obj, active_clusters, genes_per_cluster, state) {
  if (length(active_clusters) == 0) {
    showNotification("No active clusters selected.", type = "warning")
    return(NULL)
  }
  
  # Convert to numeric cluster ids
  active_clusters_num <- as.numeric(active_clusters)
  
  if (length(active_clusters_num) > 0) {
    withProgress(message = "Computing general heatmap...", value = 0, detail = "Preparing analysis", {
      tryCatch({
        # Collect DE genes for each active cluster independently
        all_de_results <- list()
        
        # Get cluster labels
        cluster_labels <- state$cluster_labels
        
        # Create default labels if needed
        if (is.null(cluster_labels)) {
          unique_clusters <- sort(unique(seurat_obj$seurat_clusters))
          cluster_labels <- setNames(
            paste("Cluster", unique_clusters),
            as.character(unique_clusters)
          )
        }
        
        # Prepare cell count info
        num_cells <- ncol(seurat_obj)
        
        # Set initial status message
        state$heatmap_status_message <- paste0(
          "Computing general heatmap with ", length(active_clusters_num), 
          " clusters, ", num_cells, " cells, and ", 
          genes_per_cluster, " top genes per cluster"
        )
        
        for (i in seq_along(active_clusters_num)) {
          cluster <- active_clusters_num[i]
          
          # Get the cluster label safely
          cluster_label <- getClusterLabel(cluster, cluster_labels)
          
          incProgress(1 / length(active_clusters_num), 
                      detail = paste("Finding markers for cluster", cluster_label))
          
          # Find markers for this cluster vs all other active clusters
          res <- FindMarkers(
            seurat_obj,
            ident.1 = cluster,
            ident.2 = setdiff(active_clusters_num, cluster),
            min.pct = 0.25,
            logfc.threshold = 0.25
          )
          
          # Add gene names
          res <- addGeneNames(res, seurat_obj)
          
          # Get top N genes for this cluster
          if (nrow(res) > 0) {
            # Sort by p-value and fold change
            res <- res[order(res$p_val_adj, -abs(res$avg_log2FC)), ]
            top_cluster_genes <- rownames(res)[1:min(genes_per_cluster, nrow(res))]
            all_de_results[[as.character(cluster)]] <- top_cluster_genes
          }
        }
        
        # Store cluster order for later use in visualization
        state$general_heatmap_clusters(names(all_de_results))
        
        incProgress(0.8, detail = "Combining results")
        
        # Combine all unique genes
        unique_genes <- unique(unlist(all_de_results))
        
        if (is.null(unique_genes) || length(unique_genes) == 0) {
          showNotification("No significant genes found for selected clusters.", type = "warning")
          return(NULL)
        }
        
        # Create gene order by cluster
        gene_order <- c()
        for (cluster in names(all_de_results)) {
          gene_order <- c(gene_order, all_de_results[[cluster]])
        }
        # Remove duplicates but keep order
        gene_order <- unique(gene_order)
        
        incProgress(0.9, detail = "Building heatmap")
        
        # Store genes for heatmap
        state$general_heatmap_genes(gene_order)
        state$heatmap_type("general")
        
        # At the end, update the status message with the results
        state$heatmap_status_message <- paste0(
          "Generated heatmap with ", length(active_clusters_num), 
          " clusters, ", num_cells, " cells, and ", 
          length(gene_order), " unique top genes. Using ", 
          genes_per_cluster, " genes per cluster."
        )
        
        incProgress(1.0, detail = "Completed")
      }, error = function(e) {
        showNotification(paste("Error computing general heatmap:", e$message), type = "error")
      })
    })
  } else {
    showNotification("No active clusters selected. Please select at least one cluster.", 
                     type = "warning")
    state$general_heatmap_genes(NULL)
  }
}

#' @title Generate DE Heatmap
#' @description Generate a heatmap visualization for top differentially expressed genes.
#' @param seurat_obj Seurat object containing the data
#' @param top_n_genes Number of top genes to include in the heatmap
#' @param de_results Data frame containing differential expression results
#' @param state Reactive state list
#' @keywords internal
generateDEHeatmap <- function(seurat_obj, top_n_genes, de_results, state) {
  withProgress(message = 'Generating heatmap...', {
    # Ensure we don't request more genes than we have
    n_genes <- min(top_n_genes, nrow(de_results))
    
    if (n_genes == 0) {
      showNotification("No genes available for heatmap.", type = "warning")
      return(NULL)
    }
    
    # Get top genes - sort by p-value and then fold change
    de_results <- de_results[order(de_results$p_val_adj, -abs(de_results$avg_log2FC)), ]
    top_genes <- rownames(de_results)[1:n_genes]
    
    # Determine what grouping to use based on the comparison
    group_by <- "seurat_clusters"  # Default
    
    # Check if this is a sample comparison
    if ("comparison" %in% colnames(de_results)) {
      comparison_text <- unique(de_results$comparison)[1]
      if (grepl("Sample", comparison_text, ignore.case = TRUE)) {
        group_by <- "sample"
      } else if (grepl("Condition", comparison_text, ignore.case = TRUE)) {
        # Find a condition column - look for the one mentioned in the comparison
        condition_cols <- grep("condition|treatment|group|genotype", 
                               colnames(seurat_obj@meta.data), 
                               value = TRUE, ignore.case = TRUE)
        if (length(condition_cols) > 0) {
          # Try to find the best match based on the comparison text
          best_match <- condition_cols[1]  # Default to first
          for (col in condition_cols) {
            col_values <- unique(seurat_obj@meta.data[[col]])
            if (any(sapply(col_values, function(val) grepl(val, comparison_text, fixed = TRUE)))) {
              best_match <- col
              break
            }
          }
          group_by <- best_match
        }
      }
    }
    
    # Store the analysis type
    state$analysis_type <- switch(group_by,
                                  "seurat_clusters" = "cluster_de",
                                  "sample" = "sample_de",
                                  "condition_de")
    
    # Update state
    state$heatmap_data(top_genes)
    state$heatmap_type("specific")
    state$group_by <- group_by  # Store the grouping in state
  })
}

#' @title Render DE Results UI
#' @description Creates the UI elements for displaying differential expression results,
#'   including heatmap and boxplot visualizations.
#' @param ns Namespace function
#' @param state Reactive state list
#' @return UI elements for DE results
#' @keywords internal
renderDEResultsUI <- function(ns, state, clustered_seurat, condition_column = NULL) {
  current_state <- state$analysis_state()
  
  if (is.null(current_state) || current_state == "none") {
    return(NULL)
  }
  
  # List of all DE analysis states that show results in the DE Results UI
  de_results_states <- c("one_vs_all", "pairwise", "sample_comparison", "condition_comparison")
  
  if (current_state %in% de_results_states) {
    current_de_genes <- state$de_genes()
    
    if (is.null(current_de_genes) || nrow(current_de_genes) == 0) {
      return(div(
        class = "alert alert-info",
        "No differential expression results available yet. Run an analysis to see results."
      ))
    }
    
    # Add filter context information if available
    filter_context_ui <- NULL
    if ("filter_context" %in% colnames(current_de_genes)) {
      filter_context <- current_de_genes$filter_context[1]
      if (!is.null(filter_context) && filter_context != "") {
        filter_context_ui <- div(
          class = "alert alert-info",
          style = "margin-bottom: 15px;",
          icon("info-circle"),
          filter_context
        )
      }
    }
    
    has_heatmap <- !is.null(state$heatmap_data()) && 
      !is.null(state$heatmap_type()) && 
      state$heatmap_type() == "specific"
    
    # Create a more descriptive title based on the analysis type
    title_prefix <- switch(current_state,
                           "one_vs_all" = "Cluster",
                           "pairwise" = "Cluster",
                           "sample_comparison" = "Sample",
                           "condition_comparison" = "Condition",
                           "Differential")
    
    tagList(
      h3(paste(title_prefix, "Expression Results")),
      # Insert filter context info if available
      if (!is.null(filter_context_ui)) filter_context_ui,
      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
          h4(style = "margin: 0;", "Volcano Plot"),
          downloadButton(ns("downloadVolcanoData"), "Download Data", 
                         class = "btn-sm btn-info"),
          downloadButton(ns("downloadVolcanoPlot"), "Save Plot", 
                         class = "btn-sm btn-success")
      ),
      plotOutput(ns("volcanoPlot"), height = "400px"),
      div(style = "margin-top: 20px;"),
      DT::dataTableOutput(ns("deTable")),
      div(
        style = "margin-top: 20px;",
        wellPanel(
          h4("Generate Expression Heatmap"),
          fluidRow(
            column(6,
                   numericInput(ns("top_n_genes"), 
                                "Number of top DE genes:", 
                                value = 20, 
                                min = 5, 
                                max = 100)
            ),
            column(6,
                   actionButton(ns("generate_heatmap"), "Generate Heatmap",
                                class = "btn-primary")
            )
          )
        )
      ),
      
      # Only show heatmap save button if heatmap data exists
      div(style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0 10px 0;",
          h4(style = "margin: 0;", "Expression Heatmap"),
          if (has_heatmap) {
            div(
              style = "display: flex; gap: 10px;",
              downloadButton(ns("downloadHeatmapData"), "Download Data", 
                             class = "btn-sm btn-info"),
              downloadButton(ns("downloadHeatmapPlot"), "Save Plot", 
                             class = "btn-sm btn-success")
            )
          }
      ),
      plotOutput(ns("heatmapPlot"), height = "600px"),
      
      createBoxplotUI(ns, current_de_genes, clustered_seurat, state$active_clusters, state$cluster_labels, condition_column)
    )
  } else if (current_state %in% c("general_heatmap", "sample_heatmap", "condition_heatmap")) {
    # For heatmap states, return NULL since that UI is handled by generalHeatmapUI
    return(NULL)
  } else {
    # For unknown states, return a warning
    return(div(
      class = "alert alert-warning",
      paste("Unknown analysis state:", current_state)
    ))
  }
}

#' @title Render General Heatmap UI
#' @description Creates the UI elements for displaying the general heatmap.
#' @param ns Namespace function
#' @param state Reactive state list
#' @param active_cluster_list Reactive expression containing active clusters
#' @return UI elements for general heatmap
#' @keywords internal
renderGeneralHeatmapUI <- function(ns, state, active_cluster_list) {
  current_state <- state$analysis_state()
  current_genes <- state$general_heatmap_genes()
  
  # Update this condition to include all heatmap states
  if (is.null(current_state) || 
      !(current_state %in% c("general_heatmap", "sample_heatmap", "condition_heatmap"))) {
    return(NULL)
  }
  
  if (is.null(current_genes) || length(current_genes) == 0) {
    return(div(
      class = "alert alert-info",
      "No genes found for heatmap. Try adjusting parameters or selecting different groups."
    ))
  }
  
  # Determine which active items to check based on the analysis state
  has_active_items <- FALSE
  analysis_type <- ""
  
  if (current_state == "general_heatmap") {
    current_active <- active_cluster_list()
    has_active_items <- length(current_active) > 0
    analysis_type <- "clusters"
  } else if (current_state == "sample_heatmap") {
    has_active_items <- TRUE  # We already filtered for active samples during the analysis
    analysis_type <- "samples"
  } else if (current_state == "condition_heatmap") {
    has_active_items <- TRUE  # We already filtered for active conditions during the analysis
    analysis_type <- "conditions"
  }
  
  if (!has_active_items) {
    return(div(
      class = "alert alert-warning",
      paste0("No active ", analysis_type, " selected for heatmap visualization.")
    ))
  }
  
  # Create a title based on the analysis type
  heatmap_title <- switch(current_state,
                          "general_heatmap" = "General Cluster Heatmap",
                          "sample_heatmap" = "Sample Comparison Heatmap",
                          "condition_heatmap" = "Condition Comparison Heatmap",
                          "General Heatmap")
  
  # Get status message - properly handle either a reactive or regular variable
  status_message <- NULL
  if (!is.null(state$heatmap_status_message)) {
    if (is.reactive(state$heatmap_status_message)) {
      status_message <- state$heatmap_status_message()
    } else {
      status_message <- state$heatmap_status_message
    }
  }
  
  # Creating UI elements
  ui_elements <- list(
    div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
        h3(style = "margin: 0;", heatmap_title),
        downloadButton(ns("downloadGeneralHeatmapData"), "Download Data", 
                       class = "btn-sm btn-info"),
        downloadButton(ns("downloadGeneralHeatmapPlot"), "Save Plot", 
                       class = "btn-sm btn-success")
    )
  )
  
  # Add status message if available
  if (!is.null(status_message)) {
    ui_elements <- c(ui_elements, 
                     list(div(
                       class = "alert alert-info",
                       style = "margin-top: 10px; margin-bottom: 10px;",
                       icon("info-circle"),
                       HTML(status_message)
                     )))
  }
  
  # Add plot
  ui_elements <- c(ui_elements, 
                   list(plotOutput(ns("generalHeatmapPlot"), height = "800px")))
  
  # Return as tagList
  do.call(tagList, ui_elements)
}

#' @title Setup Heatmap Download Handlers
#' @description Sets up download handlers for heatmap visualizations.
#' @param output Shiny output object
#' @param state Reactive state list
#' @param clustered_seurat Clustered Seurat object
#' @param active_cluster_list Reactive expression containing active clusters
#' @keywords internal
setupHeatmapDownloadHandlers <- function(output, state, clustered_seurat, active_cluster_list) {
  # Download handler for the heatmap
  output$downloadHeatmapPlot <- downloadHandler(
    filename = function() {
      group_by <- if (!is.null(state$group_by)) state$group_by else "seurat_clusters"
      prefix <- switch(group_by,
                       "seurat_clusters" = "cluster",
                       "sample" = "sample",
                       "condition")
      paste0(prefix, "_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      # Create high-resolution PNG device
      png(file, width = 3000, height = 3600, res = 300)
      
      # Generate the heatmap directly
      current_heatmap_data <- state$heatmap_data()
      
      # Determine what grouping to use
      group_by <- if (!is.null(state$group_by)) state$group_by else "seurat_clusters"
      
      if (group_by == "sample" && !is.null(sample_management)) {
        # Use sample management for sample-based heatmap
        active_items <- sample_management$getActiveSampleIds()
        sample_labels <- sample_management$getSampleLabels()
        
        createExpressionHeatmap(
          clustered_seurat(),
          current_heatmap_data,
          sample_labels,
          active_items,
          group_by = "sample"
        )
      } else if (group_by != "seurat_clusters" && group_by != "sample" && 
                 !is.null(condition_management) && 
                 group_by == condition_management$getConditionColumn()) {
        # Use condition management for condition-based heatmap
        active_items <- condition_management$getActiveConditions()
        condition_labels <- condition_management$getConditionLabels()
        
        createExpressionHeatmap(
          clustered_seurat(),
          current_heatmap_data,
          condition_labels,
          active_items,
          group_by = group_by
        )
      } else {
        # Default cluster-based heatmap
        # Get cluster labels
        current_labels <- NULL
        tryCatch({
          if (!is.null(state$cluster_labels)) {
            if (is.function(state$cluster_labels)) {
              current_labels <- state$cluster_labels()
            } else {
              current_labels <- state$cluster_labels
            }
          }
        }, error = function(e) {
          # Continue with null if error
        })
        
        if (is.null(current_labels)) {
          unique_clusters <- sort(unique(clustered_seurat()$seurat_clusters))
          current_labels <- setNames(
            paste("Cluster", unique_clusters),
            as.character(unique_clusters)
          )
        }
        
        active_clusters_num <- as.numeric(active_cluster_list())
        
        createExpressionHeatmap(
          clustered_seurat(),
          current_heatmap_data,
          current_labels,
          active_clusters_num,
          group_by = "seurat_clusters"
        )
      }
      
      # Close the device to save the file
      dev.off()
    }
  )
  
  # General heatmap download
  output$downloadGeneralHeatmapPlot <- downloadHandler(
    filename = function() {
      current_state <- state$analysis_state()
      prefix <- switch(current_state,
                       "general_heatmap" = "cluster",
                       "sample_heatmap" = "sample",
                       "condition_heatmap" = "condition",
                       "general")
      paste0(prefix, "_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      # Create high-resolution PNG device
      png(file, width = 3600, height = 4200, res = 300)
      
      # Generate the heatmap based on analysis type
      current_state <- state$analysis_state()
      
      if (current_state == "general_heatmap") {
        # Cluster-based heatmap (existing code)
        current_labels <- NULL
        tryCatch({
          if (!is.null(state$cluster_labels)) {
            if (is.function(state$cluster_labels)) {
              current_labels <- state$cluster_labels()
            } else {
              current_labels <- state$cluster_labels
            }
          }
        }, error = function(e) {
          # Continue with null if error
        })
        
        if (is.null(current_labels)) {
          unique_clusters <- sort(unique(clustered_seurat()$seurat_clusters))
          current_labels <- setNames(
            paste("Cluster", unique_clusters),
            as.character(unique_clusters)
          )
        }
        
        current_active <- active_cluster_list()
        active_clusters_num <- as.numeric(current_active)
        ordered_genes <- state$general_heatmap_genes()
        
        # Create the heatmap with ordered genes
        createGeneralHeatmap(clustered_seurat(), 
                             ordered_genes, 
                             current_labels, 
                             active_clusters_num,
                             cluster_order = state$general_heatmap_clusters())
      } else if (current_state == "sample_heatmap") {
        # Sample-based heatmap
        sample_labels <- NULL
        if (!is.null(sample_management) && is.function(sample_management$getSampleLabels)) {
          sample_labels <- sample_management$getSampleLabels()
        }
        
        seurat_obj <- clustered_seurat()
        Idents(seurat_obj) <- "sample"
        createGeneralHeatmap(
          seurat_obj,
          state$general_heatmap_genes(),
          sample_labels,
          NULL,
          NULL
        )
      } else if (current_state == "condition_heatmap") {
        # Condition-based heatmap
        condition_column <- NULL
        condition_labels <- NULL
        
        if (!is.null(condition_management)) {
          if (is.function(condition_management$getConditionColumn)) {
            condition_column <- condition_management$getConditionColumn()
          }
          if (is.function(condition_management$getConditionLabels)) {
            condition_labels <- condition_management$getConditionLabels()
          }
        }
        
        seurat_obj <- clustered_seurat()
        if (!is.null(condition_column)) {
          Idents(seurat_obj) <- condition_column
          createGeneralHeatmap(
            seurat_obj,
            state$general_heatmap_genes(),
            condition_labels,
            NULL,
            NULL
          )
        }
      }
      
      # Close the device to save the file
      dev.off()
    }
  )
  
  # Download handler for volcano plot data
  output$downloadVolcanoData <- downloadHandler(
    filename = function() {
      paste("volcano_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
    },
    content = function(file) {
      # Extract the volcano plot data
      req(state$de_genes())
      
      volcano_data <- extractVolcanoData(state$de_genes())
      
      if (!is.null(volcano_data)) {
        # Write to CSV
        write.csv(volcano_data, file, row.names = FALSE)
      } else {
        # Create a simple message if no data available
        write.csv(data.frame(message = "No volcano plot data available"), file, row.names = FALSE)
      }
    }
  )
  
  # Download handler for heatmap data
  output$downloadHeatmapData <- downloadHandler(
    filename = function() {
      paste("heatmap_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip", sep = "")
    },
    content = function(file) {
      req(state$heatmap_data(), state$heatmap_type() == "specific", clustered_seurat())
      
      # Create temp directory for files
      temp_dir <- tempdir()
      raw_file <- file.path(temp_dir, "heatmap_raw_expression.csv")
      scaled_file <- file.path(temp_dir, "heatmap_scaled_expression.csv")
      genes_file <- file.path(temp_dir, "heatmap_genes_list.csv")
      readme_file <- file.path(temp_dir, "README.txt")
      
      # Get cluster labels
      cluster_labels <- if (!is.null(state$cluster_labels)) {
        state$cluster_labels
      } else if (!is.null(cluster_management)) {
        tryCatch({
          cluster_management$getClusterLabels()
        }, error = function(e) {
          NULL
        })
      } else {
        NULL
      }
      
      # Get active clusters
      active_clusters <- tryCatch({
        if (is.function(active_cluster_list)) {
          as.numeric(active_cluster_list())
        } else {
          NULL
        }
      }, error = function(e) {
        NULL
      })
      
      # Determine grouping variable
      group_by <- "seurat_clusters"  # Default
      if (!is.null(state$group_by)) {
        group_by <- state$group_by
      }
      
      # Extract the heatmap data
      heatmap_data <- extractHeatmapData(
        clustered_seurat(),
        state$heatmap_data(),
        group_by = group_by,
        cluster_labels = cluster_labels,
        active_groups = active_clusters
      )
      
      if (!is.null(heatmap_data)) {
        # Write raw expression data
        write.csv(heatmap_data$raw_data, raw_file)
        
        # Write scaled expression data
        write.csv(heatmap_data$scaled_data, scaled_file)
        
        # Write genes list with symbols
        gene_info <- data.frame(
          gene_id = heatmap_data$genes,
          gene_symbol = heatmap_data$gene_symbols
        )
        write.csv(gene_info, genes_file, row.names = FALSE)
        
        # Create README file
        readme_text <- c(
          "Heatmap Data Export",
          "===================",
          "",
          paste("Date:", Sys.time()),
          paste("Group by:", heatmap_data$group_by),
          paste("Number of genes:", length(heatmap_data$genes)),
          paste("Number of groups:", length(heatmap_data$groups)),
          "",
          "Files included:",
          "- heatmap_raw_expression.csv: Raw expression values (averages by group)",
          "- heatmap_scaled_expression.csv: Scaled expression values used in heatmap visualization",
          "- heatmap_genes_list.csv: List of genes with their symbols"
        )
        writeLines(readme_text, readme_file)
        
        # Create ZIP file
        zip(file, c(raw_file, scaled_file, genes_file, readme_file), flags = "-j")
      } else {
        # Create a simple message if no data available
        message_file <- file.path(temp_dir, "message.txt")
        writeLines("No heatmap data available", message_file)
        zip(file, message_file, flags = "-j")
      }
    },
    contentType = "application/zip"
  )
  
  # Download handler for general heatmap data
  output$downloadGeneralHeatmapData <- downloadHandler(
    filename = function() {
      prefix <- switch(state$analysis_state(),
                       "general_heatmap" = "cluster",
                       "sample_heatmap" = "sample",
                       "condition_heatmap" = "condition",
                       "general")
      
      paste(prefix, "_heatmap_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip", sep = "")
    },
    content = function(file) {
      req(state$general_heatmap_genes(), clustered_seurat())
      
      # Create temp directory for files
      temp_dir <- tempdir()
      raw_file <- file.path(temp_dir, "heatmap_raw_expression.csv")
      scaled_file <- file.path(temp_dir, "heatmap_scaled_expression.csv")
      genes_file <- file.path(temp_dir, "heatmap_genes_list.csv")
      readme_file <- file.path(temp_dir, "README.txt")
      
      # Determine grouping and get appropriate labels and active items
      group_by <- "seurat_clusters"
      labels <- NULL
      active_items <- NULL
      current_state <- state$analysis_state()
      
      if (current_state == "general_heatmap") {
        # Cluster-based heatmap
        group_by <- "seurat_clusters"
        labels <- cluster_management$getClusterLabels()
        active_items <- as.numeric(active_cluster_list())
      } else if (current_state == "sample_heatmap") {
        # Sample-based heatmap
        group_by <- "sample"
        if (!is.null(sample_management) && is.function(sample_management$getSampleLabels)) {
          labels <- sample_management$getSampleLabels()
        }
        if (!is.null(sample_management) && is.function(sample_management$getActiveSampleIds)) {
          active_items <- sample_management$getActiveSampleIds()
        }
      } else if (current_state == "condition_heatmap") {
        # Condition-based heatmap
        if (!is.null(condition_management) && is.function(condition_management$getConditionColumn)) {
          group_by <- condition_management$getConditionColumn()
        }
        if (!is.null(condition_management) && is.function(condition_management$getConditionLabels)) {
          labels <- condition_management$getConditionLabels()
        }
        if (!is.null(condition_management) && is.function(condition_management$getActiveConditions)) {
          active_items <- condition_management$getActiveConditions()
        }
      }
      
      # Extract the heatmap data
      heatmap_data <- extractHeatmapData(
        clustered_seurat(),
        state$general_heatmap_genes(),
        group_by = group_by,
        cluster_labels = labels,
        active_groups = active_items
      )
      
      if (!is.null(heatmap_data)) {
        # Write raw expression data
        write.csv(heatmap_data$raw_data, raw_file)
        
        # Write scaled expression data
        write.csv(heatmap_data$scaled_data, scaled_file)
        
        # Write genes list with symbols
        gene_info <- data.frame(
          gene_id = heatmap_data$genes,
          gene_symbol = heatmap_data$gene_symbols
        )
        write.csv(gene_info, genes_file, row.names = FALSE)
        
        # Create README file
        group_type <- switch(current_state,
                             "general_heatmap" = "clusters",
                             "sample_heatmap" = "samples",
                             "condition_heatmap" = "conditions",
                             "groups")
        
        readme_text <- c(
          paste("General", toTitleCase(group_type), "Heatmap Data Export"),
          "===========================================",
          "",
          paste("Date:", Sys.time()),
          paste("Group by:", heatmap_data$group_by),
          paste("Number of genes:", length(heatmap_data$genes)),
          paste("Number of groups:", length(heatmap_data$groups)),
          "",
          "Files included:",
          "- heatmap_raw_expression.csv: Raw expression values (averages by group)",
          "- heatmap_scaled_expression.csv: Scaled expression values used in heatmap visualization",
          "- heatmap_genes_list.csv: List of genes with their symbols"
        )
        writeLines(readme_text, readme_file)
        
        # Create ZIP file
        zip(file, c(raw_file, scaled_file, genes_file, readme_file), flags = "-j")
      } else {
        # Create a simple message if no data available
        message_file <- file.path(temp_dir, "message.txt")
        writeLines("No heatmap data available", message_file)
        zip(file, message_file, flags = "-j")
      }
    },
    contentType = "application/zip"
  )
  
  # Helper function for title case conversion
  toTitleCase <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
}

# Additional helper functions for DE analysis UI components

#' @title Create DE Analysis UI
#' @description Creates the main UI elements for differential expression analysis.
#' @param ns Namespace function
#' @param cluster_choices Named vector of cluster choices for UI elements
#' @param sample_choices Named vector of sample choices (optional)
#' @param condition_choices Named vector of condition choices (optional)
#' @param condition_column Name of the condition column (optional)
#' @return UI elements for DE analysis
#' @keywords internal
createDEAnalysisUI <- function(ns, cluster_choices, sample_choices = NULL, condition_choices = NULL, condition_column = NULL) {
  # Check if we have enough data for different types of analysis
  has_clusters <- length(cluster_choices) > 0
  has_samples <- !is.null(sample_choices) && length(sample_choices) > 0
  has_conditions <- !is.null(condition_choices) && length(condition_choices) > 0
  
  # Create analysis type choices based on what's available
  analysis_choices <- c("Cluster" = "cluster")
  if (has_samples) {
    analysis_choices <- c(analysis_choices, "Sample" = "sample")
  }
  if (has_conditions) {
    analysis_choices <- c(analysis_choices, "Condition" = "condition")
  }
  
  # Analysis mode selector UI
  analysis_selector <- wellPanel(
    fluidRow(
      column(12,
             selectInput(ns("analysisMode"), 
                         "Select analysis type:",
                         choices = analysis_choices,
                         selected = "cluster")
      )
    )
  )
  
  # The cluster-based analysis UI - unified approach
  cluster_ui <- tagList(
    fluidRow(
      column(4,
             wellPanel(
               h4("Cluster Comparison"),
               # First group multi-select
               selectizeInput(ns("targetClusters1"), 
                              "Select first group of clusters:", 
                              choices = cluster_choices,
                              selected = if (length(cluster_choices) > 0) cluster_choices[1] else NULL,
                              multiple = TRUE),
               # Second group multi-select
               selectizeInput(ns("targetClusters2"), 
                              "Select second group of clusters:", 
                              choices = c("All" = "All", cluster_choices),
                              selected = "All",
                              multiple = TRUE),
               # Info text
               helpText("Select multiple clusters by clicking or typing. When 'All' is selected for the second group, all non-selected clusters from the first group will be used."),
               # Run button
               actionButton(ns("runDECluster"), "Run Cluster Comparison", 
                            class = "btn-primary",
                            disabled = !has_clusters)
             )
      ),
      column(4,
             wellPanel(
               h4("General Cluster Map"),
               numericInput(ns("genesPerCluster"),
                            "Top genes per cluster:",
                            value = 5,
                            min = 1,
                            max = 50),
               actionButton(ns("runGeneralHeatmap"), "Generate General Heatmap", 
                            class = "btn-primary",
                            disabled = !has_clusters)
             )
      ),
      column(4,
             wellPanel(
               h4("Cluster Analysis Information"),
               p("Cluster-based DE analysis compares gene expression between different cell clusters."),
               p("Use this to identify marker genes that distinguish cell types or states."),
               p("You can now select multiple clusters for each comparison group."),
               p("Select 'All' as the second group to compare the first group against all other clusters.")
             )
      )
    )
  )
  
  # The sample-based analysis UI
  # The sample-based analysis UI with multi-selection support
  sample_ui <- tagList(
    fluidRow(
      column(4,
             wellPanel(
               h4("Sample Comparison"),
               # First group multi-select
               selectizeInput(ns("targetSamples1"), 
                              "Select first group of samples:", 
                              choices = sample_choices,
                              selected = if (length(sample_choices) > 0) sample_choices[1] else NULL,
                              multiple = TRUE),
               # Second group multi-select  
               selectizeInput(ns("targetSamples2"), 
                              "Select second group of samples:", 
                              choices = c("All" = "All", sample_choices),
                              selected = "All",
                              multiple = TRUE),
               # Info text
               helpText("Select multiple samples by clicking or typing. When 'All' is selected for the second group, all non-selected samples from the first group will be used."),
               actionButton(ns("runDESample"), "Run Sample Comparison", 
                            class = "btn-primary",
                            disabled = !has_samples)
             )
      ),
      column(4,
             wellPanel(
               h4("General Sample Map"),
               numericInput(ns("genesPerSample"),
                            "Top genes per sample:",
                            value = 5,
                            min = 1,
                            max = 50),
               actionButton(ns("runSampleHeatmap"), "Generate Sample Heatmap", 
                            class = "btn-primary",
                            disabled = !has_samples)
             )
      ),
      column(4,
             wellPanel(
               h4("Sample Analysis Information"),
               p("Sample-based DE analysis compares gene expression between different samples, regardless of cell clustering."),
               p("Use this to identify differences between experimental groups or treatments across the entire dataset."),
               p("You can now select multiple samples for each comparison group."),
               p("Select 'All' as the second group to compare the first group against all other samples.")
             )
      )
    )
  )
  
  # The condition-based analysis UI
  condition_ui <- tagList(
    fluidRow(
      column(4,
             wellPanel(
               h4("Condition Comparison"),
               p(paste("Using condition column:", ifelse(!is.null(condition_column), condition_column, "None selected"))),
               # First group multi-select
               selectizeInput(ns("targetConditions1"), 
                              "Select first group of conditions:", 
                              choices = condition_choices,
                              selected = if (length(condition_choices) > 0) condition_choices[1] else NULL,
                              multiple = TRUE),
               # Second group multi-select
               selectizeInput(ns("targetConditions2"), 
                              "Select second group of conditions:", 
                              choices = c("All" = "All", condition_choices),
                              selected = "All",
                              multiple = TRUE),
               # Info text
               helpText("Select multiple conditions by clicking or typing. When 'All' is selected for the second group, all non-selected conditions from the first group will be used."),
               actionButton(ns("runDECondition"), "Run Condition Comparison", 
                            class = "btn-primary",
                            disabled = !has_conditions)
             )
      ),
      column(4,
             wellPanel(
               h4("General Condition Map"),
               numericInput(ns("genesPerCondition"),
                            "Top genes per condition:",
                            value = 5,
                            min = 1,
                            max = 50),
               actionButton(ns("runConditionHeatmap"), "Generate Condition Heatmap", 
                            class = "btn-primary",
                            disabled = !has_conditions)
             )
      ),
      column(4,
             wellPanel(
               h4("Condition Analysis Information"),
               p(paste("Using condition column:", ifelse(!is.null(condition_column), condition_column, "None selected"))),
               p("Condition-based DE analysis compares gene expression between different experimental conditions or metadata groups."),
               p("Use this to identify genes that respond to treatments or biological factors across the entire dataset."),
               p("You can now select multiple conditions for each comparison group."),
               p("Select 'All' as the second group to compare the first group against all other conditions.")
             )
      )
    )
  )
  
  # Main content - the UI changes based on the selected analysis mode
  # This will be controlled by JavaScript
  main_ui <- tagList(
    tags$div(id = ns("cluster_analysis_ui"), cluster_ui),
    tags$div(id = ns("sample_analysis_ui"), style = "display: none;", sample_ui),
    tags$div(id = ns("condition_analysis_ui"), style = "display: none;", condition_ui)
  )
  
  # Add JavaScript to switch between UIs based on dropdown selection
  js_code <- paste0("
    $(document).ready(function() {
      $('#", ns("analysisMode"), "').on('change', function() {
        var mode = $(this).val();
        if (mode === 'cluster') {
          $('#", ns("cluster_analysis_ui"), "').show();
          $('#", ns("sample_analysis_ui"), "').hide();
          $('#", ns("condition_analysis_ui"), "').hide();
        } else if (mode === 'sample') {
          $('#", ns("cluster_analysis_ui"), "').hide();
          $('#", ns("sample_analysis_ui"), "').show();
          $('#", ns("condition_analysis_ui"), "').hide();
        } else if (mode === 'condition') {
          $('#", ns("cluster_analysis_ui"), "').hide();
          $('#", ns("sample_analysis_ui"), "').hide();
          $('#", ns("condition_analysis_ui"), "').show();
        }
      });
    });
  ")
  
  tagList(
    tags$script(HTML(js_code)),
    if (!has_clusters && !has_samples && !has_conditions) {
      # No data available message
      div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        "No active clusters, samples, or conditions available. Please activate at least one item in the management sections."
      )
    } else {
      tagList(
        # Analysis type selector
        analysis_selector,
        # Main content that changes based on selection
        main_ui
      )
    },
    
    # Results container for all analysis types
    div(id = ns("deResults"),
        uiOutput(ns("deResultsUI"))
    ),
    # Container for general heatmap
    div(id = ns("generalHeatmapResults"),
        uiOutput(ns("generalHeatmapUI"))
    )
  )
}

# Visualization functions

#' @title Create Boxplot UI
#' @description Creates UI elements for the gene expression boxplot with searchable gene input
#' @param ns Namespace function
#' @param de_results Differential expression results
#' @param clustered_seurat Seurat object with expression data
#' @param active_clusters List of active clusters
#' @param cluster_labels Named vector of cluster labels
#' @return UI elements for boxplot visualization
#' @keywords internal
createBoxplotUI <- function(ns, de_results, seurat_obj, active_clusters, cluster_labels, condition_column = NULL) {
  # Get top genes from DE results for initial suggestions
  top_genes <- NULL
  gene_choices <- NULL
  
  if (!is.null(de_results) && nrow(de_results) > 0) {
    # Sort by significance and fold change
    sorted_results <- de_results[order(de_results$p_val_adj, -abs(de_results$avg_log2FC)), ]
    
    # Get top genes (up to 50 for initial display)
    top_n <- min(50, nrow(sorted_results))
    top_genes <- rownames(sorted_results)[1:top_n]
    
    # Create choices for selectizeInput with gene symbols as display names
    gene_choices <- setNames(
      top_genes,
      ifelse(!is.na(sorted_results$gene[1:top_n]), 
             paste0(sorted_results$gene[1:top_n], " (", top_genes, ")"),
             top_genes)
    )
  } else {
    # If no DE results, provide some initial suggestions from the full dataset
    if (!is.null(seurat_obj) && !is.null(rownames(seurat_obj))) {
      # Take some random genes as examples
      sample_size <- min(50, length(rownames(seurat_obj)))
      if (sample_size > 0) {
        # Randomly sample genes to show initially
        set.seed(42)  # For reproducibility
        sample_genes <- sample(rownames(seurat_obj), sample_size)
        
        # Try to get gene symbols if available
        gene_mapping <- seurat_obj@misc$gene_mapping
        if (!is.null(gene_mapping)) {
          gene_labels <- sapply(sample_genes, function(gene) {
            if (gene %in% names(gene_mapping) && !is.na(gene_mapping[gene])) {
              paste0(gene_mapping[gene], " (", gene, ")")
            } else {
              gene
            }
          })
          gene_choices <- setNames(sample_genes, gene_labels)
        } else {
          gene_choices <- setNames(sample_genes, sample_genes)
        }
      }
    }
  }
  
  if (is.null(gene_choices) || length(gene_choices) == 0) {
    # Fallback if we still don't have gene choices
    gene_choices <- setNames(
      "No genes available",
      "No genes available"
    )
  }
  
  # Start with basic grouping options
  group_choices <- c("Cluster" = "seurat_clusters", "Sample" = "sample")
  
  # Add condition column if available from condition management
  if (!is.null(condition_column) && condition_column %in% colnames(seurat_obj@meta.data)) {
    # Create a nice display name
    display_name <- gsub("_", " ", condition_column)
    display_name <- paste0(toupper(substr(display_name, 1, 1)), substr(display_name, 2, nchar(display_name)))
    group_choices[display_name] <- condition_column
  }
  
  # Also scan for other possible condition columns in the metadata
  if (!is.null(seurat_obj)) {
    meta_cols <- colnames(seurat_obj@meta.data)
    condition_cols <- grep("condition|treatment|group|genotype|timepoint", 
                           meta_cols, value = TRUE, ignore.case = TRUE)
    
    if (length(condition_cols) > 0) {
      for (col in condition_cols) {
        if (col != "seurat_clusters" && col != "sample" && col != condition_column) {
          # Create a nice display name
          display_name <- gsub("_", " ", col)
          display_name <- paste0(toupper(substr(display_name, 1, 1)), substr(display_name, 2, nchar(display_name)))
          group_choices[display_name] <- col
        }
      }
    }
  }
  
  # Create UI elements
  tagList(
    fluidRow(
      column(12,
             h4("Gene Expression Boxplot"),
             p("Search and select any gene to visualize expression distribution. Type to search among all genes.")
      )
    ),
    fluidRow(
      column(6,
             # Enhanced selectizeInput with server-side search - starting with an empty selection
             selectizeInput(ns("boxplot_gene"), 
                            "Select gene:", 
                            choices = gene_choices,
                            selected = NULL,  # Start with empty selection
                            options = list(
                              placeholder = 'Type to search for any gene',
                              onInitialize = I('function() { this.setValue(""); }'),
                              closeAfterSelect = TRUE,
                              searchField = c('value', 'label'),
                              load = I('function(query, callback) {
                              if (!query.length) return callback();
                              
                              // Send query to server for processing
                              Shiny.setInputValue("gene_search_query", query);
                            }'),
                              render = I('{
                              option: function(item, escape) {
                                return "<div>" + escape(item.label) + "</div>";
                              }
                            }')
                            ))
      ),
      column(6,
             selectInput(ns("boxplot_group"), 
                         "Group by:", 
                         choices = group_choices,
                         selected = "seurat_clusters")
      )
    ),
    fluidRow(
      column(4,
             checkboxInput(ns("boxplot_stats"), 
                           "Show statistical comparisons", 
                           value = TRUE)
      ),
      column(8,
             div(
               style = "display: flex; justify-content: flex-end; gap: 10px;",
               downloadButton(ns("download_boxplot_data"), 
                              "Download Data", 
                              class = "btn-sm btn-info"),
               downloadButton(ns("download_boxplot"), 
                              "Download Plot", 
                              class = "btn-sm btn-success")
             )
      )
    ),
    fluidRow(
      column(12,
             plotOutput(ns("gene_boxplot"), height = "400px")
      )
    )
  )
}

setupGeneSearch <- function(input, output, session, clustered_seurat) {
  ns <- session$ns
  
  # Create an observer for gene search queries
  observeEvent(input$gene_search_query, {
    req(input$gene_search_query, clustered_seurat())
    query <- input$gene_search_query
    
    if (nchar(query) < 2) return() # Require at least 2 characters
    
    # Get the Seurat object
    seurat_obj <- clustered_seurat()
    
    # Get gene mapping if available
    gene_mapping <- seurat_obj@misc$gene_mapping
    
    # Find matches in gene IDs - search across ALL genes in the dataset
    id_matches <- grep(query, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
    
    # Find matches in gene symbols if mapping exists
    symbol_matches <- character(0)
    if (!is.null(gene_mapping)) {
      symbol_match_indices <- grep(query, gene_mapping, ignore.case = TRUE)
      if (length(symbol_match_indices) > 0) {
        symbol_matches <- names(gene_mapping)[symbol_match_indices]
        # Only include genes that are in the dataset
        symbol_matches <- intersect(symbol_matches, rownames(seurat_obj))
      }
    }
    
    # Combine matches, prioritizing symbol matches
    all_matches <- unique(c(symbol_matches, id_matches))
    
    # Limit to first 100 matches if there are too many
    if (length(all_matches) > 100) {
      all_matches <- all_matches[1:100]
    }
    
    # If no matches found
    if (length(all_matches) == 0) {
      # Return an empty result with a message
      updateSelectizeInput(session, "boxplot_gene", 
                           choices = list("No matching genes found" = ""), 
                           selected = "")
      return()
    }
    
    # Create labels for matches
    if (!is.null(gene_mapping)) {
      match_labels <- sapply(all_matches, function(id) {
        if (id %in% names(gene_mapping) && !is.na(gene_mapping[id])) {
          paste0(gene_mapping[id], " (", id, ")")
        } else {
          id
        }
      })
    } else {
      match_labels <- all_matches
    }
    
    # Update the selectize input
    match_choices <- setNames(all_matches, match_labels)
    updateSelectizeInput(session, "boxplot_gene", choices = match_choices, server = TRUE)
  })
}

#' @title Run Sample Heatmap Analysis
#' @description Generates a heatmap of top marker genes for each sample.
#' @param seurat_obj Seurat object containing the data
#' @param active_samples List of active samples to include in analysis
#' @param genes_per_sample Number of top genes to include per sample
#' @param state Reactive state list
#' @keywords internal
runSampleHeatmapAnalysis <- function(seurat_obj, active_samples, genes_per_sample, state) {
  if (length(active_samples) == 0) {
    showNotification("No active samples selected.", type = "warning")
    return(NULL)
  }
  
  withProgress(message = "Computing sample heatmap...", value = 0, detail = "Preparing analysis", {
    tryCatch({
      # Set the identity to 'sample'
      Idents(seurat_obj) <- "sample"
      
      # Collect DE genes for each active sample independently
      all_de_results <- list()
      
      for (i in seq_along(active_samples)) {
        sample <- active_samples[i]
        
        incProgress(1 / length(active_samples), 
                    detail = paste("Finding markers for sample", sample))
        
        # Find markers for this sample vs all other active samples
        res <- FindMarkers(
          seurat_obj,
          ident.1 = sample,
          ident.2 = setdiff(active_samples, sample),
          min.pct = 0.25,
          logfc.threshold = 0.25
        )
        
        # Add gene names
        res <- addGeneNames(res, seurat_obj)
        
        # Get top N genes for this sample
        if (nrow(res) > 0) {
          # Sort by p-value and fold change
          res <- res[order(res$p_val_adj, -abs(res$avg_log2FC)), ]
          top_sample_genes <- rownames(res)[1:min(genes_per_sample, nrow(res))]
          all_de_results[[sample]] <- top_sample_genes
        }
      }
      
      incProgress(0.8, detail = "Combining results")
      
      # Combine all unique genes
      unique_genes <- unique(unlist(all_de_results))
      
      if (is.null(unique_genes) || length(unique_genes) == 0) {
        showNotification("No significant genes found for selected samples.", type = "warning")
        return(NULL)
      }
      
      # Create gene order by sample
      gene_order <- c()
      for (sample in names(all_de_results)) {
        gene_order <- c(gene_order, all_de_results[[sample]])
      }
      # Remove duplicates but keep order
      gene_order <- unique(gene_order)
      
      incProgress(0.9, detail = "Building heatmap")
      
      # Store genes for heatmap and update state
      state$general_heatmap_genes(gene_order)
      state$heatmap_type("general")
      state$analysis_state("sample_heatmap")
      
      incProgress(1.0, detail = "Completed")
    }, error = function(e) {
      showNotification(paste("Error computing sample heatmap:", e$message), type = "error")
    })
  })
}

#' @title Run Condition Heatmap Analysis
#' @description Generates a heatmap of top marker genes for each condition value.
#' @param seurat_obj Seurat object containing the data
#' @param active_conditions List of active conditions to include in analysis
#' @param condition_column Name of the metadata column containing condition values
#' @param genes_per_condition Number of top genes to include per condition
#' @param state Reactive state list
#' @keywords internal
runConditionHeatmapAnalysis <- function(seurat_obj, active_conditions, condition_column, genes_per_condition, state) {
  if (length(active_conditions) == 0) {
    showNotification("No active conditions selected.", type = "warning")
    return(NULL)
  }
  
  if (is.null(condition_column)) {
    showNotification("No condition column selected.", type = "error")
    return(NULL)
  }
  
  withProgress(message = "Computing condition heatmap...", value = 0, detail = "Preparing analysis", {
    tryCatch({
      # Set the identity to the condition column
      Idents(seurat_obj) <- condition_column
      
      # Collect DE genes for each active condition independently
      all_de_results <- list()
      
      for (i in seq_along(active_conditions)) {
        condition <- active_conditions[i]
        
        incProgress(1 / length(active_conditions), 
                    detail = paste("Finding markers for condition", condition))
        
        # Find markers for this condition vs all other active conditions
        res <- FindMarkers(
          seurat_obj,
          ident.1 = condition,
          ident.2 = setdiff(active_conditions, condition),
          min.pct = 0.25,
          logfc.threshold = 0.25
        )
        
        # Add gene names
        res <- addGeneNames(res, seurat_obj)
        
        # Get top N genes for this condition
        if (nrow(res) > 0) {
          # Sort by p-value and fold change
          res <- res[order(res$p_val_adj, -abs(res$avg_log2FC)), ]
          top_condition_genes <- rownames(res)[1:min(genes_per_condition, nrow(res))]
          all_de_results[[condition]] <- top_condition_genes
        }
      }
      
      incProgress(0.8, detail = "Combining results")
      
      # Combine all unique genes
      unique_genes <- unique(unlist(all_de_results))
      
      if (is.null(unique_genes) || length(unique_genes) == 0) {
        showNotification("No significant genes found for selected conditions.", type = "warning")
        return(NULL)
      }
      
      # Create gene order by condition
      gene_order <- c()
      for (condition in names(all_de_results)) {
        gene_order <- c(gene_order, all_de_results[[condition]])
      }
      # Remove duplicates but keep order
      gene_order <- unique(gene_order)
      
      incProgress(0.9, detail = "Building heatmap")
      
      # Store genes for heatmap and update state
      state$general_heatmap_genes(gene_order)
      state$heatmap_type("general")
      state$analysis_state("condition_heatmap")
      
      incProgress(1.0, detail = "Completed")
    }, error = function(e) {
      showNotification(paste("Error computing condition heatmap:", e$message), type = "error")
    })
  })
}