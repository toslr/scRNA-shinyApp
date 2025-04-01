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
        setupAnalysisUI(ns, clustered_seurat, cluster_management, state)
      })
    })
    
    # Render Analysis UI
    output$analysisUI <- renderUI({
      setupAnalysisUI(ns, clustered_seurat, cluster_management, state)
    })
    
    # Setup DE analysis handlers
    setupAnalysisHandlers(input, output, clustered_seurat, cluster_management, state, session)
    
    # Return results, status, and labels
    return(list(
      results = state$de_genes,
      status = state$de_status,
      cluster_labels = reactive({ state$cluster_labels })
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
setupAnalysisUI <- function(ns, clustered_seurat, cluster_management, state) {
  if (is.null(cluster_management)) {
    return(div(
      class = "alert alert-warning",
      "Cluster management not initialized. Please refresh the application."
    ))
  }
  
  active_cluster_ids <- cluster_management$getActiveClusterIds()
  
  # Check if we have any active clusters
  if (is.null(active_cluster_ids) || length(active_cluster_ids) == 0) {
    return(div(
      class = "alert alert-warning",
      "No active clusters selected. Please activate at least one cluster in the Cluster Management section in the sidebar."
    ))
  }
  
  # Get current labels either from state or directly from cluster_management
  current_labels <- if (!is.null(state$cluster_labels)) {
    state$cluster_labels
  } else {
    cluster_management$getClusterLabels()
  }
  
  # Create cluster choices for select inputs
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
  
  # Create analysis UI
  createAnalysisUI(ns, cluster_choices)
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
setupAnalysisHandlers <- function(input, output, clustered_seurat, cluster_management, state, session) {
  ns <- session$ns
  
  observe({
    req(cluster_management)
    state$cluster_labels <- cluster_management$getClusterLabels()
  })
  
  # Helper to clear state when starting a new analysis
  clear_state <- function(new_state) {
    state$analysis_state(new_state)
    
    # Clear visualizations that aren't needed for the new state
    if (new_state == "one_vs_all" || new_state == "pairwise") {
      # Reset general heatmap when switching to DE analysis
      state$general_heatmap_genes(NULL)
      state$general_heatmap_clusters(NULL)
    } else if (new_state == "general_heatmap") {
      # Clear specific heatmap data when switching to general
      state$heatmap_data(NULL)
    }
  }
  
  # Get currently active clusters
  active_cluster_list <- reactive({
    cluster_management$getActiveClusterList()
  })
  
  # Handle One vs All analysis
  observeEvent(input$runDEAll, {
    req(clustered_seurat(), input$targetClusterAll)
    current_active_list <- active_cluster_list()
    
    if (length(current_active_list) == 0) {
      showNotification("No active clusters selected.", type = "warning")
      return(NULL)
    }
    
    clear_state("one_vs_all")
    
    # Run DE analysis
    runOneVsAllAnalysis(
      clustered_seurat(), 
      input$targetClusterAll, 
      current_active_list, 
      state,
      cluster_management
    )
  })
  
  # Handle Pairwise analysis
  observeEvent(input$runDEPair, {
    req(clustered_seurat(), input$targetCluster1, input$targetCluster2)
    
    if (input$targetCluster1 == input$targetCluster2) {
      showNotification("Please select different clusters for comparison.", type = "warning")
      return(NULL)
    }
    
    clear_state("pairwise")
    
    # Run pairwise analysis
    runPairwiseAnalysis(
      clustered_seurat(), 
      input$targetCluster1, 
      input$targetCluster2, 
      state,
      cluster_management
    )
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
    renderDEResultsUI(ns, state)
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
    
    DT::datatable(results,
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
    
    # Now call the function with robust error handling
    tryCatch({
      createExpressionHeatmap(
        clustered_seurat(),
        state$heatmap_data(),
        cluster_labels,
        active_clusters
      )
    }, error = function(e) {
      # Return a simple error message plot
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Error creating heatmap:", e$message)) + 
        theme_void()
    })
  })
  
  # General heatmap plot with comprehensive error handling
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
    
    # Get active clusters
    current_active <- tryCatch({
      if (is.function(active_cluster_list)) {
        active_cluster_list()
      } else {
        NULL
      }
    }, error = function(e) {
      NULL
    })
    
    if (is.null(current_active) || length(current_active) == 0) {
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
        as.numeric(current_active),
        state$general_heatmap_clusters()
      )
    }, error = function(e) {
      # Return a simple error message plot
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Error creating general heatmap:", e$message)) + 
        theme_void()
    })
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
}

#' @title Run One vs All DE Analysis
#' @description Performs differential expression analysis comparing one cluster against all others.
#' @param seurat_obj Seurat object containing the clustered data
#' @param target_cluster Target cluster to compare against others
#' @param active_clusters List of active clusters to include in analysis
#' @param state Reactive state list
#' @param cluster_management Cluster management module instance
#' @keywords internal
runOneVsAllAnalysis <- function(seurat_obj, target_cluster, active_clusters, state, cluster_management) {
  withProgress(message = 'Computing one vs all differential expression...', {
    # Subset to active clusters
    active_cells <- seurat_obj$seurat_clusters %in% active_clusters
    seurat_subset <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
    
    # Check if enough active clusters
    if (length(unique(seurat_subset$seurat_clusters)) < 2) {
      showNotification("Need at least two active clusters for comparison.", type = "warning")
      return(NULL)
    }
    
    # Run DE analysis
    de_results <- performDEanalysis(seurat_subset, target_cluster)
    
    # Check results
    if (is.null(de_results) || nrow(de_results) == 0) {
      showNotification("No differential expression results found.", type = "warning")
      return(NULL)
    }
    
    # Add comparison label - use the cluster labels from state first, then fallback to getting them directly
    current_labels <- if (!is.null(state$cluster_labels)) {
      state$cluster_labels
    } else if (!is.null(cluster_management)) {
      cluster_management$getClusterLabels()
    } else {
      NULL
    }
    
    target_cluster_str <- as.character(target_cluster)
    cluster_label <- if (!is.null(current_labels) && target_cluster_str %in% names(current_labels)) {
      current_labels[[target_cluster_str]]
    } else {
      paste("Cluster", target_cluster)
    }
    
    de_results$comparison <- paste(cluster_label, "vs All Active")
    
    # Update state
    state$de_genes(de_results)
    state$de_status("completed")
    state$heatmap_data(NULL)
    state$heatmap_type(NULL)
  })
}

#' @title Run Pairwise DE Analysis
#' @description Performs differential expression analysis comparing two specific clusters.
#' @param seurat_obj Seurat object containing the clustered data
#' @param cluster1 First cluster for comparison
#' @param cluster2 Second cluster for comparison
#' @param state Reactive state list
#' @param cluster_management Cluster management module instance
#' @keywords internal
runPairwiseAnalysis <- function(seurat_obj, cluster1, cluster2, state, cluster_management) {
  withProgress(message = 'Computing pairwise differential expression...', {
    de_results <- performDEanalysis(
      seurat_obj, 
      cluster1, 
      cluster2
    )
    
    # Check results
    if (is.null(de_results) || nrow(de_results) == 0) {
      showNotification("No differential expression results found.", type = "warning")
      return(NULL)
    }
    
    # Add comparison label - use the cluster labels from state first, then fallback to getting them directly
    current_labels <- if (!is.null(state$cluster_labels)) {
      state$cluster_labels
    } else if (!is.null(cluster_management)) {
      cluster_management$getClusterLabels()
    } else {
      NULL
    }
    
    cluster1_str <- as.character(cluster1)
    cluster2_str <- as.character(cluster2)
    
    cluster1_label <- if (!is.null(current_labels) && cluster1_str %in% names(current_labels)) {
      current_labels[[cluster1_str]]
    } else {
      paste("Cluster", cluster1)
    }
    
    cluster2_label <- if (!is.null(current_labels) && cluster2_str %in% names(current_labels)) {
      current_labels[[cluster2_str]]
    } else {
      paste("Cluster", cluster2)
    }
    
    de_results$comparison <- paste(cluster1_label, "vs", cluster2_label)
    
    # Update state
    state$de_genes(de_results)
    state$de_status("completed")
    state$heatmap_data(NULL)
    state$heatmap_type(NULL)
  })
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
    
    # Update state
    state$heatmap_data(top_genes)
    state$heatmap_type("specific")
  })
}

#' @title Render DE Results UI
#' @description Creates the UI elements for displaying differential expression results.
#' @param ns Namespace function
#' @param state Reactive state list
#' @return UI elements for DE results
#' @keywords internal
renderDEResultsUI <- function(ns, state) {
  current_state <- state$analysis_state()
  
  if (is.null(current_state) || current_state == "none") {
    return(NULL)
  }
  
  if (current_state %in% c("one_vs_all", "pairwise")) {
    current_de_genes <- state$de_genes()
    
    if (is.null(current_de_genes) || nrow(current_de_genes) == 0) {
      return(div(
        class = "alert alert-info",
        "No differential expression results available yet. Run an analysis to see results."
      ))
    }
    
    has_heatmap <- !is.null(state$heatmap_data()) && 
      !is.null(state$heatmap_type()) && 
      state$heatmap_type() == "specific"
    
    tagList(
      h3("Differential Expression Results"),
      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
          h4(style = "margin: 0;", "Volcano Plot"),
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
            downloadButton(ns("downloadHeatmapPlot"), "Save Plot", 
                           class = "btn-sm btn-success")
          }
      ),
      plotOutput(ns("heatmapPlot"), height = "600px")
    )
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
  
  if (is.null(current_state) || current_state != "general_heatmap") {
    return(NULL)
  }
  
  if (is.null(current_genes) || length(current_genes) == 0) {
    return(div(
      class = "alert alert-info",
      "No genes found for heatmap. Try adjusting parameters or selecting different clusters."
    ))
  }
  
  current_active <- active_cluster_list()
  has_active_clusters <- length(current_active) > 0
  
  if (!has_active_clusters) {
    return(div(
      class = "alert alert-warning",
      "No active clusters selected for heatmap visualization."
    ))
  }
  
  tagList(
    div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
        h3(style = "margin: 0;", "General Cluster Heatmap"),
        downloadButton(ns("downloadGeneralHeatmapPlot"), "Save Plot", 
                       class = "btn-sm btn-success")
    ),
    plotOutput(ns("generalHeatmapPlot"), height = "800px")
  )
}

#' @title Setup Heatmap Download Handlers
#' @description Sets up download handlers for heatmap visualizations.
#' @param output Shiny output object
#' @param state Reactive state list
#' @param clustered_seurat Clustered Seurat object
#' @param active_cluster_list Reactive expression containing active clusters
#' @keywords internal
setupHeatmapDownloadHandlers <- function(output, state, clustered_seurat, active_cluster_list) {
  # Specific DE heatmap download
  output$downloadHeatmapPlot <- downloadHandler(
    filename = function() {
      paste("heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      # Create high-resolution PNG device
      png(file, width = 3000, height = 3600, res = 300)
      
      # Generate the heatmap directly
      current_heatmap_data <- state$heatmap_data()
      
      # Get cluster labels - try to get from state or use defaults
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
      
      # For all analysis types, show expression in all active clusters
      createExpressionHeatmap(clustered_seurat(), 
                              current_heatmap_data, 
                              current_labels,
                              active_clusters_num)
      
      # Close the device to save the file
      dev.off()
    }
  )
  
  # General heatmap download
  output$downloadGeneralHeatmapPlot <- downloadHandler(
    filename = function() {
      paste("general_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      # Create high-resolution PNG device
      png(file, width = 3600, height = 4200, res = 300)
      
      # Generate the heatmap directly
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
      
      # Close the device to save the file
      dev.off()
    }
  )
}

# Additional helper functions for DE analysis UI components

#' @title Create Analysis UI
#' @description Creates the main UI elements for differential expression analysis.
#' @param ns Namespace function
#' @param cluster_choices Named vector of cluster choices for UI elements
#' @return UI elements for DE analysis
#' @keywords internal
createAnalysisUI <- function(ns, cluster_choices) {
  # Check if we have enough clusters for analysis
  has_clusters <- length(cluster_choices) > 0
  has_multiple_clusters <- length(cluster_choices) > 1
  
  tagList(
    if (!has_clusters) {
      # No active clusters message
      div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        "No active clusters available. Please activate at least one cluster in the Cluster Management section above."
      )
    } else {
      # Main analysis UI when clusters are available
      fluidRow(
        column(4,
               wellPanel(
                 h4("One vs All Analysis"),
                 selectInput(ns("targetClusterAll"), 
                             "Select cluster to compare against all others:", 
                             choices = cluster_choices,
                             selected = if (length(cluster_choices) > 0) cluster_choices[1] else NULL),
                 actionButton(ns("runDEAll"), "Run One vs All DE", 
                              class = "btn-primary",
                              disabled = !has_clusters)
               )
        ),
        column(4,
               wellPanel(
                 h4("One vs One Analysis"),
                 selectInput(ns("targetCluster1"), 
                             "Select first cluster:", 
                             choices = cluster_choices,
                             selected = if (length(cluster_choices) > 0) cluster_choices[1] else NULL),
                 selectInput(ns("targetCluster2"), 
                             "Select second cluster:", 
                             choices = cluster_choices,
                             selected = if (length(cluster_choices) > 1) cluster_choices[2] else cluster_choices[1]),
                 actionButton(ns("runDEPair"), "Run Pairwise DE", 
                              class = "btn-primary",
                              disabled = !has_multiple_clusters)
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
        )
      )
    },
    
    # Results container for One vs All and One vs One
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

#' @title Create Volcano Plot
#' @description Creates a volcano plot from differential expression results.
#' @param results Differential expression results data frame
#' @return A ggplot object with the volcano plot
#' @keywords internal
createVolcanoPlot <- function(results) {
  ggplot(results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = abs(avg_log2FC) > 0.25 & p_val_adj < 0.05)) +
    scale_color_manual(values = c("grey", "red")) +
    theme_classic() +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = unique(results$comparison),
         color = "Significant") +
    theme(legend.position = "bottom")
}

#' @title Create Expression Heatmap
#' @description Creates a heatmap of gene expression across clusters.
#' @param seurat_obj Seurat object containing the data
#' @param top_genes Vector of genes to include in the heatmap
#' @param cluster_labels Named vector of cluster labels
#' @param active_clusters Vector of active cluster IDs
#' @return A pheatmap object showing gene expression
#' @keywords internal
createExpressionHeatmap <- function(seurat_obj, top_genes, cluster_labels = NULL, active_clusters = NULL) {
  # Safety check for empty gene list
  if (length(top_genes) == 0) {
    # Create an empty plot with a message
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No genes to display") + 
             theme_void())
  }
  
  # Ensure cluster_labels is not a function (reactive value)
  if (is.function(cluster_labels)) {
    tryCatch({
      cluster_labels <- cluster_labels()
    }, error = function(e) {
      cluster_labels <- NULL
    })
  }
  
  # Create default cluster labels if not provided
  if (is.null(cluster_labels)) {
    unique_clusters <- sort(unique(seurat_obj$seurat_clusters))
    cluster_labels <- setNames(
      paste("Cluster", unique_clusters),
      as.character(unique_clusters)
    )
  }
  
  # Subset to active clusters if provided
  if (!is.null(active_clusters) && length(active_clusters) > 0) {
    active_cells <- seurat_obj$seurat_clusters %in% active_clusters
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
  }
  
  # Ensure all top genes are in the Seurat object
  genes_present <- intersect(top_genes, rownames(seurat_obj))
  if (length(genes_present) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "None of the selected genes are present in the dataset") + 
             theme_void())
  }
  
  # Get expression data for genes that exist in the dataset
  expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_present, , drop = FALSE]
  
  # Calculate cluster means
  clusters <- seurat_obj$seurat_clusters
  unique_clusters <- sort(unique(clusters))
  
  # Check if we have clusters
  if (length(unique_clusters) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No clusters to display") + 
             theme_void())
  }
  
  cluster_means <- sapply(unique_clusters, function(clust) {
    cluster_cells <- which(clusters == clust)
    if (length(cluster_cells) > 0) {
      rowMeans(expr_data[, cluster_cells, drop = FALSE])
    } else {
      rep(NA, nrow(expr_data))
    }
  })
  
  # Set column names using cluster labels
  col_names <- sapply(unique_clusters, function(x) {
    cluster_key <- as.character(x)
    if (cluster_key %in% names(cluster_labels)) {
      cluster_labels[[cluster_key]]
    } else {
      paste("Cluster", x)
    }
  })
  colnames(cluster_means) <- col_names
  
  # Scale data
  if (is.matrix(cluster_means) && nrow(cluster_means) > 1 && ncol(cluster_means) > 1) {
    # Only scale if we have enough data
    scaled_data <- t(scale(t(cluster_means)))
  } else {
    # Otherwise just use the original data
    scaled_data <- cluster_means
  }
  
  # Get gene labels
  gene_mapping <- seurat_obj@misc$gene_mapping
  gene_labels <- if (!is.null(gene_mapping) && all(rownames(scaled_data) %in% names(gene_mapping))) {
    gene_mapping[rownames(scaled_data)]
  } else {
    rownames(scaled_data)
  }
  gene_labels[is.na(gene_labels)] <- rownames(scaled_data)[is.na(gene_labels)]
  
  # Create heatmap
  tryCatch({
    pheatmap(scaled_data,
             labels_row = gene_labels,
             labels_col = colnames(cluster_means),
             main = paste("Expression Heatmap:", length(genes_present), "Genes"),
             angle_col = 45,
             fontsize_row = 10,
             cluster_cols = FALSE,
             treeheight_row = 0,
             draw = TRUE)
  }, error = function(e) {
    # If pheatmap fails, return a ggplot error message
    ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Error creating heatmap:", e$message)) + 
      theme_void()
  })
}

#' @title Create General Heatmap
#' @description Creates a heatmap of cluster-specific marker genes.
#' @param seurat_obj Seurat object containing the data
#' @param genes Vector of genes to include in the heatmap
#' @param cluster_labels Named vector of cluster labels
#' @param active_clusters Vector of active cluster IDs
#' @param cluster_order Optional vector specifying the order of clusters
#' @return A pheatmap object showing cluster-specific gene expression
#' @keywords internal
createGeneralHeatmap <- function(seurat_obj, genes, cluster_labels = NULL, active_clusters = NULL, cluster_order = NULL) {
  # Safety check for empty gene list
  if (length(genes) == 0) {
    # Create an empty plot with a message
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No genes to display") + 
             theme_void())
  }
  
  # Ensure cluster_labels is not a function (reactive value)
  if (is.function(cluster_labels)) {
    tryCatch({
      cluster_labels <- cluster_labels()
    }, error = function(e) {
      cluster_labels <- NULL
    })
  }
  
  # Create default cluster labels if not provided
  if (is.null(cluster_labels)) {
    unique_clusters <- sort(unique(seurat_obj$seurat_clusters))
    cluster_labels <- setNames(
      paste("Cluster", unique_clusters),
      as.character(unique_clusters)
    )
  }
  
  # Subset to active clusters if provided
  if (!is.null(active_clusters) && length(active_clusters) > 0) {
    existing_active_clusters <- intersect(active_clusters, unique(seurat_obj$seurat_clusters))
    
    if (length(existing_active_clusters) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "None of the active clusters exist in the current dataset") + 
               theme_void())
    }
    
    # Only use active clusters that actually exist
    active_cells <- seurat_obj$seurat_clusters %in% existing_active_clusters
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[active_cells])
  }
  
  # Ensure all genes are in the Seurat object
  genes_present <- intersect(genes, rownames(seurat_obj))
  if (length(genes_present) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "None of the selected genes are present in the dataset") + 
             theme_void())
  }
  
  # Get expression data for selected genes
  expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_present, , drop = FALSE]
  
  # Calculate cluster means
  clusters <- seurat_obj$seurat_clusters
  unique_clusters <- sort(unique(clusters))
  
  # Check if we have clusters
  if (length(unique_clusters) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No clusters to display") + 
             theme_void())
  }
  
  # Order clusters if provided
  if (!is.null(cluster_order) && length(cluster_order) > 0) {
    # Convert to numeric
    ordered_clusters <- as.numeric(cluster_order)
    # Only use clusters that exist in our data
    ordered_clusters <- intersect(ordered_clusters, unique_clusters)
    # If no matching clusters (completely new set), use the default ordering
    if (length(ordered_clusters) == 0) {
      ordered_clusters <- unique_clusters
    } else {
      # Add any clusters that weren't included in the order
      ordered_clusters <- c(ordered_clusters, setdiff(unique_clusters, ordered_clusters))
    }
    
    # Replace unique_clusters with ordered version
    unique_clusters <- ordered_clusters
  }
  
  cluster_means <- sapply(unique_clusters, function(clust) {
    cluster_cells <- which(clusters == clust)
    if (length(cluster_cells) > 0) {
      rowMeans(expr_data[, cluster_cells, drop = FALSE])
    } else {
      rep(NA, nrow(expr_data))
    }
  })
  
  # Set column names using cluster labels
  col_names <- sapply(unique_clusters, function(x) {
    cluster_key <- as.character(x)
    if (cluster_key %in% names(cluster_labels)) {
      cluster_labels[[cluster_key]]
    } else {
      paste("Cluster", x)
    }
  })
  colnames(cluster_means) <- col_names
  
  # Scale data
  if (is.matrix(cluster_means) && nrow(cluster_means) > 1 && ncol(cluster_means) > 1) {
    # Only scale if we have enough data
    scaled_data <- t(scale(t(cluster_means)))
  } else {
    # Otherwise just use the original data
    scaled_data <- cluster_means
  }
  
  # Get gene labels
  gene_mapping <- seurat_obj@misc$gene_mapping
  gene_labels <- if (!is.null(gene_mapping) && all(rownames(scaled_data) %in% names(gene_mapping))) {
    gene_mapping[rownames(scaled_data)]
  } else {
    rownames(scaled_data)
  }
  gene_labels[is.na(gene_labels)] <- rownames(scaled_data)[is.na(gene_labels)]
  
  # Create heatmap with ordered rows (genes already ordered by cluster)
  tryCatch({
    pheatmap(scaled_data,
             labels_row = gene_labels,
             labels_col = colnames(cluster_means),
             main = "Top Cluster-Specific Genes",
             angle_col = 45,
             fontsize_row = 8,
             cluster_cols = FALSE,
             cluster_rows = FALSE,  # Don't cluster rows to maintain diagonal pattern
             treeheight_row = 0,
             draw = TRUE)
  }, error = function(e) {
    # If pheatmap fails, return a ggplot error message
    ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Error creating heatmap:", e$message)) + 
      theme_void()
  })
}