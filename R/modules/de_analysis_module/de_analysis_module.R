# R/modules/de_analysis_module/de_analysis_module.R

deAnalysisUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Cluster management
    wellPanel(
      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
          h4(style = "margin: 0;", "Cluster Management"),
          div(
            actionButton(ns("selectAllClusters"), "Select All", 
                         class = "btn-sm btn-success", 
                         style = "margin-right: 10px;"),
            actionButton(ns("deselectAllClusters"), "Deselect All", 
                         class = "btn-sm btn-danger",
                         style = "margin-right: 10px;"),
            actionButton(ns("updateAllLabels"), "Update All Cluster Labels",
                         class = "btn-primary")
          )
      ),
      div(id = ns("clusterManagement"),
          style = "max-height: 300px; overflow-y: auto;",
          uiOutput(ns("clusterControls"))
      )
    ),
    # Analysis UI
    uiOutput(ns("analysisUI"))
  )
}

deAnalysisServer <- function(id, clustered_seurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize state management
    state <- setupReactiveState()
    
    # Initialize cluster labels and active status when clusters change
    initializeClusterState(clustered_seurat, state)
    
    # Setup UI for cluster controls
    output$clusterControls <- renderUI({
      setupClusterControlsUI(ns, clustered_seurat, state)
    })
    
    # Handle select/deselect all buttons
    setupSelectionHandlers(input, clustered_seurat, state)
    
    # Handle label updates
    setupLabelHandlers(input, clustered_seurat, state)
    
    # Render Analysis UI
    output$analysisUI <- renderUI({
      setupAnalysisUI(ns, clustered_seurat, state)
    })
    
    # Setup DE analysis handlers
    setupAnalysisHandlers(input, output, clustered_seurat, state, session)
    
    # Return results, status, and labels
    return(list(
      results = state$de_genes,
      status = state$de_status,
      labels = state$cluster_labels,
      active = state$active_clusters
    ))
  })
}

# Module state initialization function
setupReactiveState <- function() {
  # Initialize state
  list(
    # Cluster labels and status
    cluster_labels = reactiveVal(NULL),
    active_clusters = reactiveVal(NULL),
    temp_labels = reactiveVal(NULL),
    
    # Heatmap state
    heatmap_data = reactiveVal(NULL),
    heatmap_type = reactiveVal(NULL),
    
    # DE analysis state
    de_genes = reactiveVal(NULL),
    de_status = reactiveVal(NULL),
    general_heatmap_genes = reactiveVal(NULL),
    general_heatmap_clusters = reactiveVal(NULL),
    analysis_state = reactiveVal("none"),
    
    # UI state trackers
    label_inputs = reactiveValues()
  )
}

# Initialize cluster state based on Seurat object
initializeClusterState <- function(clustered_seurat, state) {
  # Get clusters from Seurat object
  observe({
    available_clusters <- getAvailableClusters(clustered_seurat())
    
    # Skip if no clusters available
    if (is.null(available_clusters) || length(available_clusters) == 0) {
      return(NULL)
    }
    
    # Initialize or update labels if needed
    current_labels <- state$cluster_labels()
    if (is.null(current_labels) || 
        !all(as.character(available_clusters) %in% names(current_labels))) {
      
      new_labels <- initializeClusterLabels(available_clusters, current_labels)
      
      state$cluster_labels(new_labels)
      state$temp_labels(new_labels)
    }
    
    # Initialize active status if needed
    current_active <- state$active_clusters()
    if (is.null(current_active) || 
        !all(as.character(available_clusters) %in% names(current_active))) {
      state$active_clusters(initializeActiveStatus(available_clusters, current_active))
    }
  })
  
  # Handle initial population when clusters change
  observeEvent(clustered_seurat(), {
    # Reset heatmap data when Seurat object changes
    state$heatmap_data(NULL)
    state$heatmap_type(NULL)
    state$general_heatmap_genes(NULL)
    state$general_heatmap_clusters(NULL)
    state$analysis_state("none")
    
    available_clusters <- getAvailableClusters(clustered_seurat())
    
    if (!is.null(available_clusters) && length(available_clusters) > 0) {
      # Initialize new labels and active status
      new_labels <- initializeClusterLabels(available_clusters)
      state$cluster_labels(new_labels)
      state$temp_labels(new_labels)
      
      new_active <- initializeActiveStatus(available_clusters)
      state$active_clusters(new_active)
    }
  }, ignoreInit = TRUE)
}

# Get available clusters from Seurat object
getAvailableClusters <- function(seurat_obj) {
  if (is.null(seurat_obj) || !"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }
  sort(unique(seurat_obj$seurat_clusters))
}

# Setup UI for cluster controls
setupClusterControlsUI <- function(ns, clustered_seurat, state) {
  available_clusters <- getAvailableClusters(clustered_seurat())
  
  # Show message if no clusters
  if (is.null(available_clusters) || length(available_clusters) == 0) {
    return(div(
      class = "alert alert-info",
      "No clusters available. Please run the clustering step first."
    ))
  }
  
  # Get current labels and active status
  current_temp_labels <- state$temp_labels()
  current_active <- state$active_clusters()
  
  # Return UI
  createClusterControls(ns, available_clusters, current_temp_labels, current_active)
}

# Setup selection handlers for select/deselect all buttons
setupSelectionHandlers <- function(input, clustered_seurat, state) {
  # Handle Select All button
  observeEvent(input$selectAllClusters, {
    available_clusters <- getAvailableClusters(clustered_seurat())
    if (is.null(available_clusters) || length(available_clusters) == 0) {
      return(NULL)
    }
    
    current_active <- state$active_clusters()
    for (cluster in available_clusters) {
      cluster_key <- as.character(cluster)
      current_active[cluster_key] <- TRUE
    }
    state$active_clusters(current_active)
    
    # Clear general heatmap if it exists
    if (!is.null(state$heatmap_type()) && state$heatmap_type() == "general") {
      state$general_heatmap_genes(NULL)
    }
  })
  
  # Handle Deselect All button
  observeEvent(input$deselectAllClusters, {
    available_clusters <- getAvailableClusters(clustered_seurat())
    if (is.null(available_clusters) || length(available_clusters) == 0) {
      return(NULL)
    }
    
    current_active <- state$active_clusters()
    for (cluster in available_clusters) {
      cluster_key <- as.character(cluster)
      current_active[cluster_key] <- FALSE
    }
    state$active_clusters(current_active)
    
    # Clear general heatmap if it exists
    if (!is.null(state$heatmap_type()) && state$heatmap_type() == "general") {
      state$general_heatmap_genes(NULL)
    }
  })
  
  # Handle active status updates for individual clusters
  observe({
    available_clusters <- getAvailableClusters(clustered_seurat())
    
    # Skip empty clusters
    if (is.null(available_clusters) || length(available_clusters) == 0) {
      return(NULL)
    }
    
    lapply(available_clusters, function(cluster) {
      input_id <- paste0("active_", cluster)
      
      # Only observe if this input exists
      if (!is.null(input[[input_id]])) {
        observeEvent(input[[input_id]], {
          current_active <- state$active_clusters()
          if (is.null(current_active)) {
            return(NULL)
          }
          
          cluster_key <- as.character(cluster)
          
          # Check if this cluster exists in our active_clusters
          if (cluster_key %in% names(current_active)) {
            current_active[cluster_key] <- input[[input_id]]
            state$active_clusters(current_active)
            
            # Reset heatmap data when active clusters change
            if (!is.null(state$heatmap_type()) && state$heatmap_type() == "general") {
              # Only clear general heatmap as it's affected by active clusters
              state$general_heatmap_genes(NULL)
            }
          }
        }, ignoreInit = TRUE)
      }
    })
  })
}

# Setup handlers for label updates
setupLabelHandlers <- function(input, clustered_seurat, state) {
  # Update label_inputs when text changes
  observe({
    available_clusters <- getAvailableClusters(clustered_seurat())
    if (is.null(available_clusters) || length(available_clusters) == 0) {
      return(NULL)
    }
    
    for(cluster in available_clusters) {
      local({
        local_cluster <- cluster
        input_id <- paste0("label_", local_cluster)
        
        # Only observe if this input exists
        if (!is.null(input[[input_id]])) {
          observeEvent(input[[input_id]], {
            state$label_inputs[[input_id]] <- input[[input_id]]
          }, ignoreInit = TRUE)
        }
      })
    }
  })
  
  # Update temp_labels only when update button is clicked
  observeEvent(input$updateAllLabels, {
    available_clusters <- getAvailableClusters(clustered_seurat())
    if (is.null(available_clusters) || length(available_clusters) == 0) {
      return(NULL)
    }
    
    current_temp <- state$temp_labels()
    
    for(cluster in available_clusters) {
      input_id <- paste0("label_", cluster)
      if(!is.null(state$label_inputs[[input_id]])) {
        current_temp[as.character(cluster)] <- state$label_inputs[[input_id]]
      }
    }
    
    state$temp_labels(current_temp)
    state$cluster_labels(current_temp)
  })
}

# Setup Analysis UI
setupAnalysisUI <- function(ns, clustered_seurat, state) {
  available_clusters <- getAvailableClusters(clustered_seurat())
  if (is.null(available_clusters) || length(available_clusters) == 0) {
    return(div(
      class = "alert alert-info",
      "No clusters available. Please run the clustering step first."
    ))
  }
  
  # Get current labels and active status
  current_labels <- state$cluster_labels()
  current_active <- state$active_clusters()
  
  if (is.null(current_labels) || is.null(current_active)) {
    return(div(
      class = "alert alert-warning",
      "Cluster information is being loaded. Please wait..."
    ))
  }
  
  # Get active cluster list
  active_cluster_ids <- as.numeric(names(current_active[current_active == TRUE]))
  
  # Check if we have any active clusters
  if (length(active_cluster_ids) == 0) {
    return(div(
      class = "alert alert-warning",
      "No active clusters selected. Please activate at least one cluster in the Cluster Management section above."
    ))
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

# Setup analysis handlers and results UI
setupAnalysisHandlers <- function(input, output, clustered_seurat, state, session) {
  ns <- session$ns
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
    current_active <- state$active_clusters()
    if (is.null(current_active) || length(current_active) == 0) {
      return(character(0))
    }
    names(current_active[current_active == TRUE])
  })
  
  # Handle One vs All analysis
  observeEvent(input$runDEAll, {
    req(clustered_seurat(), input$targetClusterAll)
    current_active_list <- active_cluster_list()
    
    if (length(current_active_list) == 0) {
      showNotification("No active clusters selected.", type = "warning")
      return(NULL)
    }
    
    print("Running one vs all analysis")
    clear_state("one_vs_all")
    
    # Run DE analysis
    runOneVsAllAnalysis(
      clustered_seurat(), 
      input$targetClusterAll, 
      current_active_list, 
      state
    )
  })
  
  # Handle Pairwise analysis
  observeEvent(input$runDEPair, {
    req(clustered_seurat(), input$targetCluster1, input$targetCluster2)
    
    if (input$targetCluster1 == input$targetCluster2) {
      showNotification("Please select different clusters for comparison.", type = "warning")
      return(NULL)
    }
    
    print("Running pairwise analysis")
    clear_state("pairwise")
    
    # Run pairwise analysis
    runPairwiseAnalysis(
      clustered_seurat(), 
      input$targetCluster1, 
      input$targetCluster2, 
      state
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
    
    print("Running general heatmap analysis")
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
  
  # Specific heatmap plot
  specific_heatmap_plot <- reactive({
    req(state$heatmap_data())
    req(state$heatmap_type() == "specific")
    req(clustered_seurat())
    
    createExpressionHeatmap(
      clustered_seurat(),
      state$heatmap_data(),
      state$cluster_labels(),
      as.numeric(active_cluster_list())
    )
  })
  
  # General heatmap plot
  general_heatmap_plot <- reactive({
    req(state$general_heatmap_genes())
    req(clustered_seurat())
    
    current_active <- active_cluster_list()
    if (length(current_active) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No active clusters selected") + 
               theme_void())
    }
    
    createGeneralHeatmap(
      clustered_seurat(),
      state$general_heatmap_genes(),
      state$cluster_labels(),
      as.numeric(current_active),
      state$general_heatmap_clusters()
    )
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

# Run One vs All DE analysis
runOneVsAllAnalysis <- function(seurat_obj, target_cluster, active_clusters, state) {
  withProgress(message = 'Computing one vs all differential expression...', {
    # Subset to active clusters - fix: seurat_obj is an object, not a function
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
    
    # Add comparison label
    cluster_label <- getClusterLabel(target_cluster, state$cluster_labels())
    de_results$comparison <- paste(cluster_label, "vs All Active")
    
    # Update state
    state$de_genes(de_results)
    state$de_status("completed")
    state$heatmap_data(NULL)
    state$heatmap_type(NULL)
    
    print("One vs all analysis completed")
  })
}

# Run Pairwise DE analysis
runPairwiseAnalysis <- function(seurat_obj, cluster1, cluster2, state) {
  withProgress(message = 'Computing pairwise differential expression...', {
    # Fix: seurat_obj is an object, not a function
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
    
    # Add comparison label
    cluster1_label <- getClusterLabel(cluster1, state$cluster_labels())
    cluster2_label <- getClusterLabel(cluster2, state$cluster_labels())
    
    de_results$comparison <- paste(cluster1_label, "vs", cluster2_label)
    
    # Update state
    state$de_genes(de_results)
    state$de_status("completed")
    state$heatmap_data(NULL)
    state$heatmap_type(NULL)
    
    print("Pairwise analysis completed")
  })
}

# Run General Heatmap analysis
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
        
        for (i in seq_along(active_clusters_num)) {
          cluster <- active_clusters_num[i]
          
          incProgress(1 / length(active_clusters_num), 
                      detail = paste("Finding markers for cluster", 
                                     getClusterLabel(cluster, state$cluster_labels())))
          
          # Find markers for this cluster vs all other active clusters
          # Fix: seurat_obj is an object, not a function
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
          } else {
            print(paste("No DE genes found for cluster", cluster))
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
        print("General Heatmap analysis completed")
      }, error = function(e) {
        showNotification(paste("Error computing general heatmap:", e$message), type = "error")
        print(paste("Error in general heatmap:", e$message))
      })
    })
  } else {
    showNotification("No active clusters selected. Please select at least one cluster.", 
                     type = "warning")
    state$general_heatmap_genes(NULL)
  }
}

# Generate DE heatmap
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

# Render DE Results UI
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

# Render General Heatmap UI
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

# Setup download handlers for heatmaps
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
      current_labels <- state$cluster_labels()
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
      current_labels <- state$cluster_labels()
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