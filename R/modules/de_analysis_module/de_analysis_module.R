# de_analysis_module.R

deAnalysisUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Static container for cluster management
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
    # Container for analysis UI
    uiOutput(ns("analysisUI"))
  )
}

deAnalysisServer <- function(id, clustered_seurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize reactive values
    cluster_labels <- reactiveVal(NULL)
    active_clusters <- reactiveVal(NULL)
    temp_labels <- reactiveVal(NULL)
    
    # Add reactive values to track heatmap state
    heatmap_data <- reactiveVal(NULL)
    heatmap_type <- reactiveVal(NULL)
    
    label_inputs <- reactiveValues()
    
    # DE analysis reactive values
    de_genes <- reactiveVal(NULL)
    de_status <- reactiveVal(NULL)
    analysis_state <- reactiveVal("none")
    general_heatmap_genes <- reactiveVal(NULL)
    general_heatmap_clusters <- reactiveVal(NULL)
    
    # Get available clusters from Seurat object
    clusters <- reactive({
      req(clustered_seurat())
      if (!"seurat_clusters" %in% colnames(clustered_seurat()@meta.data)) {
        return(NULL)
      }
      sort(unique(clustered_seurat()$seurat_clusters))
    })
    
    observeEvent(clustered_seurat(), {
      # Reset all heatmap-related data
      heatmap_data(NULL)
      heatmap_type(NULL)
      general_heatmap_genes(NULL)
      general_heatmap_clusters(NULL)
      analysis_state("none")
    }, ignoreInit = TRUE)
    
    # Initialize cluster labels and active status when clusters change
    observe({
      available_clusters <- clusters()
      
      # Skip if no clusters available
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      # Initialize or update labels if needed
      current_labels <- cluster_labels()
      if (is.null(current_labels) || 
          !all(as.character(available_clusters) %in% names(current_labels))) {
        
        new_labels <- initializeClusterLabels(available_clusters, current_labels)
        
        cluster_labels(new_labels)
        temp_labels(new_labels)
      }
      
      # Initialize active status if needed
      current_active <- active_clusters()
      if (is.null(current_active) || 
          !all(as.character(available_clusters) %in% names(current_active))) {
        active_clusters(initializeActiveStatus(available_clusters, current_active))
      }
    })
    
    # Get currently active clusters
    active_cluster_list <- reactive({
      current_active <- active_clusters()
      if (is.null(current_active) || length(current_active) == 0) {
        return(character(0))
      }
      names(current_active[current_active == TRUE])
    })
    
    # Render cluster controls
    output$clusterControls <- renderUI({
      available_clusters <- clusters()
      
      # Show message if no clusters
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(div(
          class = "alert alert-info",
          "No clusters available. Please run the clustering step first."
        ))
      }
      
      # Get current labels and active status
      current_temp_labels <- temp_labels()
      current_active <- active_clusters()
      
      # Return UI
      createClusterControls(ns, available_clusters, current_temp_labels, current_active)
    })
    
    # Handle Select All button
    observeEvent(input$selectAllClusters, {
      available_clusters <- clusters()
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      current_active <- active_clusters()
      for (cluster in available_clusters) {
        cluster_key <- as.character(cluster)
        current_active[cluster_key] <- TRUE
      }
      active_clusters(current_active)
      
      # Clear general heatmap if it exists
      if (!is.null(heatmap_type()) && heatmap_type() == "general") {
        general_heatmap_genes(NULL)
      }
    })
    
    # Handle Deselect All button
    observeEvent(input$deselectAllClusters, {
      available_clusters <- clusters()
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      current_active <- active_clusters()
      for (cluster in available_clusters) {
        cluster_key <- as.character(cluster)
        current_active[cluster_key] <- FALSE
      }
      active_clusters(current_active)
      
      # Clear general heatmap if it exists
      if (!is.null(heatmap_type()) && heatmap_type() == "general") {
        general_heatmap_genes(NULL)
      }
    })
    
    # Update label_inputs when text changes
    observe({
      available_clusters <- clusters()
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
              label_inputs[[input_id]] <- input[[input_id]]
            }, ignoreInit = TRUE)
          }
        })
      }
    })
    
    # Update temp_labels only when update button is clicked
    observeEvent(input$updateAllLabels, {
      available_clusters <- clusters()
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      current_temp <- temp_labels()
      
      for(cluster in available_clusters) {
        input_id <- paste0("label_", cluster)
        if(!is.null(label_inputs[[input_id]])) {
          current_temp[as.character(cluster)] <- label_inputs[[input_id]]
        }
      }
      
      temp_labels(current_temp)
      cluster_labels(current_temp)
    })
    
    # Handle active status updates - reset heatmaps when active clusters change
    observe({
      available_clusters <- clusters()
      
      # Skip empty clusters
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(NULL)
      }
      
      lapply(available_clusters, function(cluster) {
        input_id <- paste0("active_", cluster)
        
        # Only observe if this input exists
        if (!is.null(input[[input_id]])) {
          observeEvent(input[[input_id]], {
            current_active <- active_clusters()
            if (is.null(current_active)) {
              return(NULL)
            }
            
            cluster_key <- as.character(cluster)
            
            # Check if this cluster exists in our active_clusters
            if (cluster_key %in% names(current_active)) {
              current_active[cluster_key] <- input[[input_id]]
              active_clusters(current_active)
              
              # Reset heatmap data when active clusters change
              if (!is.null(heatmap_type()) && heatmap_type() == "general") {
                # Only clear general heatmap as it's affected by active clusters
                general_heatmap_genes(NULL)
              }
            }
          }, ignoreInit = TRUE)
        }
      })
    })
    
    # Render the analysis UI separately
    output$analysisUI <- renderUI({
      available_clusters <- clusters()
      if (is.null(available_clusters) || length(available_clusters) == 0) {
        return(div(
          class = "alert alert-info",
          "No clusters available. Please run the clustering step first."
        ))
      }
      
      # Get current labels and active status
      current_labels <- cluster_labels()
      current_active <- active_clusters()
      
      if (is.null(current_labels) || is.null(current_active)) {
        return(div(
          class = "alert alert-warning",
          "Cluster information is being loaded. Please wait..."
        ))
      }
      
      # Create choices for select inputs (only active clusters)
      active_cluster_ids <- as.numeric(names(current_active[current_active == TRUE]))
      
      # Check if we have any active clusters
      if (length(active_cluster_ids) == 0) {
        return(div(
          class = "alert alert-warning",
          "No active clusters selected. Please activate at least one cluster in the Cluster Management section above."
        ))
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
      
      createAnalysisUI(ns, cluster_choices)
    })
    
    # Clear state when starting a new analysis
    clear_state <- function(new_state) {
      analysis_state(new_state)
      
      # Clear visualizations that aren't needed for the new state
      if (new_state == "one_vs_all" || new_state == "pairwise") {
        # Reset general heatmap when switching to DE analysis
        general_heatmap_genes(NULL)
        general_heatmap_clusters(NULL)
      } else if (new_state == "general_heatmap") {
        # Clear specific heatmap data when switching to general
        heatmap_data(NULL)
      }
    }
    
    # One vs All analysis
    observeEvent(input$runDEAll, {
      req(clustered_seurat(), input$targetClusterAll)
      current_active_list <- active_cluster_list()
      
      if (length(current_active_list) == 0) {
        showNotification("No active clusters selected.", type = "warning")
        return(NULL)
      }
      
      print("Running one vs all analysis")
      clear_state("one_vs_all")
      
      withProgress(message = 'Computing one vs all differential expression...', {
        # Subset Seurat object to only include active clusters
        active_cells <- clustered_seurat()$seurat_clusters %in% current_active_list
        seurat_subset <- subset(clustered_seurat(), cells = colnames(clustered_seurat())[active_cells])
        
        # Check if we have enough active clusters
        if (length(unique(seurat_subset$seurat_clusters)) < 2) {
          showNotification("Need at least two active clusters for comparison.", type = "warning")
          return(NULL)
        }
        
        de_results <- performDEanalysis(seurat_subset, input$targetClusterAll)
        
        # Check if we got valid results
        if (is.null(de_results) || nrow(de_results) == 0) {
          showNotification("No differential expression results found.", type = "warning")
          return(NULL)
        }
        
        # Add comparison label
        cluster_label <- if (!is.null(cluster_labels())) {
          getClusterLabel(input$targetClusterAll, cluster_labels())
        } else {
          paste("Cluster", input$targetClusterAll)
        }
        
        de_results$comparison <- paste(cluster_label, "vs All Active")
        de_genes(de_results)
        de_status("completed")
        
        # Clear any existing heatmap when running new DE
        heatmap_data(NULL)
        heatmap_type(NULL)
        
        print("One vs all analysis completed")
      })
    })
    
    # Pairwise analysis
    observeEvent(input$runDEPair, {
      req(clustered_seurat(), input$targetCluster1, input$targetCluster2)
      
      if (input$targetCluster1 == input$targetCluster2) {
        showNotification("Please select different clusters for comparison.", type = "warning")
        return(NULL)
      }
      
      print("Running pairwise analysis")
      clear_state("pairwise")
      
      withProgress(message = 'Computing pairwise differential expression...', {
        de_results <- performDEanalysis(
          clustered_seurat(), 
          input$targetCluster1, 
          input$targetCluster2
        )
        
        # Check if we got valid results
        if (is.null(de_results) || nrow(de_results) == 0) {
          showNotification("No differential expression results found.", type = "warning")
          return(NULL)
        }
        
        # Add comparison label
        cluster1_label <- if (!is.null(cluster_labels())) {
          getClusterLabel(input$targetCluster1, cluster_labels())
        } else {
          paste("Cluster", input$targetCluster1)
        }
        
        cluster2_label <- if (!is.null(cluster_labels())) {
          getClusterLabel(input$targetCluster2, cluster_labels())
        } else {
          paste("Cluster", input$targetCluster2)
        }
        
        de_results$comparison <- paste(cluster1_label, "vs", cluster2_label)
        de_genes(de_results)
        de_status("completed")
        
        # Clear any existing heatmap when running new DE
        heatmap_data(NULL)
        heatmap_type(NULL)
        
        print("Pairwise analysis completed")
      })
    })
    
    # General heatmap analysis
    observeEvent(input$runGeneralHeatmap, {
      req(clustered_seurat(), input$genesPerCluster)
      current_active_list <- active_cluster_list()
      
      if (length(current_active_list) == 0) {
        showNotification("No active clusters selected.", type = "warning")
        # Clear any previous heatmap data to prevent errors
        general_heatmap_genes(NULL)
        general_heatmap_clusters(NULL)
        return(NULL)
      }
      
      print("Running general heatmap analysis")
      clear_state("general_heatmap")
      
      # Always compute from scratch to reflect current active clusters
      active_clusters_num <- as.numeric(current_active_list)
      
      # Only proceed if we have active clusters
      if (length(active_clusters_num) > 0) {
        withProgress(message = "Computing general heatmap...", value = 0, detail = "Preparing analysis", {
          tryCatch({
            # Collect DE genes for each active cluster independently
            all_de_results <- list()
            
            # Get cluster info
            total_clusters <- length(active_clusters_num)
            
            for (i in seq_along(active_clusters_num)) {
              cluster <- active_clusters_num[i]
              
              incProgress(i / active_clusters_num, 
                          detail = paste("Finding markers for cluster", 
                                         getClusterLabel(cluster, cluster_labels())))
              
              # For each active cluster, find markers vs all other active clusters
              res <- FindMarkers(
                clustered_seurat(),
                ident.1 = cluster,
                ident.2 = setdiff(active_clusters_num, cluster),
                min.pct = 0.25,
                logfc.threshold = 0.25
              )
              
              # Add gene names
              res <- addGeneNames(res, clustered_seurat())
              
              # Get top N genes for this cluster
              if (nrow(res) > 0) {
                # Sort by p-value and fold change
                res <- res[order(res$p_val_adj, -abs(res$avg_log2FC)), ]
                top_cluster_genes <- rownames(res)[1:min(input$genesPerCluster, nrow(res))]
                all_de_results[[as.character(cluster)]] <- top_cluster_genes
              } else {
                print(paste("No DE genes found for cluster", cluster))
              }
            }
            
            # Store cluster order for later use in visualization
            general_heatmap_clusters(names(all_de_results))
            
            incProgress(0.8, detail = "Combining results")
            
            # Combine all unique genes from every cluster
            unique_genes <- unique(unlist(all_de_results))
            
            if (is.null(unique_genes) || length(unique_genes) == 0) {
              showNotification("No significant genes found for selected clusters.", type = "warning")
              return(NULL)
            }
            
            # Create order for genes, first by cluster then by significance
            gene_order <- c()
            for (cluster in names(all_de_results)) {
              gene_order <- c(gene_order, all_de_results[[cluster]])
            }
            # Remove duplicates but keep order
            gene_order <- unique(gene_order)
            
            incProgress(0.9, detail = "Building heatmap")
            
            print(paste("Computed", length(unique_genes), "total unique genes"))
            print(paste("Gene order has", length(gene_order), "genes"))
            
            # Store genes in order by cluster for diagonal-like heatmap
            general_heatmap_genes(gene_order)
            heatmap_type("general")
            
            incProgress(1.0, detail = "Completed")
            print("General Heatmap analysis completed")
          }, error = function(e) {
            showNotification(paste("Error computing general heatmap:", e$message), type = "error")
            print(paste("Error in general heatmap:", e$message))
          })
        })
      } else {
        # Handle the case with no active clusters
        showNotification("No active clusters selected. Please select at least one cluster.", 
                         type = "warning")
        general_heatmap_genes(NULL)
      }
    })
    
    # Render DE results UI
    output$deResultsUI <- renderUI({
      current_state <- analysis_state()
      
      if (is.null(current_state) || current_state == "none") {
        return(NULL)
      }
      
      if (current_state %in% c("one_vs_all", "pairwise")) {
        current_de_genes <- de_genes()
        
        if (is.null(current_de_genes) || nrow(current_de_genes) == 0) {
          return(div(
            class = "alert alert-info",
            "No differential expression results available yet. Run an analysis to see results."
          ))
        }
        
        has_heatmap <- !is.null(heatmap_data()) && 
          !is.null(heatmap_type()) && 
          heatmap_type() == "specific"
        
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
          uiOutput(ns("heatmapControls")),
          
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
    })
    
    output$generalHeatmapUI <- renderUI({
      current_state <- analysis_state()
      current_genes <- general_heatmap_genes()
      
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
    })
    
    # General heatmap - dependent on general_heatmap_genes and active clusters
    output$generalHeatmapPlot <- renderPlot({
      req(general_heatmap_genes())
      req(clustered_seurat())
      current_labels <- cluster_labels()
      current_active <- active_cluster_list()
      
      if (length(current_active) == 0) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "No active clusters selected") + 
                 theme_void())
      }
      
      print("Rendering general heatmap...")
      # Only include active clusters in the heatmap
      active_clusters_num <- as.numeric(current_active)
      
      # Get the ordered gene list
      ordered_genes <- general_heatmap_genes()
      
      # Create the heatmap with ordered genes
      createGeneralHeatmap(clustered_seurat(), 
                           ordered_genes, 
                           current_labels, 
                           active_clusters_num,
                           cluster_order = general_heatmap_clusters())
    })
    
    # Single volcano plot
    volcano_plot <- reactive({
      req(de_genes())
      results <- de_genes()
      if (nrow(results) == 0) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "No significant DE genes found") + 
                 theme_void())
      }
      
      createVolcanoPlot(results)
    })
    
    output$volcanoPlot <- renderPlot({
      volcano_plot()
    })
    
    output$downloadVolcanoPlot <- downloadHandler(
      filename = function() {
        paste("volcano_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        # Create high-resolution PNG device directly
        png(file, width = 3000, height = 2400, res = 300)
        
        # Generate the plot directly to the device
        results <- de_genes()
        print(createVolcanoPlot(results))
        
        # Close the device to save the file
        dev.off()
      }
    )
    
    # Single results table - with only CSV export
    output$deTable <- DT::renderDataTable({
      req(de_genes())
      results <- de_genes()
      
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
    
    # Show heatmap controls after DE analysis
    output$heatmapControls <- renderUI({
      req(de_genes())
      
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
      )
    })
    
    # Generate specific DE heatmap
    observeEvent(input$generate_heatmap, {
      req(clustered_seurat(), input$top_n_genes, de_genes())
      
      withProgress(message = 'Generating heatmap...', {
        de_results <- de_genes()
        
        # Ensure we don't request more genes than we have
        n_genes <- min(input$top_n_genes, nrow(de_results))
        
        if (n_genes == 0) {
          showNotification("No genes available for heatmap.", type = "warning")
          return(NULL)
        }
        
        # Get top genes - sort by p-value and then fold change
        de_results <- de_results[order(de_results$p_val_adj, -abs(de_results$avg_log2FC)), ]
        top_genes <- rownames(de_results)[1:n_genes]
        heatmap_data(top_genes)
        heatmap_type("specific")
        
        # Get active clusters for this heatmap
        current_active_list <- active_cluster_list()
        
        # Always use the active clusters for the heatmap
        if (length(current_active_list) > 0) {
          # For both 1vsAll and pairwise, show all active clusters
          print("Using all active clusters for heatmap")
        } else {
          print("No active clusters for heatmap")
        }
      })
    })
    
    # Render the specific heatmap - only when heatmap_data is available
    heatmap_plot <- reactive({
      current_heatmap_data <- heatmap_data()
      current_heatmap_type <- heatmap_type()
      
      if (is.null(current_heatmap_data) || 
          is.null(current_heatmap_type) || 
          current_heatmap_type != "specific") {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "Click 'Generate Heatmap' to visualize gene expression") + 
                 theme_void())
      }
      
      req(clustered_seurat())
      current_labels <- cluster_labels()
      
      # Use active clusters for the heatmap
      active_clusters_num <- as.numeric(active_cluster_list())
      
      # For all analysis types, show expression in all active clusters
      createExpressionHeatmap(clustered_seurat(), 
                              current_heatmap_data, 
                              current_labels,
                              active_clusters_num)
    })
    
    output$heatmapPlot <- renderPlot({
      heatmap_plot()
    })
    
    output$downloadHeatmapPlot <- downloadHandler(
      filename = function() {
        paste("heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        # Create high-resolution PNG device directly for pheatmap
        png(file, width = 3000, height = 3600, res = 300)
        
        # Generate the heatmap directly
        current_heatmap_data <- heatmap_data()
        current_labels <- cluster_labels()
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
    
    output$downloadGeneralHeatmapPlot <- downloadHandler(
      filename = function() {
        paste("general_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        # Create high-resolution PNG device directly for pheatmap
        png(file, width = 3600, height = 4200, res = 300)
        
        # Generate the heatmap directly
        current_labels <- cluster_labels()
        current_active <- active_cluster_list()
        active_clusters_num <- as.numeric(current_active)
        ordered_genes <- general_heatmap_genes()
        
        # Create the heatmap with ordered genes
        createGeneralHeatmap(clustered_seurat(), 
                             ordered_genes, 
                             current_labels, 
                             active_clusters_num,
                             cluster_order = general_heatmap_clusters())
        
        # Close the device to save the file
        dev.off()
      }
    )
    
    # Return results, status, and labels
    return(list(
      results = de_genes,
      status = de_status,
      labels = cluster_labels,
      active = active_clusters
    ))
  })
}