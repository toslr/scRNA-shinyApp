# de_analysis_module.R

deAnalysisUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Static container for cluster management
    wellPanel(
      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
          h4(style = "margin: 0;", "Cluster Management"),
          actionButton(ns("updateAllLabels"), "Update All Cluster Labels",
                       class = "btn-primary")
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
    
    cluster_labels <- reactiveVal(NULL)
    active_clusters <- reactiveVal(NULL)
    temp_labels <- reactiveVal(NULL)
    
    
    label_inputs <- reactiveValues()
    
    clusters <- reactive({
      req(clustered_seurat())
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      sort(unique(clustered_seurat()$seurat_clusters))
    })
    
    # Initialize cluster labels and active status when clusters change
    observe({
      req(clusters())
      current_clusters <- clusters()
      
      # Initialize or update labels if needed
      current_labels <- cluster_labels()
      if (is.null(current_labels) || 
          !all(as.character(current_clusters) %in% names(current_labels))) {
        
        new_labels <- initializeClusterLabels(current_clusters, current_labels)
        
        cluster_labels(new_labels)
        temp_labels(new_labels)
      }
      
      # Initialize active status if needed
      current_active <- active_clusters()
      if (is.null(current_active) || 
          !all(as.character(current_clusters) %in% names(current_active))) {
        active_clusters(initializeActiveStatus(current_clusters, current_active))
      }
      
    })
    
    # Get currently active clusters
    active_cluster_list <- reactive({
      req(active_clusters())
      names(active_clusters()[active_clusters() == TRUE])
    })
    
    # Render cluster controls
    output$clusterControls <- renderUI({
      req(clustered_seurat())
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      req(clusters())
      req(temp_labels())
      req(active_clusters())
      
      available_clusters <- clusters()
      current_temp_labels <- temp_labels()
      current_active <- active_clusters()
      
      createClusterControls(ns, available_clusters, current_temp_labels, current_active)
    })
    
    # Update label_inputs when text changes
    observe({
      req(clusters())
      for(cluster in clusters()) {
        local({
          local_cluster <- cluster
          input_id <- paste0("label_", local_cluster)
          
          observeEvent(input[[input_id]], {
            label_inputs[[input_id]] <- input[[input_id]]
          }, ignoreInit = TRUE)
        })
      }
    })
    
    # Update temp_labels only when update button is clicked
    observeEvent(input$updateAllLabels, {
      req(clusters())
      current_temp <- temp_labels()
      
      for(cluster in clusters()) {
        input_id <- paste0("label_", cluster)
        if(!is.null(label_inputs[[input_id]])) {
          current_temp[as.character(cluster)] <- label_inputs[[input_id]]
        }
      }
      
      temp_labels(current_temp)
      cluster_labels(current_temp)
    })
    
    # Handle active status updates
    observe({
      req(clusters())
      available_clusters <- clusters()
      
      lapply(available_clusters, function(cluster) {
        observeEvent(input[[paste0("active_", cluster)]], {
          req(active_clusters())
          current_active <- active_clusters()
          current_active[as.character(cluster)] <- input[[paste0("active_", cluster)]]
          active_clusters(current_active)
        })
      })
    })
    
    # Render the analysis UI separately
    output$analysisUI <- renderUI({
      req(clustered_seurat())
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      req(clusters())
      req(cluster_labels())
      req(active_clusters())
      
      available_clusters <- clusters()
      current_labels <- cluster_labels()
      current_active <- active_clusters()
      
      # Create choices for select inputs (only active clusters)
      active_cluster_ids <- as.numeric(names(current_active[current_active == TRUE]))
      cluster_choices <- setNames(
        active_cluster_ids,
        vapply(active_cluster_ids, function(x) current_labels[[as.character(x)]], character(1))
      )
      
      createAnalysisUI(ns, cluster_choices)
    })
    
    de_genes <- reactiveVal()
    de_status <- reactiveVal(NULL)
    
    # One vs All analysis (modified to handle gene names properly)
    observeEvent(input$runDEAll, {
      req(clustered_seurat(), input$targetClusterAll)
      req(active_cluster_list())
      
      withProgress(message = 'Computing one vs all differential expression...', {
        # Subset Seurat object to only include active clusters
        active_cells <- clustered_seurat()$seurat_clusters %in% active_cluster_list()
        seurat_subset <- subset(clustered_seurat(), cells = colnames(clustered_seurat())[active_cells])
        
        de_results <- performDEanalysis(seurat_subset, input$targetClusterAll)
        
        de_results$comparison <- paste(getClusterLabel(input$targetClusterAll,cluster_labels()), "vs All Active")
        de_genes(de_results)
        de_status("completed")
      })
    })
    
    # Pairwise analysis (modified to handle gene names properly)
    observeEvent(input$runDEPair, {
      req(clustered_seurat(), input$targetCluster1, input$targetCluster2)
      req(input$targetCluster1 != input$targetCluster2)
      
      withProgress(message = 'Computing pairwise differential expression...', {
        
        de_results <- performDEanalysis(clustered_seurat(), input$targetCluster1, input$targetCluster2)
        
        de_results$comparison <- paste(getClusterLabel(input$targetCluster1, cluster_labels()),
                                       "vs",
                                       getClusterLabel(input$targetCluster2, cluster_labels()))
        de_genes(de_results)
        de_status("completed")
      })
    })
    
    # Single volcano plot
    output$volcanoPlot <- renderPlot({
      req(de_genes())
      results <- de_genes()
      createVolcanoPlot(results)
    })
    
    # Single results table
    output$deTable <- DT::renderDataTable({
      req(de_genes())
      DT::datatable(de_genes(),
                    options = list(pageLength = 10,
                                   scrollX = TRUE),
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
                   actionButton(ns("generate_heatmap"), "Generate Heatmap")
            )
          )
        )
      )
    })
    
    # Generate heatmap
    observeEvent(input$generate_heatmap, {
      req(clustered_seurat(), input$top_n_genes, de_genes())
      
      withProgress(message = 'Generating heatmap...', {
        de_results <- de_genes()
        # Get top genes
        top_genes <- rownames(de_results[order(de_results$p_val_adj), ])[1:input$top_n_genes]
        # Get expression data
        seurat_obj <- clustered_seurat()
        # Generate heatmap
        output$heatmapPlot <- renderPlot({
          createExpressionHeatmap(seurat_obj, top_genes, cluster_labels())
        })
      })
    })
    
    # Return results, status, and labels
    return(list(
      results = de_genes,
      status = de_status,
      labels = cluster_labels,
      active=active_clusters
    ))
  })
}