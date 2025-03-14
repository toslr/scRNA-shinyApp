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
    
    get_cluster_label <- function(cluster) {
      labels <- cluster_labels()
      if (!is.null(labels) && !is.null(labels[as.character(cluster)])) {
        labels[as.character(cluster)]
      } else {
        paste("Cluster", cluster)
      }
    }
    
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
        
        # Create new labels for missing clusters
        new_labels <- setNames(
          paste("Cluster", current_clusters), 
          as.character(current_clusters)
        )
        
        # Merge with existing labels if they exist
        if (!is.null(current_labels)) {
          existing_clusters <- names(current_labels)
          new_labels[existing_clusters] <- current_labels[existing_clusters]
        }
        
        cluster_labels(new_labels)
        temp_labels(new_labels)
      }
      
      # Initialize active status if needed
      current_active <- active_clusters()
      if (is.null(current_active) || 
          !all(as.character(current_clusters) %in% names(current_active))) {
        
        # Set all clusters as active by default
        new_active <- setNames(
          rep(TRUE, length(current_clusters)),
          as.character(current_clusters)
        )
        
        # Merge with existing active status if it exists
        if (!is.null(current_active)) {
          existing_clusters <- names(current_active)
          new_active[existing_clusters] <- current_active[existing_clusters]
        }
        
        active_clusters(new_active)
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
      
      tagList(
        lapply(available_clusters, function(cluster) {
          div(style = "margin-bottom: 10px;",
              fluidRow(
                column(1,
                       checkboxInput(ns(paste0("active_", cluster)), 
                                     label=NULL,
                                     value = current_active[[as.character(cluster)]])
                ),
                column(4, 
                       tags$label(paste("Cluster", cluster)),
                       textInput(ns(paste0("label_", cluster)),
                                 label = NULL,
                                 value = current_temp_labels[[as.character(cluster)]])
                )

              )
          )
        })
      )
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
      
      tagList(
        fluidRow(
          column(6,
                 wellPanel(
                   h4("One vs All Analysis"),
                   selectInput(ns("targetClusterAll"), 
                               "Select cluster to compare against all others:", 
                               choices = cluster_choices,
                               selected = NULL),
                   actionButton(ns("runDEAll"), "Run One vs All DE")
                 )
          ),
          column(6,
                 wellPanel(
                   h4("One vs One Analysis"),
                   selectInput(ns("targetCluster1"), 
                               "Select first cluster:", 
                               choices = cluster_choices),
                   selectInput(ns("targetCluster2"), 
                               "Select second cluster:", 
                               choices = cluster_choices),
                   actionButton(ns("runDEPair"), "Run Pairwise DE")
                 )
          )
        ),
        plotOutput(ns("volcanoPlot"), height = "400px"),
        DT::dataTableOutput(ns("deTable"))
      )
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
        
        de_results <- FindMarkers(seurat_subset,
                                  ident.1 = input$targetClusterAll,
                                  ident.2 = NULL,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)
        
        # Add gene names - with error handling
        gene_mapping <- clustered_seurat()@misc$gene_mapping
        if (!is.null(gene_mapping)) {
          # Map gene names and handle missing mappings
          de_results$gene <- gene_mapping[rownames(de_results)]
          de_results$gene[is.na(de_results$gene)] <- rownames(de_results)[is.na(de_results$gene)]
        } else {
          warning("Gene mapping not found in Seurat object")
          de_results$gene <- rownames(de_results)
        }
        
        de_results$comparison <- paste(get_cluster_label(input$targetClusterAll), "vs All Active")
        de_genes(de_results)
        de_status("completed")
      })
    })
    
    # Pairwise analysis (modified to handle gene names properly)
    observeEvent(input$runDEPair, {
      req(clustered_seurat(), input$targetCluster1, input$targetCluster2)
      req(input$targetCluster1 != input$targetCluster2)
      
      withProgress(message = 'Computing pairwise differential expression...', {
        de_results <- FindMarkers(clustered_seurat(),
                                  ident.1 = input$targetCluster1,
                                  ident.2 = input$targetCluster2,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)
        
        # Add gene names - with error handling
        gene_mapping <- clustered_seurat()@misc$gene_mapping
        if (!is.null(gene_mapping)) {
          # Map gene names and handle missing mappings
          de_results$gene <- gene_mapping[rownames(de_results)]
          de_results$gene[is.na(de_results$gene)] <- rownames(de_results)[is.na(de_results$gene)]
        } else {
          warning("Gene mapping not found in Seurat object")
          de_results$gene <- rownames(de_results)
        }
        
        de_results$comparison <- paste(get_cluster_label(input$targetCluster1),
                                       "vs",
                                       get_cluster_label(input$targetCluster2))
        de_genes(de_results)
        de_status("completed")
      })
    })
    
    # Single volcano plot
    output$volcanoPlot <- renderPlot({
      req(de_genes())
      results <- de_genes()
      
      ggplot(results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(aes(color = abs(avg_log2FC) > 0.25 & p_val_adj < 0.05)) +
        scale_color_manual(values = c("grey", "red")) +
        theme_classic() +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        labs(title = unique(results$comparison),
             color = "Significant") +
        theme(legend.position = "bottom")
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
    
    # Return results, status, and labels
    return(list(
      results = de_genes,
      status = de_status,
      labels = cluster_labels,
      active=active_clusters
    ))
  })
}