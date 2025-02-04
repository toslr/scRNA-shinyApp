# de_analysis_module.R

deAnalysisUI <- function(id) {
  ns <- NS(id)
  # Initially show nothing, content will be rendered conditionally
  uiOutput(ns("deUI"))
}

deAnalysisServer <- function(id, clustered_seurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive value to store cluster labels
    cluster_labels <- reactiveVal(NULL)
    
    # Reactive to safely get cluster information
    clusters <- reactive({
      req(clustered_seurat())
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      sort(unique(clustered_seurat()$seurat_clusters))
    })
    
    # Initialize cluster labels when clusters change
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
      }
    })
    
    # Render the DE analysis UI
    output$deUI <- renderUI({
      req(clustered_seurat())
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      req(clusters())
      req(cluster_labels())
      
      available_clusters <- clusters()
      current_labels <- cluster_labels()
      
      # Create choices for select inputs
      cluster_choices <- setNames(
        available_clusters,
        vapply(available_clusters, function(x) current_labels[[as.character(x)]], character(1))
      )
      
      tagList(
        fluidRow(
          column(12,
                 wellPanel(
                   h4("Cluster Management"),
                   div(style = "max-height: 300px; overflow-y: auto;",
                       lapply(available_clusters, function(cluster) {
                         div(style = "margin-bottom: 10px;",
                             fluidRow(
                               column(4, 
                                      tags$label(paste("Cluster", cluster)),
                                      textInput(ns(paste0("label_", cluster)),
                                                label = NULL,
                                                value = current_labels[[as.character(cluster)]])
                               ),
                               column(2,
                                      div(style = "margin-top: 5px;",
                                          actionButton(ns(paste0("update_", cluster)), "Update",
                                                       class = "btn-sm")
                                      )
                               )
                             )
                         )
                       })
                   )
                 )
          )
        ),
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
    
    # Handle cluster label updates
    observe({
      req(clusters())
      available_clusters <- clusters()
      
      lapply(available_clusters, function(cluster) {
        observeEvent(input[[paste0("update_", cluster)]], {
          req(cluster_labels())
          new_label <- input[[paste0("label_", cluster)]]
          current_labels <- cluster_labels()
          current_labels[as.character(cluster)] <- new_label
          cluster_labels(current_labels)
        })
      })
    })
    
    # Reactive value to store current DE results
    de_genes <- reactiveVal()
    # Reactive value to track whether DE was run
    de_status <- reactiveVal(NULL)
    
    # Get cluster label for plot titles
    get_cluster_label <- function(cluster) {
      labels <- cluster_labels()
      if (!is.null(labels) && !is.null(labels[as.character(cluster)])) {
        labels[as.character(cluster)]
      } else {
        paste("Cluster", cluster)
      }
    }
    
    # One vs All analysis
    observeEvent(input$runDEAll, {
      req(clustered_seurat(), input$targetClusterAll)
      withProgress(message = 'Computing one vs all differential expression...', {
        de_results <- FindMarkers(clustered_seurat(),
                                  ident.1 = input$targetClusterAll,
                                  ident.2 = NULL,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)
        
        de_results$gene <- clustered_seurat()@misc$gene_mapping[rownames(de_results)]
        de_results$comparison <- paste(get_cluster_label(input$targetClusterAll), "vs All")
        de_genes(de_results)
        de_status("completed")
      })
    })
    
    # Pairwise analysis
    observeEvent(input$runDEPair, {
      req(clustered_seurat(), input$targetCluster1, input$targetCluster2)
      req(input$targetCluster1 != input$targetCluster2)
      
      withProgress(message = 'Computing pairwise differential expression...', {
        de_results <- FindMarkers(clustered_seurat(),
                                  ident.1 = input$targetCluster1,
                                  ident.2 = input$targetCluster2,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)
        
        de_results$gene <- clustered_seurat()@misc$gene_mapping[rownames(de_results)]
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
      labels = cluster_labels
    ))
  })
}