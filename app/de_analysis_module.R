# de_analysis_module.R

deAnalysisUI <- function(id) {
  ns <- NS(id)
  # Initially show nothing, content will be rendered conditionally
  uiOutput(ns("deUI"))
}

deAnalysisServer <- function(id, clustered_seurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Render the DE analysis UI only when clustering is done
    output$deUI <- renderUI({
      req(clustered_seurat())
      req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
      
      tagList(
        fluidRow(
          column(6,
                 wellPanel(
                   h4("One vs All Analysis"),
                   selectInput(ns("targetClusterAll"), 
                               "Select cluster to compare against all others:", 
                               choices = sort(unique(clustered_seurat()$seurat_clusters))),
                   actionButton(ns("runDEAll"), "Run One vs All DE")
                 )
          ),
          column(6,
                 wellPanel(
                   h4("One vs One Analysis"),
                   selectInput(ns("targetCluster1"), 
                               "Select first cluster:", 
                               choices = sort(unique(clustered_seurat()$seurat_clusters))),
                   selectInput(ns("targetCluster2"), 
                               "Select second cluster:", 
                               choices = sort(unique(clustered_seurat()$seurat_clusters))),
                   actionButton(ns("runDEPair"), "Run Pairwise DE")
                 )
          )
        ),
        plotOutput(ns("volcanoPlot"), height = "400px"),
        DT::dataTableOutput(ns("deTable"))
      )
    })
    
    # Reactive value to store current DE results
    de_genes <- reactiveVal()
    # Reactive value to track whether DE was run
    de_status <- reactiveVal(NULL)
    
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
        de_results$comparison <- paste("Cluster", input$targetClusterAll, "vs All")
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
        de_results$comparison <- paste("Cluster", input$targetCluster1, 
                                       "vs Cluster", input$targetCluster2)
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
    
    # Return both the results and status
    return(list(
      results = de_genes,
      status = de_status
    ))
  })
}