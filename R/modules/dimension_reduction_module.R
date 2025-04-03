# R/modules/dimension_reduction_module.R

#' @title Dimension Reduction Module UI
#' @description Creates the UI for the dimension reduction module which allows users to 
#'   visualize scRNA-seq data using PCA, UMAP, and clustering techniques.
#' @param id The module ID
#' @return A Shiny UI element containing the dimension reduction interface
#' @export
dimensionReductionUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Elbow plot section
    div(
      uiOutput(ns("elbowSaveButton")),
      plotOutput(ns("elbowPlot"), height = "400px"),
      textOutput(ns("suggestedDims")),
    ),
    
    # Dimension confirmation section
    div(
      style = "margin-top: 15px; margin-bottom: 20px;",
      fluidRow(
        column(6, 
               numericInput(ns("nDims"), 
                            "Number of dimensions for analysis:", 
                            value = 15,
                            min = 2,
                            max = 50)
        ),
        column(6,
               div(style = "margin-top: 25px;",
                   actionButton(ns("confirmDims"), "Confirm Dimensions", 
                                class = "btn-primary")
               )
        )
      )
    ),
    
    # Clustering controls (shown after dimensions are confirmed)
    uiOutput(ns("clusteringSection")),
    
    # UMAP visualization section (shown after clustering)
    uiOutput(ns("umapSection"))
  )
}

#' @title Dimension Reduction Module Server
#' @description Server logic for the dimension reduction module that performs PCA, UMAP,
#'   and clustering analyses, and creates interactive visualizations.
#' @param id The module ID
#' @param processed_seurat Reactive expression returning Seurat object with PCA
#' @param sample_management Optional sample management module for filtering by samples
#' @param condition_management Optional condition management module for filtering by conditions
#' @return Reactive expression returning the processed Seurat object with UMAP and clustering
#' @export
dimensionReductionServer <- function(id, processed_seurat, sample_management = NULL, condition_management = NULL, cluster_management = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # State management
    values <- reactiveValues(
      dims_confirmed = FALSE,
      clustering_done = FALSE,
      seurat_with_umap = NULL,
      clustered_seurat = NULL,
      left_umap_type = "2D",
      right_umap_type = "2D",
      left_color_by = "sample",
      right_color_by = "cluster",
      filtered_clusters = NULL,
      filtered_samples = NULL,
      filtered_condition = NULL,
      condition_column = NULL,
      left_plot_update_trigger = runif(1),
      right_plot_update_trigger = runif(1)
    )
    
    # Method to restore state from loaded data
    restoreState <- function(n_dims, dims_confirmed) {
      # Update UI to reflect loaded state
      if (!is.null(n_dims)) {
        updateNumericInput(session, "nDims", value = n_dims)
      }
      
      # If dimensions were confirmed in the saved state, update reactiveValues
      if (!is.null(dims_confirmed) && dims_confirmed == TRUE) {
        values$dims_confirmed <- TRUE
        
        # Make sure we have the processed Seurat object
        if (!is.null(processed_seurat())) {
          # Trigger UMAP computation silently
          withProgress(message = 'Computing UMAPs...', {
            # Apply sample and condition filtering if needed
            filtered_seurat <- processed_seurat()
            
            # Process UMAPs
            values$seurat_with_umap <- process_umaps(filtered_seurat, n_dims)
          })
        }
      }
    }
    
    # Watch for changes in sample management active status
    observe({
      req(sample_management)
      
      # Get active samples
      active_samples <- sample_management$getActiveSampleIds()
      
      # Update filtered samples
      if (!identical(values$filtered_samples, active_samples)) {
        values$filtered_samples <- active_samples
        
        # Force re-rendering of plots
        values$left_plot_update_trigger <- runif(1)
        values$right_plot_update_trigger <- runif(1)
      }
    })
    
    # Watch for changes in condition management active status
    observe({
      req(condition_management)
      
      # Get active conditions and the condition column
      active_conditions <- condition_management$getActiveConditions()
      condition_column <- condition_management$getConditionColumn()
      
      # Update filtered conditions and condition column
      if (!identical(values$filtered_conditions, active_conditions) || 
          !identical(values$condition_column, condition_column)) {
        
        values$filtered_conditions <- active_conditions
        values$condition_column <- condition_column
        
        # Force re-rendering of plots
        values$left_plot_update_trigger <- runif(1)
        values$right_plot_update_trigger <- runif(1)
      }
    })
    
    # Watch for changes in cluster management active status
    observe({
      req(cluster_management)
      
      # Get active clusters
      active_clusters <- cluster_management$getActiveClusterIds()
      
      # Update filtered clusters in the module's state
      values$filtered_clusters <- active_clusters
      
      # Force re-rendering of all plots by updating triggers
      values$left_plot_update_trigger <- runif(1)
      values$right_plot_update_trigger <- runif(1)
    })
    
    # Find elbow point for optimal dimensions
    suggested_dims <- reactive({
      req(processed_seurat())
      compute_suggested_dims(processed_seurat())
    })
    
    # Update suggested dimensions
    observe({
      req(suggested_dims())
      updateNumericInput(session, "nDims", value = suggested_dims())
    })
    
    # Create elbow plot
    elbow_plot <- reactive({
      req(processed_seurat(), suggested_dims())
      ElbowPlot(processed_seurat(), 
                ndims = ncol(Embeddings(processed_seurat(), "pca"))) +
        geom_vline(xintercept = suggested_dims(), 
                   color = "red", 
                   linetype = "dashed")
    })
    
    # Render elbow plot
    output$elbowPlot <- renderPlot({
      elbow_plot()
    })
    
    output$elbowSaveButton <- renderUI({
      req(processed_seurat(), suggested_dims())
      div(style = "margin-top: 10px; text-align: right;",
          downloadButton(ns("downloadElbowPlot"), "Save Plot", 
                         class = "btn-sm btn-success"))
    })
    
    # Download handler for elbow plot
    output$downloadElbowPlot <- downloadHandler(
      filename = function() {
        paste("pca_elbow_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
      },
      content = function(file) {
        save_elbow_plot(file, processed_seurat(), suggested_dims())
      }
    )
    
    output$suggestedDims <- renderText({
      req(suggested_dims())
      paste("Suggested number of dimensions (based on elbow point):", suggested_dims(),
            ". Please adjust the number of PC for reduction if needed.")
    })
    
    # Handle dimension confirmation
    observeEvent(input$confirmDims, {
      req(processed_seurat(), input$nDims)
      
      # Update state
      values$dims_confirmed <- TRUE
      
      # Pre-compute UMAPs for both 2D and 3D
      withProgress(message = 'Computing UMAPs...', {
        # Start with the original seurat object
        filtered_seurat <- processed_seurat()
        
        # Apply sample filtering if available
        if (!is.null(sample_management) && is.function(sample_management$getActiveSampleIds)) {
          active_samples <- sample_management$getActiveSampleIds()
          if (!is.null(active_samples) && length(active_samples) > 0) {
            filtered_seurat <- filter_by_samples(filtered_seurat, active_samples)
          }
        }
        
        # Apply condition filtering if available
        if (!is.null(condition_management) && 
            is.function(condition_management$getConditionColumn) && 
            is.function(condition_management$getActiveConditions)) {
          
          condition_column <- condition_management$getConditionColumn()
          active_conditions <- condition_management$getActiveConditions()
          
          if (!is.null(condition_column) && !is.null(active_conditions) && length(active_conditions) > 0) {
            filtered_seurat <- filter_by_conditions(filtered_seurat, condition_column, active_conditions)
          }
        }
        
        # Process UMAPs
        values$seurat_with_umap <- process_umaps(filtered_seurat, input$nDims)
      })
    })
    
    # Clustering section UI
    output$clusteringSection <- renderUI({
      req(values$dims_confirmed)
      
      div(
        style = "margin-top: 20px; margin-bottom: 20px;",
        h3("Clustering"),
        fluidRow(
          column(6, 
                 numericInput(ns("resolution"), 
                              "Clustering resolution:", 
                              0.5, 
                              min = 0, 
                              max = 2, 
                              step = 0.1)
          ),
          column(6,
                 div(style = "margin-top: 25px;",
                     actionButton(ns("runClustering"), "Run Clustering", 
                                  class = "btn-primary")
                 )
          )
        ),
        # Show only after clustering is done
        if (values$clustering_done) {
          div(
            style = "margin-top: 15px;",
            downloadButton(ns("downloadClusterStats"), "Download Cluster Statistics", 
                           class = "btn-sm btn-success")
          )
        }
      )
    })
    
    # Run clustering
    observeEvent(input$runClustering, {
      req(values$seurat_with_umap, input$nDims, input$resolution)
      
      withProgress(message = 'Clustering...', {
        # Run clustering
        values$clustered_seurat <- run_clustering(
          values$seurat_with_umap, 
          input$nDims, 
          input$resolution
        )
        values$clustering_done <- TRUE
      })
    })
    
    #' @title Get Unique Values
    #' @description Helper to get unique values for a metadata column
    #' @param column_name Name of the metadata column
    #' @return Vector of unique values for the specified column
    #' @keywords internal
    get_unique_values <- function(column_name) {
      req(values$clustered_seurat)
      
      if (column_name == "cluster" && "seurat_clusters" %in% colnames(values$clustered_seurat@meta.data)) {
        # For clusters, get them as characters but sort numerically
        cluster_values <- unique(as.character(values$clustered_seurat$seurat_clusters))
        # Convert to numeric, sort, then back to character
        return(as.character(sort(as.numeric(cluster_values))))
      } else if (column_name %in% colnames(values$clustered_seurat@meta.data)) {
        # For other metadata columns
        metadata_values <- unique(as.character(values$clustered_seurat@meta.data[[column_name]]))
        
        # Check if values are numeric and sort appropriately
        if (all(grepl("^\\d+$", metadata_values)) || 
            all(grepl("^\\d+\\.\\d+$", metadata_values))) {
          # If all values appear to be numeric, sort numerically
          return(as.character(sort(as.numeric(metadata_values))))
        } else {
          # Otherwise sort alphabetically
          return(sort(metadata_values))
        }
      } else {
        return(character(0))
      }
    }
    
    # Reactive values to track active items for each plot
    left_active_items <- reactiveVal(NULL)
    right_active_items <- reactiveVal(NULL)
    
    # Update active items when color_by changes
    observe({
      if (values$left_color_by == "cluster") {
        left_active_items(get_unique_values("cluster"))
      } else if (values$left_color_by == "sample") {
        left_active_items(get_unique_values("sample"))
      } else if (values$left_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        left_active_items(get_unique_values(values$left_color_by))
      }
    })
    
    observe({
      if (values$right_color_by == "cluster") {
        right_active_items(get_unique_values("cluster"))
      } else if (values$right_color_by == "sample") {
        right_active_items(get_unique_values("sample"))
      } else if (values$right_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        right_active_items(get_unique_values(values$right_color_by))
      }
    })
    
    # Toggle handlers for left plot
    observeEvent(input$left_items_select_all, {
      current_items <- if (values$left_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$left_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$left_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$left_color_by)
      } else {
        character(0)
      }
      
      for (item in current_items) {
        input_name <- paste0("left_items_", gsub("[^a-zA-Z0-9]", "_", item))
        updateCheckboxInput(session, input_name, value = TRUE)
      }
    })
    
    observeEvent(input$left_items_deselect_all, {
      current_items <- if (values$left_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$left_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$left_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$left_color_by)
      } else {
        character(0)
      }
      
      for (item in current_items) {
        input_name <- paste0("left_items_", gsub("[^a-zA-Z0-9]", "_", item))
        updateCheckboxInput(session, input_name, value = FALSE)
      }
    })
    
    # Toggle handlers for right plot
    observeEvent(input$right_items_select_all, {
      current_items <- if (values$right_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$right_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$right_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$right_color_by)
      } else {
        character(0)
      }
      
      for (item in current_items) {
        input_name <- paste0("right_items_", gsub("[^a-zA-Z0-9]", "_", item))
        updateCheckboxInput(session, input_name, value = TRUE)
      }
    })
    
    observeEvent(input$right_items_deselect_all, {
      current_items <- if (values$right_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$right_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$right_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$right_color_by)
      } else {
        character(0)
      }
      
      for (item in current_items) {
        input_name <- paste0("right_items_", gsub("[^a-zA-Z0-9]", "_", item))
        updateCheckboxInput(session, input_name, value = FALSE)
      }
    })
    
    # Item legend toggle UI for left plot
    output$leftLegendToggle <- renderUI({
      req(values$clustered_seurat)
      
      current_items <- if (values$left_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$left_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$left_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$left_color_by)
      } else {
        character(0)
      }
      
      if (length(current_items) == 0) {
        return(NULL)
      }
      
      toggle_title <- if (values$left_color_by == "cluster") {
        "Cluster Options"
      } else if (values$left_color_by == "sample") {
        "Sample Options"
      } else {
        paste(values$left_color_by, "Options")
      }
      
      create_legend_toggle_ui(ns, current_items, "left_items", toggle_title)
    })
    
    # Item legend toggle UI for right plot
    output$rightLegendToggle <- renderUI({
      req(values$clustered_seurat)
      
      current_items <- if (values$right_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$right_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$right_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$right_color_by)
      } else {
        character(0)
      }
      
      if (length(current_items) == 0) {
        return(NULL)
      }
      
      toggle_title <- if (values$right_color_by == "cluster") {
        "Cluster Options"
      } else if (values$right_color_by == "sample") {
        "Sample Options"
      } else {
        paste(values$right_color_by, "Options")
      }
      
      create_legend_toggle_ui(ns, current_items, "right_items", toggle_title)
    })
    
    # UMAP visualization section
    output$umapSection <- renderUI({
      req(values$clustering_done)
      
      div(
        style = "margin-top: 30px;",
        h3("UMAP Visualization"),
        
        # Side-by-side UMAP layout
        fluidRow(
          # Left UMAP
          column(6,
                 wellPanel(
                   div(
                     style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                     h4("UMAP 1"),
                     div(
                       style = "display: flex; gap: 10px;",
                       actionButton(ns("leftUMAP2D"), "2D", 
                                    class = if(values$left_umap_type == "2D") "btn-primary" else "btn-default"),
                       actionButton(ns("leftUMAP3D"), "3D", 
                                    class = if(values$left_umap_type == "3D") "btn-primary" else "btn-default")
                     )
                   ),
                   div(
                     style = "margin-bottom: 15px;",
                     selectInput(ns("leftColorBy"), "Color by:", 
                                 choices = get_coloring_options(values$clustered_seurat),
                                 selected = values$left_color_by)
                   ),
                   # Gene search for left UMAP (shown only when "gene" is selected)
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'gene'", ns("leftColorBy")),
                     div(
                       style = "margin-bottom: 15px;",
                       textInput(ns("leftGeneQuery"), "Search gene:", placeholder = "e.g. Sox10"),
                       actionButton(ns("searchLeftGene"), "Search", class = "btn-sm btn-primary"),
                       uiOutput(ns("leftGeneSelection"))
                     )
                   ),
                   # Legend toggle for clusters/samples (not for gene expression)
                   conditionalPanel(
                     condition = sprintf("input['%s'] != 'gene'", ns("leftColorBy")),
                     uiOutput(ns("leftLegendToggle"))
                   ),
                   # UMAP Visualization
                   div(
                     style = "text-align: right; margin-bottom: 5px;",
                     downloadButton(ns("downloadLeftUMAP"), "Save Plot", class = "btn-sm btn-success")
                   ),
                   # Plot output based on 2D/3D selection
                   if (values$left_umap_type == "2D") {
                     plotOutput(ns("leftUMAPPlot"), height = "500px")
                   } else {
                     plotlyOutput(ns("leftUMAPPlot3D"), height = "500px")
                   }
                 )
          ),
          
          # Right UMAP
          column(6,
                 wellPanel(
                   div(
                     style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                     h4("UMAP 2"),
                     div(
                       style = "display: flex; gap: 10px;",
                       actionButton(ns("rightUMAP2D"), "2D", 
                                    class = if(values$right_umap_type == "2D") "btn-primary" else "btn-default"),
                       actionButton(ns("rightUMAP3D"), "3D", 
                                    class = if(values$right_umap_type == "3D") "btn-primary" else "btn-default")
                     )
                   ),
                   div(
                     style = "margin-bottom: 15px;",
                     selectInput(ns("rightColorBy"), "Color by:", 
                                 choices = get_coloring_options(values$clustered_seurat),
                                 selected = values$right_color_by)
                   ),
                   # Gene search for right UMAP (shown only when "gene" is selected)
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'gene'", ns("rightColorBy")),
                     div(
                       style = "margin-bottom: 15px;",
                       textInput(ns("rightGeneQuery"), "Search gene:", placeholder = "e.g. Sox10"),
                       actionButton(ns("searchRightGene"), "Search", class = "btn-sm btn-primary"),
                       uiOutput(ns("rightGeneSelection"))
                     )
                   ),
                   # Legend toggle for clusters/samples (not for gene expression)
                   conditionalPanel(
                     condition = sprintf("input['%s'] != 'gene'", ns("rightColorBy")),
                     uiOutput(ns("rightLegendToggle"))
                   ),
                   # UMAP Visualization
                   div(
                     style = "text-align: right; margin-bottom: 5px;",
                     downloadButton(ns("downloadRightUMAP"), "Save Plot", class = "btn-sm btn-success")
                   ),
                   # Plot output based on 2D/3D selection
                   if (values$right_umap_type == "2D") {
                     plotOutput(ns("rightUMAPPlot"), height = "500px")
                   } else {
                     plotlyOutput(ns("rightUMAPPlot3D"), height = "500px")
                   }
                 )
          )
        )
      )
    })
    
    # Toggle UMAP types
    observeEvent(input$leftUMAP2D, {
      values$left_umap_type <- "2D"
    })
    
    observeEvent(input$leftUMAP3D, {
      values$left_umap_type <- "3D"
    })
    
    observeEvent(input$rightUMAP2D, {
      values$right_umap_type <- "2D"
    })
    
    observeEvent(input$rightUMAP3D, {
      values$right_umap_type <- "3D"
    })
    
    # Update color by selections
    observeEvent(input$leftColorBy, {
      values$left_color_by <- input$leftColorBy
    })
    
    observeEvent(input$rightColorBy, {
      values$right_color_by <- input$rightColorBy
    })
    
    # Gene search for left UMAP
    left_gene_search_results <- reactiveVal(NULL)
    left_selected_gene <- reactiveVal(NULL)
    
    observeEvent(input$searchLeftGene, {
      req(values$clustered_seurat, input$leftGeneQuery)
      left_gene_search_results(search_genes(values$clustered_seurat, input$leftGeneQuery))
    })
    
    output$leftGeneSelection <- renderUI({
      results <- left_gene_search_results()
      if (is.null(results) || nrow(results) == 0) {
        div(p("No matching genes found"))
      } else {
        div(
          style = "margin-top: 10px;",
          selectInput(ns("leftSelectedGene"), 
                      paste("Select from", nrow(results), "genes:"),
                      choices = setNames(results$ensembl_id, results$gene_symbol),
                      selected = results$ensembl_id[1])
        )
      }
    })
    
    observeEvent(input$leftSelectedGene, {
      left_selected_gene(input$leftSelectedGene)
    })
    
    # Gene search for right UMAP
    right_gene_search_results <- reactiveVal(NULL)
    right_selected_gene <- reactiveVal(NULL)
    
    observeEvent(input$searchRightGene, {
      req(values$clustered_seurat, input$rightGeneQuery)
      right_gene_search_results(search_genes(values$clustered_seurat, input$rightGeneQuery))
    })
    
    output$rightGeneSelection <- renderUI({
      results <- right_gene_search_results()
      if (is.null(results) || nrow(results) == 0) {
        div(p("No matching genes found"))
      } else {
        div(
          style = "margin-top: 10px;",
          selectInput(ns("rightSelectedGene"), 
                      paste("Select from", nrow(results), "genes:"),
                      choices = setNames(results$ensembl_id, results$gene_symbol),
                      selected = results$ensembl_id[1])
        )
      }
    })
    
    observeEvent(input$rightSelectedGene, {
      right_selected_gene(input$rightSelectedGene)
    })
    
    # Item toggle UI for left plot
    output$leftItemToggles <- renderUI({
      req(values$clustered_seurat)
      
      current_items <- if (values$left_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$left_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$left_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$left_color_by)
      } else {
        character(0)
      }
      
      if (length(current_items) == 0) {
        return(NULL)
      }
      
      toggle_title <- if (values$left_color_by == "cluster") {
        "Toggle Clusters"
      } else if (values$left_color_by == "sample") {
        "Toggle Samples"
      } else {
        paste("Toggle", values$left_color_by)
      }
      
      create_item_toggles(ns, current_items, "left_items", toggle_title)
    })
    
    # Item toggle UI for right plot
    output$rightItemToggles <- renderUI({
      req(values$clustered_seurat)
      
      current_items <- if (values$right_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$right_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$right_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$right_color_by)
      } else {
        character(0)
      }
      
      if (length(current_items) == 0) {
        return(NULL)
      }
      
      toggle_title <- if (values$right_color_by == "cluster") {
        "Toggle Clusters"
      } else if (values$right_color_by == "sample") {
        "Toggle Samples"
      } else {
        paste("Toggle", values$right_color_by)
      }
      
      create_item_toggles(ns, current_items, "right_items", toggle_title)
    })
    
    #' @title Get Active Toggle Items
    #' @description Gets the active items from toggle checkboxes
    #' @param input Shiny input object
    #' @param items Vector of all possible items
    #' @param input_id Base input ID for the toggle checkboxes
    #' @return Vector of active item values
    #' @keywords internal
    get_active_toggle_items <- function(input, items, input_id) {
      active_items <- character(0)
      
      for (item in items) {
        safe_id <- gsub("[^a-zA-Z0-9]", "_", item)
        input_name <- paste0(input_id, "_", safe_id)
        if (!is.null(input[[input_name]]) && input[[input_name]]) {
          active_items <- c(active_items, item)
        }
      }
      
      return(active_items)
    }
    
    # Get active items for left plot
    get_left_active_items <- reactive({
      req(values$clustered_seurat)
      
      if (values$left_color_by == "gene") {
        return(NULL)  # No filtering for gene expression
      }
      
      current_items <- if (values$left_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$left_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$left_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$left_color_by)
      } else {
        character(0)
      }
      
      get_active_toggle_items(input, current_items, "left_items")
    })
    
    # Get active items for right plot
    get_right_active_items <- reactive({
      req(values$clustered_seurat)
      
      if (values$right_color_by == "gene") {
        return(NULL)  # No filtering for gene expression
      }
      
      current_items <- if (values$right_color_by == "cluster") {
        get_unique_values("cluster")
      } else if (values$right_color_by == "sample") {
        get_unique_values("sample")
      } else if (values$right_color_by %in% colnames(values$clustered_seurat@meta.data)) {
        get_unique_values(values$right_color_by)
      } else {
        character(0)
      }
      
      get_active_toggle_items(input, current_items, "right_items")
    })
    
    # Left UMAP plot functions
    output$leftUMAPPlot <- renderPlot({
      req(values$clustered_seurat, values$left_umap_type == "2D")
      
      # Add this line to make the rendering reactive to the update trigger
      trigger <- values$left_plot_update_trigger
      
      # First, create a filtered Seurat object based on active clusters
      filtered_seurat <- values$clustered_seurat
      
      # Apply global sample filtering if there are active samples defined
      if (!is.null(values$filtered_samples) && length(values$filtered_samples) > 0) {
        # Filter to show only cells from active samples
        cells_to_keep <- filtered_seurat$sample %in% values$filtered_samples
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "No cells match the active sample selection") + 
                   theme_void())
        }
      }
      
      # Apply global condition filtering if there are active conditions defined
      if (!is.null(values$condition_column) && !is.null(values$filtered_conditions) && 
          length(values$filtered_conditions) > 0 && 
          values$condition_column %in% colnames(filtered_seurat@meta.data)) {
        
        # Filter to show only cells from active conditions
        cells_to_keep <- filtered_seurat@meta.data[[values$condition_column]] %in% values$filtered_conditions
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "No cells match the active condition selection") + 
                   theme_void())
        }
      }
      
      # Apply global cluster filtering if there are active clusters defined
      if (!is.null(values$filtered_clusters) && length(values$filtered_clusters) > 0) {
        # Filter to show only cells from active clusters
        cells_to_keep <- filtered_seurat$seurat_clusters %in% values$filtered_clusters
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        }
      }
      
      # Now get active items for the coloring-specific toggles
      active_items <- get_left_active_items()
      
      # Generate the plot based on coloring option with the already-filtered Seurat object
      if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
        create_2d_umap_plot(filtered_seurat, 
                            color_by = "gene", 
                            gene_id = left_selected_gene(), 
                            reduction = "umap2d")
      } else {
        create_2d_umap_plot(filtered_seurat, 
                            color_by = values$left_color_by, 
                            reduction = "umap2d",
                            active_items = active_items)
      }
    })
    
    output$leftUMAPPlot3D <- renderPlotly({
      req(values$clustered_seurat, values$left_umap_type == "3D")
      
      # Add this line to make the rendering reactive to the update trigger
      trigger <- values$left_plot_update_trigger
      
      # First, create a filtered Seurat object based on active clusters
      filtered_seurat <- values$clustered_seurat
      
      # Apply global sample filtering if there are active samples defined
      if (!is.null(values$filtered_samples) && length(values$filtered_samples) > 0) {
        # Filter to show only cells from active samples
        cells_to_keep <- filtered_seurat$sample %in% values$filtered_samples
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(plotly_empty() %>% 
                   layout(title = "No cells match the active sample selection"))
        }
      }
      
      # Apply global condition filtering if there are active conditions defined
      if (!is.null(values$condition_column) && !is.null(values$filtered_conditions) && 
          length(values$filtered_conditions) > 0 && 
          values$condition_column %in% colnames(filtered_seurat@meta.data)) {
        
        # Filter to show only cells from active conditions
        cells_to_keep <- filtered_seurat@meta.data[[values$condition_column]] %in% values$filtered_conditions
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(plotly_empty() %>% 
                   layout(title = "No cells match the active condition selection"))
        }
      }
      
      # Apply global cluster filtering if there are active clusters defined
      if (!is.null(values$filtered_clusters) && length(values$filtered_clusters) > 0) {
        # Filter to show only cells from active clusters
        cells_to_keep <- filtered_seurat$seurat_clusters %in% values$filtered_clusters
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(plotly_empty() %>% 
                   layout(title = "No cells match the active cluster selection"))
        }
      }
      
      # Now get active items for the coloring-specific toggles
      active_items <- get_left_active_items()
      
      # Generate the plot based on coloring option with the already-filtered Seurat object
      if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
        create_3d_gene_umap_plot(filtered_seurat, 
                                 gene_id = left_selected_gene(), 
                                 reduction = "umap3d")
      } else {
        create_3d_umap_plot(filtered_seurat, 
                            color_by = values$left_color_by, 
                            reduction = "umap3d",
                            active_items = active_items)
      }
    })
    
    output$rightUMAPPlot <- renderPlot({
      req(values$clustered_seurat, values$right_umap_type == "2D")
      
      # Add this line to make the rendering reactive to the update trigger
      trigger <- values$right_plot_update_trigger
      
      # First, create a filtered Seurat object based on active clusters
      filtered_seurat <- values$clustered_seurat
      
      # Apply global sample filtering if there are active samples defined
      if (!is.null(values$filtered_samples) && length(values$filtered_samples) > 0) {
        # Filter to show only cells from active samples
        cells_to_keep <- filtered_seurat$sample %in% values$filtered_samples
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "No cells match the active sample selection") + 
                   theme_void())
        }
      }
      
      # Apply global condition filtering if there are active conditions defined
      if (!is.null(values$condition_column) && !is.null(values$filtered_conditions) && 
          length(values$filtered_conditions) > 0 && 
          values$condition_column %in% colnames(filtered_seurat@meta.data)) {
        
        # Filter to show only cells from active conditions
        cells_to_keep <- filtered_seurat@meta.data[[values$condition_column]] %in% values$filtered_conditions
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "No cells match the active condition selection") + 
                   theme_void())
        }
      }
      
      # Apply global cluster filtering if there are active clusters defined
      if (!is.null(values$filtered_clusters) && length(values$filtered_clusters) > 0) {
        # Filter to show only cells from active clusters
        cells_to_keep <- filtered_seurat$seurat_clusters %in% values$filtered_clusters
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        }
      }
      
      # Now get active items for the coloring-specific toggles
      active_items <- get_right_active_items()
      
      # Generate the plot based on coloring option with the already-filtered Seurat object
      if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
        create_2d_umap_plot(filtered_seurat, 
                            color_by = "gene", 
                            gene_id = right_selected_gene(), 
                            reduction = "umap2d")
      } else {
        create_2d_umap_plot(filtered_seurat, 
                            color_by = values$right_color_by, 
                            reduction = "umap2d",
                            active_items = active_items)
      }
    })
    
    output$rightUMAPPlot3D <- renderPlotly({
      req(values$clustered_seurat, values$right_umap_type == "3D")
      
      # Add this line to make the rendering reactive to the update trigger
      trigger <- values$right_plot_update_trigger
      
      # First, create a filtered Seurat object based on active clusters
      filtered_seurat <- values$clustered_seurat
      
      # Apply global sample filtering if there are active samples defined
      if (!is.null(values$filtered_samples) && length(values$filtered_samples) > 0) {
        # Filter to show only cells from active samples
        cells_to_keep <- filtered_seurat$sample %in% values$filtered_samples
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(plotly_empty() %>% 
                   layout(title = "No cells match the active sample selection"))
        }
      }
      
      # Apply global condition filtering if there are active conditions defined
      if (!is.null(values$condition_column) && !is.null(values$filtered_conditions) && 
          length(values$filtered_conditions) > 0 && 
          values$condition_column %in% colnames(filtered_seurat@meta.data)) {
        
        # Filter to show only cells from active conditions
        cells_to_keep <- filtered_seurat@meta.data[[values$condition_column]] %in% values$filtered_conditions
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(plotly_empty() %>% 
                   layout(title = "No cells match the active condition selection"))
        }
      }
      
      # Apply global cluster filtering if there are active clusters defined
      if (!is.null(values$filtered_clusters) && length(values$filtered_clusters) > 0) {
        # Filter to show only cells from active clusters
        cells_to_keep <- filtered_seurat$seurat_clusters %in% values$filtered_clusters
        if (any(cells_to_keep)) {
          filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
        } else {
          return(plotly_empty() %>% 
                   layout(title = "No cells match the active cluster selection"))
        }
      }
      
      # Now get active items for the coloring-specific toggles
      active_items <- get_right_active_items()
      
      # Generate the plot based on coloring option with the already-filtered Seurat object
      if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
        create_3d_gene_umap_plot(filtered_seurat, 
                                 gene_id = right_selected_gene(), 
                                 reduction = "umap3d")
      } else {
        create_3d_umap_plot(filtered_seurat, 
                            color_by = values$right_color_by, 
                            reduction = "umap3d",
                            active_items = active_items)
      }
    })
    
    # Download handlers
    output$downloadLeftUMAP <- downloadHandler(
      filename = function() {
        suffix <- if (values$left_umap_type == "2D") "png" else "html"
        paste0("umap_left_", values$left_color_by, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", suffix)
      },
      content = function(file) {
        # First, create a filtered Seurat object based on active clusters
        filtered_seurat <- values$clustered_seurat
        
        # Apply global sample filtering if there are active samples defined
        if (!is.null(values$filtered_samples) && length(values$filtered_samples) > 0) {
          # Filter to show only cells from active samples
          cells_to_keep <- filtered_seurat$sample %in% values$filtered_samples
          if (any(cells_to_keep)) {
            filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
          }
        }
        
        # Apply global condition filtering if there are active conditions defined
        if (!is.null(values$condition_column) && !is.null(values$filtered_conditions) && 
            length(values$filtered_conditions) > 0 && 
            values$condition_column %in% colnames(filtered_seurat@meta.data)) {
          
          # Filter to show only cells from active conditions
          cells_to_keep <- filtered_seurat@meta.data[[values$condition_column]] %in% values$filtered_conditions
          if (any(cells_to_keep)) {
            filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
          }
        }
        
        # Apply global cluster filtering if there are active clusters defined
        if (!is.null(values$filtered_clusters) && length(values$filtered_clusters) > 0) {
          # Filter to show only cells from active clusters
          cells_to_keep <- filtered_seurat$seurat_clusters %in% values$filtered_clusters
          if (any(cells_to_keep)) {
            filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
          }
        }
        
        # Get active items for filtering
        active_items <- get_left_active_items()
        
        if (values$left_umap_type == "2D") {
          # Handle 2D plot
          p <- if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
            create_2d_umap_plot(filtered_seurat, 
                                color_by = "gene", 
                                gene_id = left_selected_gene(), 
                                reduction = "umap2d")
          } else {
            create_2d_umap_plot(filtered_seurat, 
                                color_by = values$left_color_by, 
                                reduction = "umap2d",
                                active_items = active_items)
          }
          save_ggplot(p, file)
        } else {
          # Handle 3D plot
          p <- if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
            create_3d_gene_umap_plot(filtered_seurat, 
                                     gene_id = left_selected_gene(), 
                                     reduction = "umap3d")
          } else {
            create_3d_umap_plot(filtered_seurat, 
                                color_by = values$left_color_by, 
                                reduction = "umap3d",
                                active_items = active_items)
          }
          save_plotly(p, file)
        }
      }
    )
    
    output$downloadRightUMAP <- downloadHandler(
      filename = function() {
        suffix <- if (values$right_umap_type == "2D") "png" else "html"
        paste0("umap_right_", values$right_color_by, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", suffix)
      },
      content = function(file) {
        # First, create a filtered Seurat object based on active clusters
        filtered_seurat <- values$clustered_seurat
        
        # Apply global sample filtering if there are active samples defined
        if (!is.null(values$filtered_samples) && length(values$filtered_samples) > 0) {
          # Filter to show only cells from active samples
          cells_to_keep <- filtered_seurat$sample %in% values$filtered_samples
          if (any(cells_to_keep)) {
            filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
          }
        }
        
        # Apply global condition filtering if there are active conditions defined
        if (!is.null(values$condition_column) && !is.null(values$filtered_conditions) && 
            length(values$filtered_conditions) > 0 && 
            values$condition_column %in% colnames(filtered_seurat@meta.data)) {
          
          # Filter to show only cells from active conditions
          cells_to_keep <- filtered_seurat@meta.data[[values$condition_column]] %in% values$filtered_conditions
          if (any(cells_to_keep)) {
            filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
          }
        }
        
        # Apply global cluster filtering if there are active clusters defined
        if (!is.null(values$filtered_clusters) && length(values$filtered_clusters) > 0) {
          # Filter to show only cells from active clusters
          cells_to_keep <- filtered_seurat$seurat_clusters %in% values$filtered_clusters
          if (any(cells_to_keep)) {
            filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
          }
        }
        
        # Get active items for filtering
        active_items <- get_right_active_items()
        
        if (values$right_umap_type == "2D") {
          # Handle 2D plot
          p <- if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
            create_2d_umap_plot(filtered_seurat, 
                                color_by = "gene", 
                                gene_id = right_selected_gene(), 
                                reduction = "umap2d")
          } else {
            create_2d_umap_plot(filtered_seurat, 
                                color_by = values$right_color_by, 
                                reduction = "umap2d",
                                active_items = active_items)
          }
          save_ggplot(p, file)
        } else {
          # Handle 3D plot
          p <- if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
            create_3d_gene_umap_plot(filtered_seurat, 
                                     gene_id = right_selected_gene(), 
                                     reduction = "umap3d")
          } else {
            create_3d_umap_plot(filtered_seurat, 
                                color_by = values$right_color_by, 
                                reduction = "umap3d",
                                active_items = active_items)
          }
          save_plotly(p, file)
        }
      }
    )
    
    # Download cluster statistics
    output$downloadClusterStats <- downloadHandler(
      filename = function() {
        paste0("cluster_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        req(values$clustered_seurat)
        
        # Get cluster statistics
        cluster_stats <- get_cluster_stats(values$clustered_seurat)
        write.csv(cluster_stats, file, row.names = FALSE)
      }
    )
    
    # Return the clustered Seurat object along with the state restoration method
    return(list(
      data = reactive({ values$clustered_seurat }),
      restoreState = restoreState
    ))
  })
}

#' @title Create Legend Toggle UI
#' @description Creates a collapsible legend toggle UI for controlling item visibility
#' @param ns Namespace function for creating input IDs
#' @param items Vector of items to toggle (clusters, samples, etc.)
#' @param input_id Base input ID to use for checkboxes
#' @param title String title for the toggle section
#' @param initial_state Logical whether all items should be selected initially
#' @return Shiny UI element for the legend toggle
#' @keywords internal
create_legend_toggle_ui <- function(ns, items, input_id, title = "Legend Options", initial_state = TRUE) {
  if (is.null(items) || length(items) == 0) {
    return(tags$div("No items available"))
  }
  
  # Create a collapsible panel
  dropdown_panel <- tags$div(
    class = "panel panel-default",
    style = "margin-bottom: 15px;",
    
    # Collapsible header
    tags$div(
      class = "panel-heading",
      style = "cursor: pointer; padding: 8px 15px;",
      `data-toggle` = "collapse",
      `data-target` = paste0("#", ns(paste0(input_id, "_collapse"))),
      tags$div(
        style = "display: flex; justify-content: space-between; align-items: center;",
        tags$span(tags$i(class = "fa fa-filter", style = "margin-right: 5px;"), title),
        tags$span(class = "caret")
      )
    ),
    
    # Collapsible content 
    tags$div(
      id = ns(paste0(input_id, "_collapse")),
      class = "panel-collapse collapse",
      tags$div(
        class = "panel-body",
        style = "padding: 10px;",
        
        # Toggle buttons
        tags$div(
          style = "display: flex; justify-content: space-between; margin-bottom: 10px;",
          actionButton(ns(paste0(input_id, "_select_all")), "Select All", 
                       class = "btn-sm btn-default", 
                       style = "flex: 1; margin-right: 5px;"),
          actionButton(ns(paste0(input_id, "_deselect_all")), "Deselect All", 
                       class = "btn-sm btn-default",
                       style = "flex: 1; margin-left: 5px;")
        ),
        
        # Item list with checkboxes
        tags$div(
          style = "max-height: 200px; overflow-y: auto; padding: 5px; border: 1px solid #ddd; border-radius: 4px;",
          lapply(items, function(item) {
            div(
              style = "margin-bottom: 5px;",
              checkboxInput(
                ns(paste0(input_id, "_", gsub("[^a-zA-Z0-9]", "_", item))),
                label = item,
                value = initial_state
              )
            )
          })
        )
      )
    )
  )
  
  return(dropdown_panel)
}

#' @title Create Item Toggles
#' @description Creates a set of checkboxes for toggling items
#' @param ns Namespace function for creating input IDs
#' @param items Vector of items to toggle (clusters, samples, etc.)
#' @param input_id Base input ID to use for checkboxes
#' @param title String title for the toggle section
#' @return Shiny UI element for the item toggles
#' @keywords internal
create_item_toggles <- function(ns, items, input_id, title = "Toggle Items") {
  if (is.null(items) || length(items) == 0) {
    return(tags$div("No items available"))
  }
  
  tags$div(
    tags$h5(title),
    tags$div(
      style = "display: flex; justify-content: space-between; margin-bottom: 10px;",
      actionButton(ns(paste0(input_id, "_select_all")), "Select All", 
                   class = "btn-sm btn-default"),
      actionButton(ns(paste0(input_id, "_deselect_all")), "Deselect All", 
                   class = "btn-sm btn-default")
    ),
    tags$div(
      style = "max-height: 300px; overflow-y: auto; border: 1px solid #ddd; padding: 8px;",
      lapply(items, function(item) {
        checkboxInput(
          ns(paste0(input_id, "_", gsub("[^a-zA-Z0-9]", "_", item))),
          label = item,
          value = TRUE
        )
      })
    )
  )
}