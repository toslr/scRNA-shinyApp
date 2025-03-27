# R/modules/dimension_reduction_module.R

# Source utility files
# These will be sourced in app.R, but are listed here for reference
# source("R/modules/dimension_reduction_utils/dimred_computation.R")
# source("R/modules/dimension_reduction_utils/dimred_visualization.R")

#' UI for dimension reduction module
#'
#' @param id Character ID for namespacing
#' @return UI elements for dimension reduction
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

#' Server function for dimension reduction module
#'
#' @param id Character ID for namespacing
#' @param processed_seurat Reactive expression returning Seurat object with PCA
#' @param sample_management Optional sample management module
#' @param condition_management Optional condition management module
#' @return Reactive expression returning processed Seurat object
dimensionReductionServer <- function(id, processed_seurat, sample_management = NULL, condition_management = NULL) {
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
      right_color_by = "cluster"
    )
    
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
    
    # Helper to get unique values for a metadata column
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
    
    # Helper function to get active toggle items
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
    
    # Left UMAP plot functions
    output$leftUMAPPlot <- renderPlot({
      req(values$clustered_seurat, values$left_umap_type == "2D")
      
      # Get active items for filtering
      active_items <- get_left_active_items()
      
      # Generate the plot based on coloring option
      if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
        create_2d_umap_plot(values$clustered_seurat, 
                            color_by = "gene", 
                            gene_id = left_selected_gene(), 
                            reduction = "umap2d")
      } else {
        create_2d_umap_plot(values$clustered_seurat, 
                            color_by = values$left_color_by, 
                            reduction = "umap2d",
                            active_items = active_items)
      }
    })
    
    output$leftUMAPPlot3D <- renderPlotly({
      req(values$clustered_seurat, values$left_umap_type == "3D")
      
      # Get active items for filtering
      active_items <- get_left_active_items()
      
      # Generate the plot based on coloring option
      if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
        create_3d_gene_umap_plot(values$clustered_seurat, 
                                 gene_id = left_selected_gene(), 
                                 reduction = "umap3d")
      } else {
        create_3d_umap_plot(values$clustered_seurat, 
                            color_by = values$left_color_by, 
                            reduction = "umap3d",
                            active_items = active_items)
      }
    })
    
    # Right UMAP plot functions
    output$rightUMAPPlot <- renderPlot({
      req(values$clustered_seurat, values$right_umap_type == "2D")
      
      # Get active items for filtering
      active_items <- get_right_active_items()
      
      # Generate the plot based on coloring option
      if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
        create_2d_umap_plot(values$clustered_seurat, 
                            color_by = "gene", 
                            gene_id = right_selected_gene(), 
                            reduction = "umap2d")
      } else {
        create_2d_umap_plot(values$clustered_seurat, 
                            color_by = values$right_color_by, 
                            reduction = "umap2d",
                            active_items = active_items)
      }
    })
    
    output$rightUMAPPlot3D <- renderPlotly({
      req(values$clustered_seurat, values$right_umap_type == "3D")
      
      # Get active items for filtering
      active_items <- get_right_active_items()
      
      # Generate the plot based on coloring option
      if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
        create_3d_gene_umap_plot(values$clustered_seurat, 
                                 gene_id = right_selected_gene(), 
                                 reduction = "umap3d")
      } else {
        create_3d_umap_plot(values$clustered_seurat, 
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
        # Get active items for filtering
        active_items <- get_left_active_items()
        
        if (values$left_umap_type == "2D") {
          # Handle 2D plot
          p <- if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
            create_2d_umap_plot(values$clustered_seurat, 
                                color_by = "gene", 
                                gene_id = left_selected_gene(), 
                                reduction = "umap2d")
          } else {
            create_2d_umap_plot(values$clustered_seurat, 
                                color_by = values$left_color_by, 
                                reduction = "umap2d",
                                active_items = active_items)
          }
          save_ggplot(p, file)
        } else {
          # Handle 3D plot
          p <- if (values$left_color_by == "gene" && !is.null(left_selected_gene())) {
            create_3d_gene_umap_plot(values$clustered_seurat, 
                                     gene_id = left_selected_gene(), 
                                     reduction = "umap3d")
          } else {
            create_3d_umap_plot(values$clustered_seurat, 
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
        # Get active items for filtering
        active_items <- get_right_active_items()
        
        if (values$right_umap_type == "2D") {
          # Handle 2D plot
          p <- if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
            create_2d_umap_plot(values$clustered_seurat, 
                                color_by = "gene", 
                                gene_id = right_selected_gene(), 
                                reduction = "umap2d")
          } else {
            create_2d_umap_plot(values$clustered_seurat, 
                                color_by = values$right_color_by, 
                                reduction = "umap2d",
                                active_items = active_items)
          }
          save_ggplot(p, file)
        } else {
          # Handle 3D plot
          p <- if (values$right_color_by == "gene" && !is.null(right_selected_gene())) {
            create_3d_gene_umap_plot(values$clustered_seurat, 
                                     gene_id = right_selected_gene(), 
                                     reduction = "umap3d")
          } else {
            create_3d_umap_plot(values$clustered_seurat, 
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
    
    # Return the clustered Seurat object
    return(reactive({ values$clustered_seurat }))
  })
}