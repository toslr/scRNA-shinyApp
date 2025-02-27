# R/modules/de_analysis_module/ui_components.R

createClusterControls <- function(ns, available_clusters, current_temp_labels, current_active) {
  tagList(
    # Individual cluster controls
    lapply(available_clusters, function(cluster) {
      cluster_key <- as.character(cluster)
      is_active <- if (cluster_key %in% names(current_active)) {
        current_active[[cluster_key]]
      } else {
        TRUE  # Default to active if not found
      }
      
      current_label <- if (cluster_key %in% names(current_temp_labels)) {
        current_temp_labels[[cluster_key]]
      } else {
        paste("Cluster", cluster)
      }
      
      div(
        style = paste0(
          "margin-bottom: 10px; padding: 8px; border-radius: 4px; ",
          if (is_active) "background-color: #f8f9fa;" else "background-color: #e9ecef; opacity: 0.8;"
        ),
        fluidRow(
          column(1,
                 checkboxInput(ns(paste0("active_", cluster)), 
                               label = NULL,
                               value = is_active)
          ),
          column(4, 
                 tags$div(
                   tags$label(
                     style = "margin-bottom: 5px; font-weight: normal;",
                     paste("Cluster", cluster)
                   ),
                   textInput(ns(paste0("label_", cluster)),
                             label = NULL,
                             value = current_label)
                 )
          )
        )
      )
    })
  )
}

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