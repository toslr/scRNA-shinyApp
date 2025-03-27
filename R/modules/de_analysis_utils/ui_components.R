#' @title Create Analysis UI
#' @description Creates the user interface components for differential expression analysis,
#'   including one vs all analysis, one vs one analysis, and general heatmap controls.
#' @param ns Namespace function for the module
#' @param cluster_choices Named vector of cluster choices for select inputs
#' @return A tagList containing the DE analysis UI components
#' @export
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