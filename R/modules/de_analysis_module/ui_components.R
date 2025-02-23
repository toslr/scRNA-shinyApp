# R/modules/de_analysis_module/ui_components.R

createClusterControls <- function(ns, available_clusters, current_temp_labels, current_active) {
  tagList(
    lapply(available_clusters, function(cluster) {
      div(style = "margin-bottom: 10px;",
          fluidRow(
            column(1,
                   checkboxInput(ns(paste0("active_", cluster)), 
                                 label = NULL,
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
}

createAnalysisUI <- function(ns, cluster_choices) {
  tagList(
    fluidRow(
      column(4,
             wellPanel(
               h4("One vs All Analysis"),
               selectInput(ns("targetClusterAll"), 
                           "Select cluster to compare against all others:", 
                           choices = cluster_choices,
                           selected = NULL),
               actionButton(ns("runDEAll"), "Run One vs All DE")
             )
      ),
      column(4,
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
      ),
      column(4,
           wellPanel(
             h4("General Cluster Map"),
             numericInput(ns("genesPerCluster"),
                          "Top genes per cluster:",
                          value=5,
                          min=1,
                          max=50),
             actionButton(ns("runGeneralHeatmap"), "Generate General Heatmap")
           )
      )
    ),
    
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