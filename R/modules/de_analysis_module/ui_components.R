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
    DT::dataTableOutput(ns("deTable")),
    uiOutput(ns("heatmapControls")),
    plotOutput(ns("heatmapPlot"), height = "600px")
  )
}