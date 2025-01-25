# app.R

# Load required libraries
library(shiny)
library(patchwork)
library(shinyFiles)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scCustomize)
library(shinyjs)

# Source all module files
source("R/ui.R")
source("R/modules/data_input_module.R")
source("R/modules/qc_module.R")
source("R/modules/dimension_reduction_module.R")
source("R/modules/de_analysis_module.R")


server <- function(input, output, session) {
  # Chain the reactive values through the modules
  seurat_data <- dataInputServer("dataInput")
  processed_seurat <- qcServer("qc", seurat_data)
  clustered_seurat <- dimensionReductionServer("dimRed", processed_seurat)
  de_module <- deAnalysisServer("de", clustered_seurat)
  
  # Create reactive values to track completion of each step
  steps_completed <- reactiveValues(
    data_input = FALSE,
    qc = FALSE,
    dimred = FALSE,
    clustering = FALSE,
    de = FALSE
  )
  
  # Data input completion
  observe({
    if (!is.null(seurat_data())) {
      steps_completed$data_input <- TRUE
    } else {
      steps_completed$data_input <- FALSE
    }
  })
  
  # QC completion
  observe({
    if (!is.null(processed_seurat())) {
      steps_completed$qc <- TRUE
    } else {
      steps_completed$qc <- FALSE
    }
  })
  
  # Dimension reduction completion
  observe({
    if (!is.null(clustered_seurat()) && 
        "umap" %in% names(clustered_seurat()@reductions)) {
      steps_completed$dimred <- TRUE
    } else {
      steps_completed$dimred <- FALSE
    }
  })
  
  # Clustering completion
  observe({
    if (!is.null(clustered_seurat()) && 
        "seurat_clusters" %in% colnames(clustered_seurat()@meta.data)) {
      steps_completed$clustering <- TRUE
    } else {
      steps_completed$clustering <- FALSE
    }
  })
  
  # DE completion
  observe({
    if (!is.null(de_module$status()) && 
        de_module$status() == "completed") {
      steps_completed$de <- TRUE
    } else {
      steps_completed$de <- FALSE
    }
  })
  
  # Conditional section renders
  output$qcSection <- renderUI({
    req(seurat_data())
    div(id = "qc-section",
        h3(class = "section-header", "Quality Control"),
        qcUI("qc")
    )
  })
  
  output$dimredSection <- renderUI({
    req(processed_seurat())
    div(id = "dimred-section",
        h3(class = "section-header", "Dimension Reduction"),
        dimensionReductionUI("dimRed")
    )
  })
  
  output$deSection <- renderUI({
    req(clustered_seurat())
    req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
    div(id = "de-section",
        h3(class = "section-header", "Differential Expression"),
        deAnalysisUI("de")
    )
  })
  
  # Dynamic task list
  output$taskList <- renderUI({
    tasks <- tags$ul(class = "nav nav-pills nav-stacked",
                     style = "margin-top: 10px;")
    
    # Show QC link if data is loaded
    if (steps_completed$data_input) {
      tasks <- tagAppendChild(tasks,
                              tags$li(class = if(steps_completed$qc) "active" else "",
                                      tags$a(href = "javascript:void(0)",
                                             onclick = "scrollToSection('qc-section')",
                                             "Quality Control",
                                             if(steps_completed$qc) tags$span(class="badge", "✓"))
                              ))
    }
    
    # Show dimension reduction link if QC is done
    if (steps_completed$qc) {
      tasks <- tagAppendChild(tasks,
                              tags$li(class = if(steps_completed$dimred) "active" else "",
                                      tags$a(href = "javascript:void(0)",
                                             onclick = "scrollToSection('dimred-section')",
                                             "Dimension Reduction",
                                             if(steps_completed$dimred) tags$span(class="badge", "✓"))
                              ))
    }
    
    # Show DE link if clustering is done
    if (steps_completed$clustering) {
      tasks <- tagAppendChild(tasks,
                              tags$li(class = if(steps_completed$de) "active" else "",
                                      tags$a(href = "javascript:void(0)",
                                             onclick = "scrollToSection('de-section')",
                                             "Differential Expression",
                                             if(steps_completed$de) tags$span(class="badge", "✓"))
                              ))
    }
    
    tasks
  })
}

shinyApp(ui = buildUI(), server = server)