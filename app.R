library(shiny)
library(patchwork)
library(shinyFiles)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scCustomize)
library(shinyjs)

source("app/data_input_module.R")
source("app/qc_module.R")
source("app/dimension_reduction_module.R")
source("app/de_analysis_module.R")


# Custom JS for smooth scrolling
jscode <- "
function scrollToSection(sectionId) {
  document.getElementById(sectionId).scrollIntoView({ 
    behavior: 'smooth', 
    block: 'start' 
  });
}
"

ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$script(HTML(jscode)),
    # CSS for fixed topbar, sidebar and scrollable main panel
    tags$style(HTML("
      #topbar {
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        height: 60px;
        background-color: #2c3e50;
        color: white;
        z-index: 1000;
        padding: 0 20px;
        display: flex;
        align-items: center;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      #sidebar {
        position: fixed;
        top: 60px;  /* Start below topbar */
        left: 0;
        width: 20%;
        height: calc(100vh - 60px);  /* Full height minus topbar */
        overflow-y: auto;
        padding: 15px;
        background-color: #f5f5f5;
        border-right: 1px solid #e3e3e3;
      }
      #main-content {
        margin-left: 20%;
        margin-top: 60px;  /* Start below topbar */
        padding: 15px;
        width: 80%;
      }
      .nav-pills > li > a {
        padding: 8px 15px;
        margin: 2px 0;
      }
      .nav-pills > li.active > a {
        background-color: #337ab7;
        color: white;
      }
      .badge {
        margin-left: 5px;
        background-color: #5cb85c;
      }
      #app-title {
        font-size: 24px;
        font-weight: 500;
        margin: 0;
      }
    "))
  ),
  
  # Top bar with title
  div(id = "topbar",
      h1(id = "app-title", "Single-Cell RNA Analysis")
  ),
  
  # Main layout
  div(
    # Fixed sidebar
    div(id = "sidebar",
        dataInputUI("dataInput"),
        tags$hr(),
        h4("Navigation"),
        div(
          id = "navigation-panel",
          uiOutput("taskList")
        )
    ),
    # Scrollable main content
    div(id = "main-content",
        uiOutput("qcSection"),
        uiOutput("dimredSection"),
        uiOutput("deSection")
    )
  )
)

server <- function(input, output, session) {
  # Chain the reactive values through the modules
  seurat_data <- dataInputServer("dataInput")
  processed_seurat <- qcServer("qc", seurat_data)
  clustered_seurat <- dimensionReductionServer("dimRed", processed_seurat)
  de_module <- deAnalysisServer("de", clustered_seurat)  # Now returns a list
  
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
        h3("Quality Control"),
        qcUI("qc")
    )
  })
  
  output$dimredSection <- renderUI({
    req(processed_seurat())
    div(id = "dimred-section",
        h3("Dimension Reduction"),
        dimensionReductionUI("dimRed")
    )
  })
  
  output$deSection <- renderUI({
    req(clustered_seurat())
    req("seurat_clusters" %in% colnames(clustered_seurat()@meta.data))
    div(id = "de-section",
        h3("Differential Expression"),
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

shinyApp(ui = ui, server = server)