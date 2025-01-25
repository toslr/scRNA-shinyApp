# R/ui.R

buildUI <- function() {
  fluidPage(
    useShinyjs(),
    includeCSS("www/styles.css"),
    includeScript("www/script.js"),
    
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
}