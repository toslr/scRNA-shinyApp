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
      # Sidebar
      div(id = "sidebar",
          div(class = "sidebar-section",
              metadataUI("metadata")
          ),
          tags$hr(),
          div(class = "sidebar-section",
              dataInputUI("dataInput")
          ),
          tags$hr(),
          div(class = "sidebar-section",
              saveLoadUI("saveLoad")
          ),
          tags$hr(),
          h4("Navigation"),
          div(
            id = "navigation-panel",
            uiOutput("taskList")
          ),
          div(
            id= "copyright",
            tags$hr(),
            p("© Tom Soulaire",
              style = "text-align: left; color: #666; font-size: 0.9em; margin-bottom: 5px;"),
            p("Zuchero Lab",
              style = "text-align: left; color: #666; font-size: 0.9em; margin-bottom: 5px;"),
            p("Stanford University",
              style = "text-align: left; color: #666; font-size: 0.9em; margin-bottom: 5px;")
          )
      ),
      
      # Scrollable main content
      div(id = "main-content",
          uiOutput("metadataSection"),
          uiOutput("qcSection"),
          uiOutput("dimredSection"),
          uiOutput("deSection")
      )
    )
  )
}