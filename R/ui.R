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
          div(class = "sidebar-section collapsible-section",
              div(class = "section-header",
                  h4("Metadata", style = "display: inline;"),
                  tags$button(
                    class = "btn btn-link btn-collapse",
                    tags$i(class = "fa fa-chevron-down")
                  )
              ),
              div(class = "section-content",
                  metadataUI("metadata")
              )
          ),
          
          div(class = "sidebar-section collapsible-section",
              div(class = "section-header",
                  h4("Data Input", style = "display: inline;"),
                  tags$button(
                    class = "btn btn-link btn-collapse",
                    tags$i(class = "fa fa-chevron-down")
                  )
              ),
              div(class = "section-content",
                  dataInputUI("dataInput")
              )
          ),
          
          div(class = "sidebar-section collapsible-section",
              div(class = "section-header",
                  h4("Save/Load Analysis", style = "display: inline;"),
                  tags$button(
                    class = "btn btn-link btn-collapse",
                    tags$i(class = "fa fa-chevron-down")
                  )
              ),
              div(class = "section-content",
                  saveLoadUI("saveLoad")
              )
          ),
          
          div(class = "sidebar-section collapsible-section",
              div(class = "section-header",
                  h4("Navigation", style = "display: inline;"),
                  tags$button(
                    class = "btn btn-link btn-collapse",
                    tags$i(class = "fa fa-chevron-down")
                  )
              ),
              div(class = "section-content",
                  div(id = "navigation-panel",
                      uiOutput("taskList")
                  )
              )
          ),
          
          div(class = "sidebar-section collapsible-section",
              div(class = "section-header",
                  h4("Sample Management", style = "display: inline;"),
                  tags$button(
                    class = "btn btn-link btn-collapse",
                    tags$i(class = "fa fa-chevron-down")
                  )
              ),
              div(class = "section-content",
                  uiOutput("sampleControls")
              )
          ),
          
          div(class = "sidebar-section collapsible-section",
              div(class = "section-header",
                  h4("Condition Management", style = "display: inline;"),
                  tags$button(
                    class = "btn btn-link btn-collapse",
                    tags$i(class = "fa fa-chevron-down")
                  )
              ),
              div(class = "section-content",
                  uiOutput("conditionControls")
              )
          ),
          
          div(class = "sidebar-section collapsible-section",
              div(class = "section-header",
                  h4("Cluster Management", style = "display: inline;"),
                  tags$button(
                    class = "btn btn-link btn-collapse",
                    tags$i(class = "fa fa-chevron-down")
                  )
              ),
              div(class = "section-content",
                  uiOutput("clusterControls")
              )
          ),
          
          div(
            id= "copyright",
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