# R/ui.R

#' @title Build Application UI
#' @description Constructs the complete user interface for the single-cell RNA-seq analysis application,
#'   including the sidebar navigation, module UIs, and main content area.
#' @return A Shiny UI object containing the complete application interface
#' @export
buildUI <- function() {
  fluidPage(
    useShinyjs(),
    includeCSS("www/styles.css"),
    includeScript("www/script.js"),
    
    # Top bar with title and navigation
    div(id = "topbar",
        h1(id = "app-title", "Single-Cell RNA Analysis"),
        # Add the navigation panel here
        div(id = "top-navigation",
            uiOutput("metroNavigation")
        )
    ),
    
    # Main layout
    div(
      # Sidebar
      div(id = "sidebar",
          # Metadata section
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
          
          # Data Input section
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
          
          # Save/Load Analysis section
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
          
          # Sample Management section
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
          
          # Condition Management section
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
          
          # Cluster Management section
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
          
          # Copyright footer
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
      
      # Scrollable main content area
      div(id = "main-content",
          uiOutput("metadataSection"),
          uiOutput("qcSection"),
          uiOutput("dimredSection"),
          uiOutput("deSection")
      )
    )
  )
}