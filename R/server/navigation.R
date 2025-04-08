# R/navigation.R (modified)

setupNavigation <- function(output, steps_completed) {
  # Metro-style navigation
  output$metroNavigation <- renderUI({
    # Define our modules/steps in order
    steps <- list(
      list(id = "metadata-section", name = "Metadata", completed = steps_completed$metadata),
      list(id = "qc-section", name = "Quality Control", completed = steps_completed$qc, 
           available = steps_completed$data_input),
      list(id = "dimred-section", name = "Dimension Reduction", completed = steps_completed$dimred, 
           available = steps_completed$qc),
      list(id = "de-section", name = "Differential Expression", completed = steps_completed$de, 
           available = steps_completed$clustering)
    )
    
    # Create the metro line UI
    div(class = "metro-navigation",
        div(class = "metro-line",
            # Add each station
            lapply(1:length(steps), function(i) {
              step <- steps[[i]]
              
              # Determine the status class
              status_class <- if (!is.null(step$available) && !step$available) {
                "disabled"
              } else if (step$completed) {
                "completed"
              } else {
                "available"
              }
              
              # Create station
              div(class = paste("metro-station", status_class),
                  # Add clickable functionality
                  if (!is.null(step$available) && !step$available) {
                    # Non-clickable if not available
                    div(class = "metro-station-content",
                        div(class = "metro-dot"),
                        div(class = "metro-label", step$name)
                    )
                  } else {
                    # Clickable if available
                    tags$a(
                      href = "javascript:void(0)",
                      onclick = paste0("scrollToSection('", step$id, "')"),
                      class = "metro-station-content",
                      div(class = "metro-dot"),
                      div(class = "metro-label", 
                          step$name,
                          if(step$completed) tags$span(class="metro-check", "✓"))
                    )
                  }
              )
            })
        )
    )
  })
  
  # Keep the original task list empty
  output$taskList <- renderUI({
    NULL
  })
}