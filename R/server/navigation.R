# R/navigation.R

setupNavigation <- function(output, steps_completed) {
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