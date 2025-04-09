# R/modules/sample_management_module.R

#' @title Sample Management Module UI
#' @description Creates the UI for the sample management module which allows users to 
#'   toggle sample visibility and edit sample labels.
#' @param id The module ID
#' @return A Shiny UI element containing the sample management interface
#' @export
sampleManagementUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("sampleRows"))
  )
}

#' @title Sample Management Module Server
#' @description Server logic for the sample management module which handles sample
#'   activation/deactivation and label editing.
#' @param id The module ID
#' @param seurat_data Reactive expression containing the Seurat object
#' @return A list of reactive expressions for accessing sample labels and active status
#' @export
sampleManagementServer <- function(id, seurat_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Add loading locks to prevent competing updates
    is_loading_state <- reactiveVal(FALSE)
    ignore_input_changes <- reactiveVal(FALSE)
    
    # Add flags to control update flow and prevent cascades
    updating_from_checkbox <- reactiveVal(FALSE)
    updating_select_all <- reactiveVal(FALSE)
    
    # Use dedicated reactiveVals for stable storage
    stable_labels <- reactiveVal(list())
    stable_active <- reactiveVal(list())
    
    # Initialize state
    state <- reactiveValues(
      all_samples = NULL,
      last_update = NULL,
      temp_labels = NULL  # Add this to store temporary edits
    )
    
    # Flag to control UI updates
    inhibit_ui_updates <- reactiveVal(FALSE)
    
    # Add debouncing for UI updates to prevent multiple rapid redraws
    observe({
      # This will throttle UI updates when typing
      if (inhibit_ui_updates()) {
        # If inhibited, don't trigger UI updates for a short time
        invalidateLater(300) # Wait 300ms before allowing updates again
        inhibit_ui_updates(FALSE)
      }
    })
    
    # Create a deep copy of a list
    deep_copy <- function(list_obj) {
      if (is.null(list_obj)) return(list())
      result <- list()
      for (name in names(list_obj)) {
        result[[name]] <- list_obj[[name]]
      }
      return(result)
    }
    
    # Get samples from Seurat object and initialize state
    observe({
      # Skip if we're in loading state
      if (is_loading_state()) {
        return(NULL)
      }
      
      available_samples <- getAvailableSamples(seurat_data())
      
      # Skip if no samples available
      if (is.null(available_samples) || length(available_samples) == 0) {
        return(NULL)
      }
      
      # Update all_samples
      state$all_samples <- available_samples
      
      # Initialize labels if needed
      current_labels <- stable_labels()
      if (is.null(current_labels) || !all(available_samples %in% names(current_labels))) {
        new_labels <- initializeSampleLabels(available_samples, current_labels)
        stable_labels(new_labels)
        # Initialize temp_labels too
        state$temp_labels <- new_labels
      }
      
      # Initialize active status if needed
      current_active <- stable_active()
      if (is.null(current_active) || !all(available_samples %in% names(current_active))) {
        new_active <- initializeActiveStatus(available_samples, current_active)
        stable_active(new_active)
      }
    })
    
    # Setup UI for sample rows
    output$sampleRows <- renderUI({
      # Skip if no samples
      if (is.null(state$all_samples) || length(state$all_samples) == 0) {
        return(div(
          class = "alert alert-info",
          "No samples available. Please load data first."
        ))
      }
      
      # Use isolate here to prevent UI rebuilding during typing
      isolate({
        # Get current stable values
        current_labels <- stable_labels()
        current_active <- stable_active()
        current_temp <- state$temp_labels
        
        # Generate UI controls - NO SELECT ALL BUTTON HERE
        div(
          style = "max-height: 400px; overflow-y: auto; border: 1px solid #e3e3e3; padding: 5px; border-radius: 4px;",
          lapply(state$all_samples, function(sample) {
            safe_id <- make_safe_id(sample)
            
            # Get active status
            is_active <- if (!is.null(current_active) && sample %in% names(current_active)) {
              current_active[[sample]]
            } else {
              TRUE  # Default to active
            }
            
            # Get current label - Try to get from temp_labels first, then fall back to stable_labels
            current_label <- if (!is.null(current_temp) && sample %in% names(current_temp)) {
              current_temp[[sample]]
            } else if (!is.null(current_labels) && sample %in% names(current_labels)) {
              current_labels[[sample]]
            } else {
              sample  # Default to sample name
            }
            
            div(
              id = paste0("sample-row-", safe_id),
              style = paste0(
                "margin-bottom: 10px; padding: 8px; border-radius: 4px; ",
                if (is_active) "background-color: #f8f9fa;" else "background-color: #e9ecef; opacity: 0.8;"
              ),
              fluidRow(
                column(2,
                       tags$div(
                         style = "margin-top: 5px;",
                         tags$div(
                           class = "checkbox",
                           tags$label(
                             tags$input(
                               type = "checkbox",
                               id = ns(paste0("active_", safe_id)),
                               checked = if(is_active) "checked" else NULL,
                               class = "sample-management-checkbox",
                               `data-gsm` = sample  # Add data attribute for sample ID
                             )
                           )
                         )
                       )
                ),
                column(10, 
                       tags$div(
                         textInput(ns(paste0("label_", safe_id)),
                                   label = NULL,
                                   value = current_label,
                                   width = "100%")
                       )
                )
              )
            )
          })
        )
      })
    })
    
    # Track label changes with a dedicated observer
    observe({
      req(state$all_samples)
      
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      # Setup observers for each sample's label
      lapply(state$all_samples, function(sample) {
        safe_id <- make_safe_id(sample)
        input_id <- paste0("label_", safe_id)
        
        # Create observer for this specific input
        observeEvent(input[[input_id]], {
          # Skip if in loading state
          if (ignore_input_changes()) {
            return(NULL)
          }
          
          # Get current value
          current_value <- input[[input_id]]
          
          # Make sure temp_labels is initialized
          if (is.null(state$temp_labels)) {
            state$temp_labels <- list()
          }
          
          # Store in temp_labels without triggering UI refresh
          isolate({
            state$temp_labels[[sample]] <- current_value
          })
          
          # Don't update the UI or trigger other reactivity here
          # Just log the change for debugging
          isolate({
            print(paste("Temporary label for sample", sample, "set to", current_value))
          })
          
        }, ignoreInit = TRUE, ignoreNULL = TRUE)
      })
    })
    
    # Add JavaScript handler for individual checkboxes
    session$sendCustomMessage(type = "initializeSampleCheckboxes", message = list())
    
    # Handle selectAllSamples checkbox
    observeEvent(input$selectAllSamples, {
      # Skip if we're already updating from a checkbox change
      if (updating_from_checkbox()) {
        return(NULL)
      }
      
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      # Mark that we're updating from Select All
      updating_select_all(TRUE)
      
      req(state$all_samples)
      available_samples <- state$all_samples
      
      # Create a new active_samples list
      new_active <- setNames(
        rep(input$selectAllSamples, length(available_samples)),
        available_samples
      )
      
      # Block individual checkbox observers
      ignore_input_changes(TRUE)
      
      # Update stable storage with our copy
      stable_active(deep_copy(new_active))
      
      # Update individual checkboxes using JavaScript
      selected_samples <- if(input$selectAllSamples) available_samples else character(0)
      session$sendCustomMessage(
        type = "updateSampleCheckboxes",
        message = list(samples = selected_samples)
      )
      
      # Now it's safe to allow changes - slight delay to ensure UI updates first
      shinyjs::delay(200, {
        ignore_input_changes(FALSE)
        updating_select_all(FALSE)
        
        # Trigger state update
        state$last_update <- Sys.time()
      })
    }, ignoreInit = TRUE, priority = 10)  # Higher priority
    
    # Handler for individual sample checkbox updates from JavaScript
    observeEvent(input$sample_checkbox_changed, {
      # Skip if in loading state or during Select All update
      if (ignore_input_changes() || updating_select_all()) {
        return(NULL)
      }
      
      # Extract the sample ID and checked state from the input
      sample_id <- input$sample_checkbox_changed$sample
      is_checked <- input$sample_checkbox_changed$checked
      
      # Mark that we're updating from a checkbox
      updating_from_checkbox(TRUE)
      
      # Get current active status - create a complete copy to avoid reference issues
      current_active <- deep_copy(stable_active())
      
      # Get previous value with default
      old_value <- if (!is.null(current_active) && sample_id %in% names(current_active)) {
        current_active[[sample_id]]
      } else {
        TRUE
      }
      
      # Check if value actually changed
      if (old_value != is_checked) {
        # Update our deep copy with the new value
        current_active[[sample_id]] <- is_checked
        
        # Update stable active status
        stable_active(current_active)
        
        # Update "Select All" checkbox state
        all_active <- all(unlist(current_active))
        updateCheckboxInput(session, "selectAllSamples", value = all_active)
        
        # Signal update
        state$last_update <- Sys.time()
      }
      
      # Clear updating flag after a delay
      shinyjs::delay(100, {
        updating_from_checkbox(FALSE)
      })
    })
    
    # Update labels when the update button is clicked
    observeEvent(input$updateSampleLabels, {
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      req(state$all_samples)
      
      # Only proceed if we have temporary labels to apply
      if (!is.null(state$temp_labels) && length(state$temp_labels) > 0) {
        # Get current stable labels
        current_labels <- stable_labels()
        if (is.null(current_labels)) {
          current_labels <- list()
        }
        
        # Merge in temporary labels
        for (sample in names(state$temp_labels)) {
          current_labels[[sample]] <- state$temp_labels[[sample]]
        }
        
        # Update stable labels atomically
        stable_labels(current_labels)
        
        # Trigger last update for reactivity
        state$last_update <- Sys.time()
        
        showNotification("Sample labels updated", type = "message")
      } else {
        showNotification("No label changes to save", type = "message")
      }
    })
    
    # Function to get state for saving
    getFullState = function() {
      list(
        sample_labels = stable_labels(),
        active_samples = stable_active(),
        all_samples = state$all_samples
      )
    }
    
    # Function to restore state
    setFullState = function(saved_state) {
      if (!is.null(saved_state)) {
        # Set loading flags
        is_loading_state(TRUE)
        ignore_input_changes(TRUE)
        
        print("Setting sample management full state")
        
        # Block updates during loading
        shinyjs::delay(100, {
          # Restore all_samples if provided
          if (!is.null(saved_state$all_samples) && length(saved_state$all_samples) > 0) {
            state$all_samples <- saved_state$all_samples
          }
          
          # Restore labels atomically
          if (!is.null(saved_state$sample_labels)) {
            print(paste("Restoring", length(saved_state$sample_labels), "sample labels"))
            stable_labels(saved_state$sample_labels)
            
            # Also set temp_labels to match
            state$temp_labels <- saved_state$sample_labels
            
            # Print all labels being restored
            for (sample in names(saved_state$sample_labels)) {
              print(paste("Label for sample", sample, "=", saved_state$sample_labels[[sample]]))
            }
          }
          
          # Restore active status
          if (!is.null(saved_state$active_samples)) {
            print(paste("Restoring", length(saved_state$active_samples), "active status values"))
            stable_active(saved_state$active_samples)
            
            # Print active status
            for (sample in names(saved_state$active_samples)) {
              print(paste("Sample", sample, "active =", saved_state$active_samples[[sample]]))
            }
          }
          
          # Update UI after a short delay
          shinyjs::delay(200, {
            # Update sample labels
            if (!is.null(saved_state$sample_labels)) {
              for (sample in names(saved_state$sample_labels)) {
                safe_id <- make_safe_id(sample)
                updateTextInput(session, paste0("label_", safe_id), 
                                value = saved_state$sample_labels[[sample]])
              }
            }
            
            # Update checkboxes using JavaScript for better consistency
            if (!is.null(saved_state$active_samples)) {
              active_samples <- names(saved_state$active_samples[unlist(saved_state$active_samples) == TRUE])
              session$sendCustomMessage(
                type = "updateSampleCheckboxes",
                message = list(samples = active_samples)
              )
              
              # Update "Select All" checkbox
              all_active <- all(unlist(saved_state$active_samples))
              updateCheckboxInput(session, "selectAllSamples", value = all_active)
            }
            
            # Force UI refresh
            state$last_update <- Sys.time()
            
            # Release locks after delay
            shinyjs::delay(500, {
              ignore_input_changes(FALSE)
              shinyjs::delay(200, {
                is_loading_state(FALSE)
                print("Sample state loading complete")
              })
            })
          })
        })
      }
    }
    
    # Return reactive expressions
    list(
      getSampleLabels = reactive({ stable_labels() }),
      getActiveStatus = reactive({ stable_active() }),
      getActiveSampleIds = reactive({
        active_samples <- stable_active()
        if (is.null(active_samples)) return(NULL)
        names(active_samples[unlist(active_samples) == TRUE])
      }),
      getActiveSampleList = reactive({
        active_samples <- stable_active()
        if (is.null(active_samples)) return(NULL)
        names(active_samples[unlist(active_samples) == TRUE])
      }),
      updateLabels = function(new_labels) {
        stable_labels(new_labels)
      },
      updateActiveStatus = function(new_active) {
        stable_active(new_active)
        
        # Update UI to reflect the new active status
        active_samples <- names(new_active[unlist(new_active) == TRUE])
        session$sendCustomMessage(
          type = "updateSampleCheckboxes",
          message = list(samples = active_samples)
        )
      },
      getFullState = getFullState,
      setFullState = setFullState,
      updateFromButton = function() {
        # Set UI update inhibitor
        inhibit_ui_updates(TRUE)
        
        # Only proceed if we have temporary labels to apply
        if (!is.null(state$temp_labels) && length(state$temp_labels) > 0) {
          # Start with current stable labels
          current_labels <- stable_labels()
          if (is.null(current_labels)) {
            current_labels <- list()
          }
          
          # Merge in temporary labels
          for (sample in names(state$temp_labels)) {
            current_labels[[sample]] <- state$temp_labels[[sample]]
          }
          
          # Apply the update
          stable_labels(current_labels)
          
          # Clear temp labels to avoid duplicate application
          state$temp_labels <- list()
          
          # Trigger update in other modules
          state$last_update <- Sys.time()
          
          showNotification("Sample labels updated", type = "message")
        } else {
          showNotification("No label changes to save", type = "message")
        }
        
        # Release UI inhibitor after a short delay
        shinyjs::delay(200, {
          inhibit_ui_updates(FALSE)
        })
      }
    )
  })
}

#' @title Get Available Samples
#' @description Helper function that extracts unique sample IDs from a Seurat object.
#' @param seurat_obj Seurat object 
#' @return Vector of unique sample IDs, or NULL if no samples are found
#' @keywords internal
getAvailableSamples <- function(seurat_obj) {
  if (is.null(seurat_obj) || !"sample" %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }
  unique(seurat_obj$sample)
}

#' @title Initialize Sample Labels
#' @description Creates default labels for samples, preserving any existing labels
#'   from the current_labels parameter.
#' @param samples Vector of sample IDs
#' @param current_labels Optional named vector of existing sample labels
#' @return A named vector mapping sample IDs to their labels
#' @keywords internal
initializeSampleLabels <- function(samples, current_labels = NULL) {
  new_labels <- setNames(
    samples,  # Use sample IDs as default labels
    samples
  )
  
  # Merge with existing labels if they exist
  if (!is.null(current_labels)) {
    existing_samples <- names(current_labels)
    for (sample in existing_samples) {
      if (sample %in% names(new_labels)) {
        new_labels[sample] <- current_labels[sample]
      }
    }
  }
  
  return(new_labels)
}

#' @title Make Safe ID
#' @description Creates a safe HTML/Shiny ID from a sample ID by replacing
#'   special characters with underscores and ensuring it starts with a letter.
#' @param id String value to convert to a safe ID
#' @return String containing a safe ID for use in HTML/Shiny
#' @keywords internal
make_safe_id <- function(id) {
  # Replace special characters with underscores
  safe_id <- gsub("[^a-zA-Z0-9]", "_", id)
  
  # Ensure it starts with a letter (Shiny input IDs must start with a letter)
  if (!grepl("^[a-zA-Z]", safe_id)) {
    safe_id <- paste0("s_", safe_id)
  }
  
  return(safe_id)
}

#' @title Initialize Active Status
#' @description Creates a named vector indicating which samples are active (TRUE/FALSE),
#'   with the default being all samples active. Preserves existing statuses if provided.
#' @param samples Vector of sample IDs
#' @param current_active Optional named vector of existing active statuses
#' @return A named logical vector indicating which samples are active
#' @keywords internal
initializeActiveStatus <- function(samples, current_active = NULL) {
  new_active <- setNames(
    rep(TRUE, length(samples)),
    samples
  )
  
  # Merge with existing active status if it exists
  if (!is.null(current_active)) {
    existing_samples <- names(current_active)
    for (sample in existing_samples) {
      if (sample %in% names(new_active)) {
        new_active[sample] <- current_active[sample]
      }
    }
  }
  
  return(new_active)
}

#' @title Filter By Samples
#' @description Filters a Seurat object to include only cells from specified samples.
#' @param seurat_obj Seurat object
#' @param active_samples Vector of sample names to keep
#' @return Filtered Seurat object
#' @export
filterBySamples <- function(seurat_obj, active_samples) {
  # Check if we have a valid Seurat object and active samples
  if (is.null(seurat_obj) || is.null(active_samples) || length(active_samples) == 0) {
    return(seurat_obj)
  }
  
  # Check if sample column exists
  if (!"sample" %in% colnames(seurat_obj@meta.data)) {
    warning("Sample column not found in Seurat object metadata")
    return(seurat_obj)
  }
  
  # Get cells that belong to active samples
  cells_to_keep <- seurat_obj$sample %in% active_samples
  
  # Only subset if necessary
  if (all(cells_to_keep)) {
    return(seurat_obj)  # All cells are already active
  } else if (!any(cells_to_keep)) {
    warning("No cells match the active samples selection")
    return(seurat_obj)  # Return original to avoid empty object
  }
  
  # Subset the Seurat object
  subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
}