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
    
    # Initialize state
    # Using reactiveValues for all state elements for more stable reactivity
    state <- reactiveValues(
      sample_labels = NULL,
      active_samples = NULL,
      all_samples = NULL,
      input_values = list()
    )
    
    # Get samples from Seurat object and initialize state
    observe({
      available_samples <- getAvailableSamples(seurat_data())
      
      # Skip if no samples available
      if (is.null(available_samples) || length(available_samples) == 0) {
        return(NULL)
      }
      
      # Update all_samples
      state$all_samples <- available_samples
      
      # Initialize labels if needed
      if (is.null(state$sample_labels) || !all(available_samples %in% names(state$sample_labels))) {
        new_labels <- initializeSampleLabels(available_samples, state$sample_labels)
        state$sample_labels <- new_labels
      }
      
      # Initialize active status if needed
      if (is.null(state$active_samples) || !all(available_samples %in% names(state$active_samples))) {
        state$active_samples <- initializeActiveStatus(available_samples, state$active_samples)
      }
    })
    
    # Setup UI for sample rows with input tracking
    output$sampleRows <- renderUI({
      # Skip if no samples
      if (is.null(state$all_samples) || length(state$all_samples) == 0) {
        return(div(
          class = "alert alert-info",
          "No samples available. Please load data first."
        ))
      }
      
      # Generate UI controls
      tagList(
        lapply(state$all_samples, function(sample) {
          safe_id <- make_safe_id(sample)
          is_active <- if (!is.null(state$active_samples) && sample %in% names(state$active_samples)) {
            state$active_samples[[sample]]
          } else {
            TRUE  # Default to active
          }
          
          # Get current label with fallback options
          current_label <- NULL
          
          # First check if we have a stored input value
          input_id <- paste0("label_", safe_id)
          if (!is.null(state$input_values[[input_id]])) {
            current_label <- state$input_values[[input_id]]
          } 
          # Then check saved labels
          else if (!is.null(state$sample_labels) && sample %in% names(state$sample_labels)) {
            current_label <- state$sample_labels[[sample]]
          } 
          # Fallback to sample name
          else {
            current_label <- sample
          }
          
          div(
            style = paste0(
              "margin-bottom: 10px; padding: 8px; border-radius: 4px; ",
              if (is_active) "background-color: #f8f9fa;" else "background-color: #e9ecef; opacity: 0.8;"
            ),
            fluidRow(
              column(2,
                     checkboxInput(ns(paste0("active_", safe_id)), 
                                   label = "",
                                   value = is_active)
              ),
              column(10, 
                     tags$div(
                       textInput(ns(paste0("label_", safe_id)),
                                 label = NULL,
                                 value = current_label)
                     )
              )
            )
          )
        })
      )
    })
    
    # Track input values for each sample's label
    observe({
      req(state$all_samples)
      
      lapply(state$all_samples, function(sample) {
        safe_id <- make_safe_id(sample)
        input_id <- paste0("label_", safe_id)
        
        # Update input_values when input changes
        if (!is.null(input[[input_id]])) {
          # This isolate prevents a circular dependency
          isolate({
            state$input_values[[input_id]] <- input[[input_id]]
          })
        }
      })
    })
    
    # Handle Select All checkbox
    observeEvent(input$selectAllSamples, {
      req(state$all_samples)
      
      # Create a new active_samples list with all samples set to the checkbox value
      new_active <- setNames(
        rep(input$selectAllSamples, length(state$all_samples)),
        state$all_samples
      )
      
      # Update state
      state$active_samples <- new_active
      
      # Update individual checkboxes to match
      for (sample in state$all_samples) {
        safe_id <- make_safe_id(sample)
        updateCheckboxInput(session, paste0("active_", safe_id), value = input$selectAllSamples)
      }
    })
    
    # Handle individual sample activation
    lapply(1:100, function(i) { # Use a large number to handle many potential samples
      observeEvent(input[[paste0("active_", i)]], {
        # Find which sample this corresponds to
        req(state$all_samples)
        
        matching_samples <- sapply(state$all_samples, function(sample) {
          make_safe_id(sample) == as.character(i)
        })
        
        sample_match <- state$all_samples[matching_samples]
        
        if (length(sample_match) == 1) {
          # Update the active status for this sample
          state$active_samples[[sample_match]] <- input[[paste0("active_", i)]]
          
          # Update "Select All" checkbox based on all samples
          all_selected <- all(unlist(state$active_samples))
          updateCheckboxInput(session, "selectAllSamples", value = all_selected)
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)
    })
    
    # Update labels when the update button is clicked
    observeEvent(input$updateAllLabels, {
      req(state$all_samples)
      
      # Collect all input values
      for (sample in state$all_samples) {
        safe_id <- make_safe_id(sample)
        input_id <- paste0("label_", safe_id)
        
        # Only update if we have a value
        if (!is.null(input[[input_id]])) {
          state$sample_labels[[sample]] <- input[[input_id]]
          state$input_values[[input_id]] <- input[[input_id]]
        }
      }
      
      showNotification("Sample labels updated", type = "message")
    })
    
    # Observe individual checkbox changes
    observe({
      req(state$all_samples)
      
      for (sample in state$all_samples) {
        local({
          local_sample <- sample
          safe_id <- make_safe_id(local_sample)
          input_id <- paste0("active_", safe_id)
          
          if (!is.null(input[[input_id]])) {
            observeEvent(input[[input_id]], {
              # Update active status for this sample
              state$active_samples[[local_sample]] <- input[[input_id]]
              
              # Update "Select All" checkbox
              all_active <- all(unlist(state$active_samples))
              updateCheckboxInput(session, "selectAllSamples", value = all_active)
            }, ignoreInit = TRUE)
          }
        })
      }
    })
    
    # Return reactive expressions
    list(
      getSampleLabels = reactive({ state$sample_labels }),
      getActiveStatus = reactive({ state$active_samples }),
      getActiveSampleIds = reactive({
        if (is.null(state$active_samples)) return(NULL)
        names(state$active_samples[unlist(state$active_samples) == TRUE])
      }),
      getFullState = function() {
        list(
          sample_labels = state$sample_labels,
          active_samples = state$active_samples,
          all_samples = state$all_samples
        )
      },
      setFullState = function(saved_state) {
        if (!is.null(saved_state)) {
          if (!is.null(saved_state$sample_labels)) {
            state$sample_labels <- saved_state$sample_labels
          }
          if (!is.null(saved_state$active_samples)) {
            state$active_samples <- saved_state$active_samples
          }
          # Force UI refresh
          state$all_samples <- state$all_samples
        }
      },
      getActiveSampleList = reactive({
        if (is.null(state$active_samples)) return(NULL)
        names(state$active_samples[unlist(state$active_samples) == TRUE])
      }),
      updateLabels = function(new_labels) {
        state$sample_labels <- new_labels
        for (sample in names(new_labels)) {
          safe_id <- make_safe_id(sample)
          input_id <- paste0("label_", safe_id)
          state$input_values[[input_id]] <- new_labels[[sample]]
        }
      },
      updateActiveStatus = function(new_active) {
        state$active_samples <- new_active
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