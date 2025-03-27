# R/modules/sample_management_module.R

sampleManagementUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("sampleRows"))
  )
}

sampleManagementServer <- function(id, seurat_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize state
    state <- list(
      sample_labels = reactiveVal(NULL),
      active_samples = reactiveVal(NULL),
      temp_labels = reactiveVal(NULL),
      label_inputs = reactiveValues()
    )
    
    # Get samples from Seurat object
    observe({
      available_samples <- getAvailableSamples(seurat_data())
      
      # Skip if no samples available
      if (is.null(available_samples) || length(available_samples) == 0) {
        return(NULL)
      }
      
      # Initialize or update labels if needed
      current_labels <- state$sample_labels()
      if (is.null(current_labels) || 
          !all(available_samples %in% names(current_labels))) {
        
        new_labels <- initializeSampleLabels(available_samples, current_labels)
        
        state$sample_labels(new_labels)
        state$temp_labels(new_labels)
      }
      
      # Initialize active status if needed
      current_active <- state$active_samples()
      if (is.null(current_active) || 
          !all(available_samples %in% names(current_active))) {
        state$active_samples(initializeActiveStatus(available_samples, current_active))
      }
    })
    
    # Handle initial population when samples change
    observeEvent(seurat_data(), {
      available_samples <- getAvailableSamples(seurat_data())
      
      if (!is.null(available_samples) && length(available_samples) > 0) {
        # Initialize new labels and active status
        new_labels <- initializeSampleLabels(available_samples)
        state$sample_labels(new_labels)
        state$temp_labels(new_labels)
        
        new_active <- initializeActiveStatus(available_samples)
        state$active_samples(new_active)
      }
    }, ignoreInit = TRUE)
    
    # Setup UI for sample rows
    output$sampleRows <- renderUI({
      available_samples <- getAvailableSamples(seurat_data())
      
      # Show message if no samples
      if (is.null(available_samples) || length(available_samples) == 0) {
        return(div(
          class = "alert alert-info",
          "No samples available. Please load data first."
        ))
      }
      
      # Get current labels and active status
      current_temp_labels <- state$temp_labels()
      current_active <- state$active_samples()
      
      # Return UI
      createSampleControls(ns, available_samples, current_temp_labels, current_active)
    })
    
    # Handle Select All checkbox
    observeEvent(input$selectAllSamples, {
      available_samples <- getAvailableSamples(seurat_data())
      if (is.null(available_samples) || length(available_samples) == 0) {
        return(NULL)
      }
      
      current_active <- state$active_samples()
      for (sample in available_samples) {
        current_active[sample] <- input$selectAllSamples
      }
      state$active_samples(current_active)
    })
    
    # Handle active status updates for individual samples
    observe({
      available_samples <- getAvailableSamples(seurat_data())
      
      # Skip empty samples
      if (is.null(available_samples) || length(available_samples) == 0) {
        return(NULL)
      }
      
      lapply(available_samples, function(sample) {
        input_id <- paste0("active_", make_safe_id(sample))
        
        # Only observe if this input exists
        if (!is.null(input[[input_id]])) {
          observeEvent(input[[input_id]], {
            current_active <- state$active_samples()
            if (is.null(current_active)) {
              return(NULL)
            }
            
            # Check if this sample exists in our active_samples
            if (sample %in% names(current_active)) {
              current_active[sample] <- input[[input_id]]
              state$active_samples(current_active)
              
              # Update select all checkbox based on all sample checkboxes
              all_selected <- all(unlist(current_active))
              updateCheckboxInput(session, "selectAllSamples", value = all_selected)
            }
          }, ignoreInit = TRUE)
        }
      })
    })
    
    # Update label_inputs when text changes
    observe({
      available_samples <- getAvailableSamples(seurat_data())
      if (is.null(available_samples) || length(available_samples) == 0) {
        return(NULL)
      }
      
      for(sample in available_samples) {
        local({
          local_sample <- sample
          input_id <- paste0("label_", make_safe_id(local_sample))
          
          # Only observe if this input exists
          if (!is.null(input[[input_id]])) {
            observeEvent(input[[input_id]], {
              state$label_inputs[[input_id]] <- input[[input_id]]
            }, ignoreInit = TRUE)
          }
        })
      }
    })
    
    # Update temp_labels only when update button is clicked
    observeEvent(input$updateAllLabels, {
      available_samples <- getAvailableSamples(seurat_data())
      if (is.null(available_samples) || length(available_samples) == 0) {
        return(NULL)
      }
      
      current_temp <- state$temp_labels()
      
      for(sample in available_samples) {
        input_id <- paste0("label_", make_safe_id(sample))
        if(!is.null(state$label_inputs[[input_id]])) {
          current_temp[sample] <- state$label_inputs[[input_id]]
        }
      }
      
      state$temp_labels(current_temp)
      state$sample_labels(current_temp)
      
      showNotification("Sample labels updated", type = "message")
    })
    
    # Return reactive expressions
    list(
      getSampleLabels = reactive({ state$sample_labels() }),
      getActiveStatus = reactive({ state$active_samples() }),
      getActiveSampleIds = reactive({
        current_active <- state$active_samples()
        if (is.null(current_active)) return(NULL)
        names(current_active[current_active == TRUE])
      }),
      getActiveSampleList = reactive({
        current_active <- state$active_samples()
        if (is.null(current_active)) return(NULL)
        names(current_active[current_active == TRUE])
      }),
      updateLabels = function(new_labels) {
        state$sample_labels(new_labels)
        state$temp_labels(new_labels)
      },
      updateActiveStatus = function(new_active) {
        state$active_samples(new_active)
      }
    )
  })
}

# Helper function: Create UI controls for each sample
createSampleControls <- function(ns, available_samples, current_temp_labels, current_active) {
  tagList(
    # Individual sample controls
    lapply(available_samples, function(sample) {
      is_active <- if (sample %in% names(current_active)) {
        current_active[[sample]]
      } else {
        TRUE  # Default to active if not found
      }
      
      current_label <- if (sample %in% names(current_temp_labels)) {
        current_temp_labels[[sample]]
      } else {
        sample
      }
      
      safe_id <- make_safe_id(sample)
      
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
}

# Helper function: Get available samples from Seurat object
getAvailableSamples <- function(seurat_obj) {
  if (is.null(seurat_obj) || !"sample" %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }
  unique(seurat_obj$sample)
}

# Initialize sample labels
initializeSampleLabels <- function(samples, current_labels = NULL) {
  new_labels <- setNames(
    samples,  # Use sample IDs as default labels
    samples
  )
  
  # Merge with existing labels if they exist
  if (!is.null(current_labels)) {
    existing_samples <- names(current_labels)
    new_labels[existing_samples] <- current_labels[existing_samples]
  }
  
  return(new_labels)
}

# Helper function to make safe HTML/Shiny IDs from GSM IDs or other sample identifiers
make_safe_id <- function(id) {
  # Replace special characters with underscores
  safe_id <- gsub("[^a-zA-Z0-9]", "_", id)
  
  # Ensure it starts with a letter (Shiny input IDs must start with a letter)
  if (!grepl("^[a-zA-Z]", safe_id)) {
    safe_id <- paste0("s_", safe_id)
  }
  
  return(safe_id)
}

# Filter Seurat object to include only active samples
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