#' @title Condition Management Module UI
#' @description Creates the UI for the condition management module which allows users to 
#'   select metadata columns as conditions and toggle condition values.
#' @param id The module ID
#' @return A Shiny UI element containing the condition management interface
#' @export
conditionManagementUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      class = "condition-management-container",
      uiOutput(ns("conditionSelector")),
      uiOutput(ns("conditionRows"))
    )
  )
}

#' @title Condition Management Module Server
#' @description Server logic for the condition management module which allows selecting
#'   metadata columns as conditions and filtering data by condition values.
#' @param id The module ID
#' @param seurat_data Reactive expression containing the Seurat object
#' @param metadata_module Metadata module instance for accessing metadata information
#' @return A list of reactive expressions for accessing condition information and active status
#' @export
conditionManagementServer <- function(id, seurat_data, metadata_module) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize state using reactiveValues for better state management
    state <- reactiveValues(
      condition_column = NULL,
      condition_labels = NULL,
      active_conditions = NULL,
      temp_labels = NULL,
      available_columns = NULL,
      input_values = list()  # Added to track input values like in sample_management
    )
    
    # Get available metadata columns
    observe({
      # Skip if no Seurat object
      if (is.null(seurat_data())) return(NULL)
      
      # Get metadata columns
      meta_cols <- colnames(seurat_data()@meta.data)
      
      # Filter out standard columns we don't want to show
      exclude_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA", 
                           "seurat_clusters", "percent.mt", "sample", "title")
      filtered_cols <- setdiff(meta_cols, exclude_columns)
      
      # Set available columns
      state$available_columns <- filtered_cols
      
      # Set default column selection if needed
      if (is.null(state$condition_column) && length(filtered_cols) > 0) {
        # Prioritize condition-like columns
        priority_cols <- grep("condition|treatment|group|genotype|timepoint", 
                              filtered_cols, value = TRUE, ignore.case = TRUE)
        
        if (length(priority_cols) > 0) {
          state$condition_column <- priority_cols[1]
        } else if (length(filtered_cols) > 0) {
          state$condition_column <- filtered_cols[1]
        }
        
        # Initialize conditions for the selected column
        updateConditionState(state$condition_column)
      }
    })
    
    # Render condition column selector
    output$conditionSelector <- renderUI({
      req(state$available_columns)
      available_cols <- state$available_columns
      
      if (length(available_cols) == 0) {
        return(div(
          class = "alert alert-info",
          "No metadata columns available. Please load data first."
        ))
      }
      
      # Find the current selected column
      current_column <- state$condition_column
      
      selectInput(ns("conditionColumn"), 
                  "Select condition column:", 
                  choices = available_cols,
                  selected = current_column)
    })
    
    # Update condition state when dropdown selection changes
    observeEvent(input$conditionColumn, {
      state$condition_column <- input$conditionColumn
      updateConditionState(input$conditionColumn)
    })
    
    #' @title Update Condition State
    #' @description Initialize or update condition management state based on the selected column
    #' @param column The selected metadata column
    #' @keywords internal
    updateConditionState <- function(column) {
      if (is.null(column) || is.null(seurat_data())) return(NULL)
      
      # Get unique condition values
      if (column %in% colnames(seurat_data()@meta.data)) {
        condition_values <- unique(as.character(seurat_data()@meta.data[[column]]))
        
        # Skip if no values
        if (length(condition_values) == 0) return(NULL)
        
        # Initialize labels and active status
        state$condition_labels <- initializeConditionLabels(condition_values, state$condition_labels)
        state$temp_labels <- state$condition_labels
        state$active_conditions <- initializeActiveStatus(condition_values, state$active_conditions)
      }
    }
    
    #' @title Get Available Conditions
    #' @description Get unique condition values from the selected column
    #' @return Vector of unique condition values or NULL if none available
    #' @keywords internal
    getAvailableConditions <- function() {
      column <- state$condition_column
      if (is.null(column) || is.null(seurat_data())) return(NULL)
      
      if (column %in% colnames(seurat_data()@meta.data)) {
        unique(as.character(seurat_data()@meta.data[[column]]))
      } else {
        NULL
      }
    }
    
    # Setup UI for condition rows
    output$conditionRows <- renderUI({
      req(state$condition_column)
      column <- state$condition_column
      
      available_conditions <- getAvailableConditions()
      
      # Show message if no conditions
      if (is.null(available_conditions) || length(available_conditions) == 0) {
        return(div(
          class = "alert alert-info",
          paste0("No condition values found in column '", column, "'. ",
                 "Please select a different column.")
        ))
      }
      
      # Get current labels and active status
      current_temp_labels <- state$temp_labels
      current_active <- state$active_conditions
      
      # Calculate if all conditions are currently selected
      all_selected <- !is.null(current_active) && all(unlist(current_active))
      
      # Return UI
      tagList(
        div(
          style = "display: flex; align-items: center; margin: 15px 0 10px 0;",
          checkboxInput(ns("selectAllConditions"), "Select All Conditions", 
                        value = all_selected),
          tags$div(style = "margin-left: 8px;",
                   actionButton(ns("updateAllLabels"), "Save Labels", 
                                class = "btn-sm btn-primary")
          )
        ),
        div(
          style = "max-height: 400px; overflow-y: auto; border: 1px solid #e3e3e3; padding: 5px; border-radius: 4px;",
          createConditionControls(ns, available_conditions, current_temp_labels, current_active)
        )
      )
    })
    
    # Track input values for each condition's checkbox
    observe({
      req(getAvailableConditions())
      available_conditions <- getAvailableConditions()
      
      # Monitor all checkboxes
      lapply(available_conditions, function(condition) {
        safe_id <- make_safe_id(condition)
        input_id <- paste0("active_", safe_id)
        
        # Update input_values when checkbox changes
        if (!is.null(input[[input_id]])) {
          # This isolate prevents a circular dependency
          isolate({
            if (!identical(state$input_values[[input_id]], input[[input_id]])) {
              state$input_values[[input_id]] <- input[[input_id]]
              
              # Also update the active_conditions directly
              state$active_conditions[[condition]] <- input[[input_id]]
            }
          })
        }
      })
    })
    
    # Handle selectAllConditions checkbox
    observeEvent(input$selectAllConditions, {
      req(getAvailableConditions())
      available_conditions <- getAvailableConditions()
      
      # Create a new active_conditions list with all conditions set to the checkbox value
      new_active <- setNames(
        rep(input$selectAllConditions, length(available_conditions)),
        available_conditions
      )
      
      # Update state
      state$active_conditions <- new_active
      
      # Update individual checkboxes to match
      for (condition in available_conditions) {
        safe_id <- make_safe_id(condition)
        updateCheckboxInput(session, paste0("active_", safe_id), value = input$selectAllConditions)
      }
    })
    
    # Track status of "Select All" checkbox based on individual selections
    observe({
      req(getAvailableConditions())
      available_conditions <- getAvailableConditions()
      
      # Skip if not initialized
      if (is.null(state$active_conditions) || length(state$active_conditions) == 0) {
        return(NULL)
      }
      
      # Check if all conditions are active
      all_active <- all(unlist(state$active_conditions))
      
      # Update the "Select All" checkbox if needed
      if (!is.null(input$selectAllConditions) && input$selectAllConditions != all_active) {
        updateCheckboxInput(session, "selectAllConditions", value = all_active)
      }
    })
    
    # Add JavaScript to handle checkbox changes
    observe({
      # Generate JavaScript for all checkboxes
      req(getAvailableConditions())
      available_conditions <- getAvailableConditions()
      
      lapply(available_conditions, function(condition) {
        safe_id <- make_safe_id(condition)
        checkbox_id <- paste0(ns("active_"), safe_id)
        
        # Add 'condition-management-checkbox' class to each checkbox
        shinyjs::addClass(id = paste0("active_", safe_id), class = "condition-management-checkbox")
      })
    })
    
    # Handle label inputs
    observe({
      req(getAvailableConditions())
      available_conditions <- getAvailableConditions()
      
      for(condition in available_conditions) {
        local({
          local_condition <- condition
          safe_id <- make_safe_id(local_condition)
          input_id <- paste0("label_", safe_id)
          
          # Create observer for this specific label input
          if (!is.null(input[[input_id]])) {
            observeEvent(input[[input_id]], {
              # Store the input value in temp_labels
              state$temp_labels[[local_condition]] <- input[[input_id]]
            }, ignoreInit = TRUE)
          }
        })
      }
    })
    
    # Update permanent labels when update button is clicked
    observeEvent(input$updateAllLabels, {
      # Update the permanent labels from the temporary ones
      state$condition_labels <- state$temp_labels
      showNotification("Condition labels updated", type = "message")
    })
    
    # Return reactive expressions
    list(
      getConditionColumn = reactive({ state$condition_column }),
      getConditionLabels = reactive({ state$condition_labels }),
      getActiveStatus = reactive({ state$active_conditions }),
      getActiveConditions = reactive({
        conditions <- state$active_conditions
        if (is.null(conditions)) return(NULL)
        names(conditions[conditions == TRUE])
      }),
      getActiveConditionList = reactive({
        conditions <- state$active_conditions
        if (is.null(conditions)) return(NULL)
        names(conditions[conditions == TRUE])
      }),
      updateLabels = function(new_labels) {
        state$condition_labels <- new_labels
        state$temp_labels <- new_labels
      },
      updateActiveStatus = function(new_active) {
        state$active_conditions <- new_active
      }
    )
  })
}

#' @title Create Condition Controls
#' @description Helper function that creates UI controls for each condition,
#'   including active/inactive toggle and label editing.
#' @param ns Namespace function
#' @param available_conditions Vector of available condition values
#' @param current_temp_labels Named vector of current condition labels
#' @param current_active Named vector of current active status (TRUE/FALSE)
#' @return A UI element with controls for each condition
#' @keywords internal
createConditionControls <- function(ns, available_conditions, current_temp_labels, current_active) {
  tagList(
    # Individual condition controls
    lapply(seq_along(available_conditions), function(i) {
      condition <- available_conditions[i]
      is_active <- if (!is.null(current_active) && condition %in% names(current_active)) {
        current_active[[condition]]
      } else {
        TRUE  # Default to active if not found
      }
      
      current_label <- if (!is.null(current_temp_labels) && condition %in% names(current_temp_labels)) {
        current_temp_labels[[condition]]
      } else {
        condition
      }
      
      safe_id <- make_safe_id(condition)
      
      div(
        id = paste0("condition-row-", i),
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
                         class = "condition-management-checkbox",
                         checked = if(is_active) "checked" else NULL
                       )
                     )
                   )
                 )
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

#' @title Make Safe ID
#' @description Creates a safe HTML/Shiny ID from a condition value by replacing
#'   special characters with underscores and ensuring it starts with a letter.
#' @param id String value to convert to a safe ID
#' @return String containing a safe ID for use in HTML/Shiny
#' @keywords internal
make_safe_id <- function(id) {
  # First, convert to character to handle numeric values
  id_str <- as.character(id)
  
  # Replace special characters with underscores
  safe_id <- gsub("[^a-zA-Z0-9]", "_", id_str)
  
  # Ensure it starts with a letter (Shiny input IDs must start with a letter)
  if (!grepl("^[a-zA-Z]", safe_id)) {
    safe_id <- paste0("c_", safe_id)
  }
  
  return(safe_id)
}

#' @title Initialize Condition Labels
#' @description Creates default labels for conditions, preserving any existing labels
#'   from the current_labels parameter.
#' @param conditions Vector of condition values
#' @param current_labels Optional named vector of existing condition labels
#' @return A named vector mapping condition values to their labels
#' @keywords internal
initializeConditionLabels <- function(conditions, current_labels = NULL) {
  new_labels <- setNames(
    conditions,  # Use condition values as default labels
    conditions
  )
  
  # Merge with existing labels if they exist
  if (!is.null(current_labels)) {
    existing_conditions <- names(current_labels)
    for (condition in existing_conditions) {
      if (condition %in% names(new_labels)) {
        new_labels[condition] <- current_labels[condition]
      }
    }
  }
  
  return(new_labels)
}

#' @title Initialize Active Status
#' @description Creates a named vector indicating which conditions are active (TRUE/FALSE),
#'   with the default being all conditions active. Preserves existing statuses if provided.
#' @param conditions Vector of condition values
#' @param current_active Optional named vector of existing active statuses
#' @return A named logical vector indicating which conditions are active
#' @keywords internal
initializeActiveStatus <- function(conditions, current_active = NULL) {
  new_active <- setNames(
    rep(TRUE, length(conditions)),
    conditions
  )
  
  if (!is.null(current_active)) {
    existing_conditions <- names(current_active)
    for (condition in existing_conditions) {
      if (condition %in% names(new_active)) {
        new_active[condition] <- current_active[condition]
      }
    }
  }
  
  return(new_active)
}

#' @title Filter By Conditions
#' @description Filters a Seurat object to include only cells that match active conditions.
#' @param seurat_obj Seurat object to filter
#' @param condition_column String name of the metadata column to use for filtering
#' @param active_conditions Vector of condition values to keep
#' @return Filtered Seurat object
#' @export
filterByConditions <- function(seurat_obj, condition_column, active_conditions) {
  # Check if we have valid inputs
  if (is.null(seurat_obj) || is.null(condition_column) || 
      is.null(active_conditions) || length(active_conditions) == 0) {
    return(seurat_obj)
  }
  
  # Check if condition column exists
  if (!(condition_column %in% colnames(seurat_obj@meta.data))) {
    warning(paste("Condition column", condition_column, "not found in Seurat object metadata"))
    return(seurat_obj)
  }
  
  # Get cells that match active conditions
  cells_to_keep <- seurat_obj@meta.data[[condition_column]] %in% active_conditions
  
  # Only subset if necessary
  if (all(cells_to_keep)) {
    return(seurat_obj)  # All cells are already active
  } else if (!any(cells_to_keep)) {
    warning("No cells match the active conditions selection")
    return(seurat_obj)  # Return original to avoid empty object
  }
  
  # Subset the Seurat object
  subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
}