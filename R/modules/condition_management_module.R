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
    
    # Use flags to control update flow and prevent cascades
    updating_from_checkbox <- reactiveVal(FALSE)
    updating_select_all <- reactiveVal(FALSE)
    is_loading_state <- reactiveVal(FALSE)
    ignore_input_changes <- reactiveVal(FALSE)
    
    # Use dedicated reactiveVals for stable storage
    stable_condition_column <- reactiveVal(NULL)
    stable_labels <- reactiveVal(list())
    stable_active <- reactiveVal(list())
    
    # Initialize state
    state <- reactiveValues(
      available_columns = NULL,
      temp_labels = NULL,
      last_update = NULL
    )
    
    # Get available metadata columns
    observe({
      # Skip if we're in loading state
      if (is_loading_state()) {
        return(NULL)
      }
      
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
      if (is.null(stable_condition_column()) && length(filtered_cols) > 0) {
        # Prioritize condition-like columns
        priority_cols <- grep("condition|treatment|group|genotype|timepoint", 
                              filtered_cols, value = TRUE, ignore.case = TRUE)
        
        if (length(priority_cols) > 0) {
          stable_condition_column(priority_cols[1])
        } else if (length(filtered_cols) > 0) {
          stable_condition_column(filtered_cols[1])
        }
        
        # Initialize conditions for the selected column
        updateConditionState(stable_condition_column())
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
      
      # Get the current selected column
      current_column <- stable_condition_column()
      
      selectInput(ns("conditionColumn"), 
                  "Select condition column:", 
                  choices = available_cols,
                  selected = current_column)
    })
    
    # Update condition state when dropdown selection changes
    observeEvent(input$conditionColumn, {
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      stable_condition_column(input$conditionColumn)
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
        
        # Initialize labels
        current_labels <- stable_labels()
        new_labels <- initializeConditionLabels(condition_values, current_labels)
        stable_labels(new_labels)
        state$temp_labels <- new_labels
        
        # Initialize active status
        current_active <- stable_active()
        new_active <- initializeActiveStatus(condition_values, current_active)
        stable_active(new_active)
      }
    }
    
    #' @title Get Available Conditions
    #' @description Get unique condition values from the selected column
    #' @return Vector of unique condition values or NULL if none available
    #' @keywords internal
    getAvailableConditions <- function() {
      column <- stable_condition_column()
      if (is.null(column) || is.null(seurat_data())) return(NULL)
      
      if (column %in% colnames(seurat_data()@meta.data)) {
        unique(as.character(seurat_data()@meta.data[[column]]))
      } else {
        NULL
      }
    }
    
    # Setup UI for condition rows
    output$conditionRows <- renderUI({
      column <- stable_condition_column()
      req(column)
      
      available_conditions <- getAvailableConditions()
      
      # Show message if no conditions
      if (is.null(available_conditions) || length(available_conditions) == 0) {
        return(div(
          class = "alert alert-info",
          paste0("No condition values found in column '", column, "'. ",
                 "Please select a different column.")
        ))
      }
      
      # Get current values from stable storage
      current_labels <- stable_labels()
      current_active <- stable_active()
      
      # Calculate if all conditions are currently selected
      all_selected <- !is.null(current_active) && length(current_active) > 0 && all(unlist(current_active))
      
      # Build UI
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
          createConditionControls(ns, available_conditions, current_labels, current_active)
        )
      )
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
    
    # Handle selectAllConditions checkbox
    observeEvent(input$selectAllConditions, {
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
      
      req(getAvailableConditions())
      available_conditions <- getAvailableConditions()
      
      # Create a new active_conditions list
      new_active <- setNames(
        rep(input$selectAllConditions, length(available_conditions)),
        available_conditions
      )
      
      # Block individual checkbox observers
      ignore_input_changes(TRUE)
      
      # Update stable storage with our copy
      stable_active(new_active)
      
      # Update individual checkboxes
      for (condition in available_conditions) {
        safe_id <- make_safe_id(condition)
        updateCheckboxInput(session, paste0("active_", safe_id), value = input$selectAllConditions)
      }
      
      # Now it's safe to allow changes - slight delay to ensure UI updates first
      shinyjs::delay(200, {
        ignore_input_changes(FALSE)
        updating_select_all(FALSE)
        
        # Trigger state update
        state$last_update <- Sys.time()
      })
    }, ignoreInit = TRUE, priority = 10)  # Higher priority
    
    # Individual condition checkbox handling
    observe({
      req(getAvailableConditions())
      available_conditions <- getAvailableConditions()
      
      # Setup observers for each condition
      lapply(available_conditions, function(condition) {
        safe_id <- make_safe_id(condition)
        input_id <- paste0("active_", safe_id)
        
        # Create observer for this specific input
        observeEvent(input[[input_id]], {
          # Skip if in loading state or during Select All update
          if (ignore_input_changes() || updating_select_all()) {
            return(NULL)
          }
          
          # Mark that we're updating from a checkbox
          updating_from_checkbox(TRUE)
          
          # Get current active status - create a complete copy to avoid reference issues
          current_active <- deep_copy(stable_active())
          
          # Get previous value with default
          old_value <- if (!is.null(current_active) && condition %in% names(current_active)) {
            current_active[[condition]]
          } else {
            TRUE
          }
          
          # Check if value actually changed
          if (old_value != input[[input_id]]) {
            # Update our deep copy with the new value
            current_active[[condition]] <- input[[input_id]]
            
            # Update stable active status
            stable_active(current_active)
            
            # Signal update
            state$last_update <- Sys.time()
          }
          
          # Clear updating flag after a delay
          shinyjs::delay(100, {
            updating_from_checkbox(FALSE)
          })
        }, ignoreInit = TRUE, priority = 20)  # Very high priority
      })
    })
    
    # Separate observer just for syncing the "Select All" checkbox
    observe({
      # Skip during active updates
      if (updating_from_checkbox() || updating_select_all() || ignore_input_changes()) {
        return(NULL)
      }
      
      # Get current active status
      current_active <- stable_active()
      
      # Skip if not ready yet
      if (is.null(current_active) || length(current_active) == 0) {
        return(NULL)
      }
      
      # Calculate whether all are selected
      all_selected <- all(unlist(current_active))
      
      # Check current value of the checkbox
      current_state <- input$selectAllConditions
      
      # Only update if different and not null
      if (!is.null(current_state) && all_selected != current_state) {
        # Temporarily block updates
        ignore_input_changes(TRUE)
        updating_select_all(TRUE)
        
        # Update the checkbox - this will NOT trigger its observer due to our flags
        updateCheckboxInput(session, "selectAllConditions", value = all_selected)
        
        # Release locks after a delay
        shinyjs::delay(200, {
          ignore_input_changes(FALSE)
          updating_select_all(FALSE)
        })
      }
    })
    
    # Handle label inputs with dedicated observers
    observe({
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      req(getAvailableConditions())
      available_conditions <- getAvailableConditions()
      
      # Setup observers for each condition's label
      lapply(available_conditions, function(condition) {
        safe_id <- make_safe_id(condition)
        input_id <- paste0("label_", safe_id)
        
        # Create observer for this specific input
        observeEvent(input[[input_id]], {
          # Skip if in loading state
          if (ignore_input_changes()) {
            return(NULL)
          }
          
          # Update temp labels
          state$temp_labels[[condition]] <- input[[input_id]]
        }, ignoreInit = TRUE)
      })
    })
    
    # Update permanent labels when update button is clicked
    observeEvent(input$updateAllLabels, {
      # Skip if in loading state
      if (ignore_input_changes()) {
        return(NULL)
      }
      
      # Update the stable labels from the temporary ones
      stable_labels(deep_copy(state$temp_labels))
      
      # Trigger update for reactivity
      state$last_update <- Sys.time()
      
      showNotification("Condition labels updated", type = "message")
    })
    
    # Function to get state for saving
    getFullState = function() {
      list(
        condition_column = stable_condition_column(),
        condition_labels = stable_labels(),
        active_conditions = stable_active(),
        temp_labels = state$temp_labels,
        available_columns = state$available_columns
      )
    }
    
    # Function to restore state
    setFullState = function(saved_state) {
      if (!is.null(saved_state)) {
        # Set loading flags
        is_loading_state(TRUE)
        ignore_input_changes(TRUE)
        
        # Block updates during loading
        shinyjs::delay(100, {
          # Restore available columns if provided
          if (!is.null(saved_state$available_columns)) {
            state$available_columns <- saved_state$available_columns
          }
          
          # Restore condition column
          if (!is.null(saved_state$condition_column)) {
            stable_condition_column(saved_state$condition_column)
            
            # Update UI dropdown
            updateSelectInput(session, "conditionColumn", selected = saved_state$condition_column)
          }
          
          # Restore labels
          if (!is.null(saved_state$condition_labels)) {
            stable_labels(deep_copy(saved_state$condition_labels))
            
            # Also restore temp labels if available
            if (!is.null(saved_state$temp_labels)) {
              state$temp_labels <- deep_copy(saved_state$temp_labels)
            } else {
              state$temp_labels <- deep_copy(saved_state$condition_labels)
            }
          }
          
          # Restore active status
          if (!is.null(saved_state$active_conditions)) {
            stable_active(deep_copy(saved_state$active_conditions))
          }
          
          # Update UI after a short delay
          shinyjs::delay(300, {
            # Get available conditions after column is set
            available_conditions <- getAvailableConditions()
            
            if (!is.null(available_conditions) && length(available_conditions) > 0) {
              updating_select_all(TRUE)
              
              # Update checkboxes
              if (!is.null(saved_state$active_conditions)) {
                for (condition in names(saved_state$active_conditions)) {
                  if (condition %in% available_conditions) {
                    safe_id <- make_safe_id(condition)
                    updateCheckboxInput(session, paste0("active_", safe_id), 
                                        value = saved_state$active_conditions[[condition]])
                  }
                }
                
                # Update "Select All" checkbox
                all_active <- all(unlist(saved_state$active_conditions))
                updateCheckboxInput(session, "selectAllConditions", value = all_active)
              }
              
              # Update label inputs
              if (!is.null(saved_state$condition_labels)) {
                for (condition in names(saved_state$condition_labels)) {
                  if (condition %in% available_conditions) {
                    safe_id <- make_safe_id(condition)
                    updateTextInput(session, paste0("label_", safe_id), 
                                    value = saved_state$condition_labels[[condition]])
                  }
                }
              }
              
              updating_select_all(FALSE)
            }
            
            # Force UI refresh
            state$last_update <- Sys.time()
            
            # Release locks after delay
            shinyjs::delay(500, {
              ignore_input_changes(FALSE)
              shinyjs::delay(200, {
                is_loading_state(FALSE)
              })
            })
          })
        })
      }
    }
    
    # Return reactive expressions
    list(
      getConditionColumn = reactive({ stable_condition_column() }),
      getConditionLabels = reactive({ stable_labels() }),
      getActiveStatus = reactive({ stable_active() }),
      getActiveConditions = reactive({
        active_conditions <- stable_active()
        if (is.null(active_conditions)) return(NULL)
        names(active_conditions[unlist(active_conditions) == TRUE])
      }),
      getActiveConditionList = reactive({
        active_conditions <- stable_active()
        if (is.null(active_conditions)) return(NULL)
        names(active_conditions[unlist(active_conditions) == TRUE])
      }),
      updateLabels = function(new_labels) {
        stable_labels(deep_copy(new_labels))
        state$temp_labels <- deep_copy(new_labels)
        
        # Trigger update
        state$last_update <- Sys.time()
      },
      updateActiveStatus = function(new_active) {
        # Skip if we're already updating internally
        if (updating_from_checkbox() || updating_select_all()) {
          return()
        }
        
        # Update with a deep copy
        stable_active(deep_copy(new_active))
        
        # Trigger update
        state$last_update <- Sys.time()
      },
      getFullState = getFullState,
      setFullState = setFullState,
      updateFromButton = function() {
        # Update the stable labels from the temporary ones
        if (!is.null(state$temp_labels) && length(state$temp_labels) > 0) {
          stable_labels(deep_copy(state$temp_labels))
          
          # Trigger update for reactivity
          state$last_update <- Sys.time()
          
          showNotification("Condition labels updated", type = "message")
        } else {
          showNotification("No label changes to save", type = "message")
        }
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