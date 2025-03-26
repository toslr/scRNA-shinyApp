# R/modules/condition_management_module.R

conditionManagementUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("conditionSelector")),
    uiOutput(ns("conditionRows"))
  )
}

conditionManagementServer <- function(id, seurat_data, metadata_module) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize state
    state <- list(
      condition_column = reactiveVal(NULL),
      condition_labels = reactiveVal(NULL),
      active_conditions = reactiveVal(NULL),
      temp_labels = reactiveVal(NULL),
      label_inputs = reactiveValues(),
      available_columns = reactiveVal(NULL)
    )
    
    # Get available metadata columns
    observe({
      # Skip if no Seurat object
      if (is.null(seurat_data())) return(NULL)
      
      # Get metadata columns
      meta_cols <- colnames(seurat_data()@meta.data)
      
      # Filter out standard columns we don't want to show
      exclude_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA", 
                           "seurat_clusters", "percent.mt")
      filtered_cols <- setdiff(meta_cols, exclude_columns)
      
      # Set available columns
      state$available_columns(filtered_cols)
    })
    
    # Render condition column selector
    output$conditionSelector <- renderUI({
      available_cols <- state$available_columns()
      
      if (is.null(available_cols) || length(available_cols) == 0) {
        return(div(
          class = "alert alert-info",
          "No metadata columns available. Please load data first."
        ))
      }
      
      # Find the current selected column or set default
      current_column <- state$condition_column()
      
      # Prioritize columns that look like conditions
      priority_cols <- grep("condition|treatment|group|genotype|timepoint", 
                            available_cols, value = TRUE, ignore.case = TRUE)
      
      selected_col <- if (!is.null(current_column) && current_column %in% available_cols) {
        current_column
      } else if (length(priority_cols) > 0) {
        priority_cols[1]
      } else if ("sample" %in% available_cols) {
        "sample"
      } else {
        available_cols[1]
      }
      
      div(
        selectInput(ns("conditionColumn"), 
                    "Select condition column:", 
                    choices = available_cols,
                    selected = selected_col),
        actionButton(ns("updateConditionColumn"), 
                     "Apply", 
                     class = "btn-sm btn-primary")
      )
    })
    
    # Update condition column when Apply button is clicked
    observeEvent(input$updateConditionColumn, {
      req(input$conditionColumn)
      # Update the column selection
      state$condition_column(input$conditionColumn)
      
      # Reset condition management state for the new column
      updateConditionState()
    })
    
    # Initialize or update condition management state
    updateConditionState <- function() {
      column <- state$condition_column()
      if (is.null(column) || is.null(seurat_data())) return(NULL)
      
      # Get unique condition values
      if (column %in% colnames(seurat_data()@meta.data)) {
        condition_values <- unique(as.character(seurat_data()@meta.data[[column]]))
        
        # Skip if no values
        if (length(condition_values) == 0) return(NULL)
        
        # Initialize labels and active status
        state$condition_labels(initializeConditionLabels(condition_values))
        state$temp_labels(state$condition_labels())
        state$active_conditions(initializeActiveStatus(condition_values))
      }
    }
    
    # Get available conditions from selected column
    getAvailableConditions <- function() {
      column <- state$condition_column()
      if (is.null(column) || is.null(seurat_data())) return(NULL)
      
      if (column %in% colnames(seurat_data()@meta.data)) {
        unique(as.character(seurat_data()@meta.data[[column]]))
      } else {
        NULL
      }
    }
    
    # Handle initial population when Seurat data changes
    observeEvent(seurat_data(), {
      # Reset available columns
      meta_cols <- if (!is.null(seurat_data())) {
        colnames(seurat_data()@meta.data)
      } else {
        NULL
      }
      
      if (!is.null(meta_cols)) {
        # Filter out standard columns
        exclude_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA", 
                             "seurat_clusters", "percent.mt")
        filtered_cols <- setdiff(meta_cols, exclude_columns)
        
        # Set available columns
        state$available_columns(filtered_cols)
        
        # Try to set a default condition column
        priority_cols <- grep("condition|treatment|group|genotype|timepoint", 
                              filtered_cols, value = TRUE, ignore.case = TRUE)
        
        if (length(priority_cols) > 0) {
          state$condition_column(priority_cols[1])
        } else if ("sample" %in% filtered_cols) {
          state$condition_column("sample")
        } else if (length(filtered_cols) > 0) {
          state$condition_column(filtered_cols[1])
        }
        
        # Initialize condition state
        updateConditionState()
      }
    }, ignoreInit = TRUE)
    
    # Setup UI for condition rows
    output$conditionRows <- renderUI({
      column <- state$condition_column()
      if (is.null(column)) {
        return(div(
          class = "alert alert-info",
          "Please select a condition column above and click Apply."
        ))
      }
      
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
      current_temp_labels <- state$temp_labels()
      current_active <- state$active_conditions()
      
      # Return UI
      tagList(
        div(
          style = "display: flex; align-items: center; margin: 15px 0 10px 0;",
          checkboxInput(ns("selectAllConditions"), "Select All Conditions", value = TRUE),
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
    
    # Handle Select All checkbox
    observeEvent(input$selectAllConditions, {
      available_conditions <- getAvailableConditions()
      if (is.null(available_conditions) || length(available_conditions) == 0) {
        return(NULL)
      }
      
      current_active <- state$active_conditions()
      for (condition in available_conditions) {
        current_active[condition] <- input$selectAllConditions
      }
      state$active_conditions(current_active)
    })
    
    # Handle active status updates for individual conditions
    observe({
      available_conditions <- getAvailableConditions()
      
      # Skip empty conditions
      if (is.null(available_conditions) || length(available_conditions) == 0) {
        return(NULL)
      }
      
      lapply(available_conditions, function(condition) {
        input_id <- paste0("active_", make_safe_id(condition))
        
        # Only observe if this input exists
        if (!is.null(input[[input_id]])) {
          observeEvent(input[[input_id]], {
            current_active <- state$active_conditions()
            if (is.null(current_active)) {
              return(NULL)
            }
            
            # Check if this condition exists in our active_conditions
            if (condition %in% names(current_active)) {
              current_active[condition] <- input[[input_id]]
              state$active_conditions(current_active)
              
              # Update select all checkbox based on all condition checkboxes
              all_selected <- all(unlist(current_active))
              updateCheckboxInput(session, "selectAllConditions", value = all_selected)
            }
          }, ignoreInit = TRUE)
        }
      })
    })
    
    # Update label_inputs when text changes
    observe({
      available_conditions <- getAvailableConditions()
      if (is.null(available_conditions) || length(available_conditions) == 0) {
        return(NULL)
      }
      
      for(condition in available_conditions) {
        local({
          local_condition <- condition
          input_id <- paste0("label_", make_safe_id(local_condition))
          
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
      available_conditions <- getAvailableConditions()
      if (is.null(available_conditions) || length(available_conditions) == 0) {
        return(NULL)
      }
      
      current_temp <- state$temp_labels()
      
      for(condition in available_conditions) {
        input_id <- paste0("label_", make_safe_id(condition))
        if(!is.null(state$label_inputs[[input_id]])) {
          current_temp[condition] <- state$label_inputs[[input_id]]
        }
      }
      
      state$temp_labels(current_temp)
      state$condition_labels(current_temp)
      
      showNotification("Condition labels updated", type = "message")
    })
    
    # Return reactive expressions
    list(
      getConditionColumn = reactive({ state$condition_column() }),
      getConditionLabels = reactive({ state$condition_labels() }),
      getActiveStatus = reactive({ state$active_conditions() }),
      getActiveConditions = reactive({
        current_active <- state$active_conditions()
        if (is.null(current_active)) return(NULL)
        names(current_active[current_active == TRUE])
      }),
      getActiveConditionList = reactive({
        current_active <- state$active_conditions()
        if (is.null(current_active)) return(NULL)
        names(current_active[current_active == TRUE])
      }),
      updateLabels = function(new_labels) {
        state$condition_labels(new_labels)
        state$temp_labels(new_labels)
      },
      updateActiveStatus = function(new_active) {
        state$active_conditions(new_active)
      }
    )
  })
}

# Helper function: Create UI controls for each condition
createConditionControls <- function(ns, available_conditions, current_temp_labels, current_active) {
  tagList(
    # Individual condition controls
    lapply(available_conditions, function(condition) {
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

# Initialize condition labels
initializeConditionLabels <- function(conditions, current_labels = NULL) {
  new_labels <- setNames(
    conditions,  # Use condition values as default labels
    conditions
  )
  
  # Merge with existing labels if they exist
  if (!is.null(current_labels)) {
    existing_conditions <- names(current_labels)
    new_labels[existing_conditions] <- current_labels[existing_conditions]
  }
  
  return(new_labels)
}

# Helper function to make safe HTML/Shiny IDs from condition values
make_safe_id <- function(id) {
  # Replace special characters with underscores
  safe_id <- gsub("[^a-zA-Z0-9]", "_", id)
  
  # Ensure it starts with a letter (Shiny input IDs must start with a letter)
  if (!grepl("^[a-zA-Z]", safe_id)) {
    safe_id <- paste0("c_", safe_id)
  }
  
  return(safe_id)
}

# Initialize active status for conditions
initializeActiveStatus <- function(conditions, current_active = NULL) {
  new_active <- setNames(
    rep(TRUE, length(conditions)),
    conditions
  )
  
  if (!is.null(current_active)) {
    existing_conditions <- names(current_active)
    new_active[existing_conditions] <- current_active[existing_conditions]
  }
  
  return(new_active)
}

# Filter Seurat object by active conditions
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