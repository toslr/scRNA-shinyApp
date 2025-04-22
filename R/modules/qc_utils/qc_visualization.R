#' @title Create QC Plot
#' @description Creates violin plots for key QC metrics with improved formatting for large datasets
#' @param seurat_obj Seurat object to visualize
#' @return A ggplot object with violin plots for QC metrics
#' @keywords internal
createQCPlot <- function(seurat_obj) {
  # Check if we have sample_label column to use instead
  group_var <- if ("sample_label" %in% colnames(seurat_obj@meta.data)) {
    "sample_label"
  } else {
    "sample"
  }
  
  # Create plots for metrics we know exist
  plot_list <- list()
  
  # nFeature plot (should always exist)
  if ("nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA", group.by = group_var, pt.size = 0.1) +
      theme(legend.position = "none") +
      ggtitle("Number of Features")
    plot_list[["nFeature"]] <- p1
  }
  
  # nCount plot (should always exist)
  if ("nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    p2 <- VlnPlot(seurat_obj, features = "nCount_RNA", group.by = group_var, pt.size = 0.1) +
      theme(legend.position = "none") +
      ggtitle("Number of UMIs")
    plot_list[["nCount"]] <- p2
  }
  
  # percent.mt plot - check if it exists and has valid values
  if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
    # Handle NA and infinite values
    if (any(is.na(seurat_obj$percent.mt))) {
      seurat_obj$percent.mt[is.na(seurat_obj$percent.mt)] <- 0
    }
    
    if (any(is.infinite(seurat_obj$percent.mt))) {
      seurat_obj$percent.mt[is.infinite(seurat_obj$percent.mt)] <- 0
    }
    
    # Create the plot
    p3 <- VlnPlot(seurat_obj, features = "percent.mt", group.by = group_var, pt.size = 0.1) +
      theme(legend.position = "none") +
      ggtitle("Mitochondrial %")
    plot_list[["percent.mt"]] <- p3
  }
  
  # Combine the plots
  if (length(plot_list) == 0) {
    # Return an empty plot with message if no plots created
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No QC metrics available to plot") + 
             theme_void())
  } else if (length(plot_list) == 1) {
    # Return the single plot
    return(plot_list[[1]])
  } else {
    # Combine multiple plots
    p_combined <- patchwork::wrap_plots(plot_list, ncol = length(plot_list))
    return(p_combined)
  }
}

#' @title Save QC Plot
#' @description Saves the QC plot to a file with error handling.
#' @param file File path to save to
#' @param plot The plot to save
#' @return None
#' @keywords internal
saveQCPlot <- function(file, plot) {
  tryCatch({
    ggsave(file, plot = plot, device = "png", width = 12, height = 8, dpi = 300)
  }, error = function(e) {
    # Fallback to base R graphics if ggsave fails
    png(file, width = 12, height = 8, units = "in", res = 300)
    print(plot)
    dev.off()
  })
}

#' @title Create Filter Controls
#' @description Creates UI elements for filtering the Seurat object based on QC metrics.
#' @param session The current Shiny session
#' @param values Reactive values containing filter parameters
#' @return A UI element with filter controls
#' @keywords internal
createFilterControls <- function(session, values) {
  ns <- session$ns
  tagList(
    tags$div(
      id = ns("filter_controls"),
      tags$br(),
      tags$strong("Please adjust filtering parameters:"),
      tags$br(),
      tags$br(),
      div(id = ns("minFeature_container"),
          numericInput(ns("minFeature"), "Minimum Features:", values$min_feature)
      ),
      div(id = ns("maxFeature_container"),
          numericInput(ns("maxFeature"), "Maximum Features:", values$max_feature)
      ),
      div(id = ns("maxMT_container"),
          numericInput(ns("maxMT"), "Maximum MT %:", values$max_mt)
      ),
      actionButton(ns("processSeurat"), "Filter and run PCA")
    )
  )
}

#' @title Get Samples to Plot
#' @description Selects a subset of samples to display in QC plots (limits to first 5 if more exist).
#' @param seurat_obj Seurat object containing the data
#' @return A Seurat object with only the selected samples
#' @keywords internal
getSamplesToPlot <- function(seurat_obj) {
  all_samples <- unique(seurat_obj$sample)
  
  if (length(all_samples) > 5) {
    samples_to_plot <- all_samples[1:5]
    
    # Check if we'll have cells left
    cells_to_keep <- seurat_obj$sample %in% samples_to_plot
    cells_count <- sum(cells_to_keep)
    
    if (cells_count == 0) {
      return(seurat_obj)
    }
    
    # Try to subset safely
    tryCatch({
      return(subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep]))
    }, error = function(e) {
      return(seurat_obj)
    })
  }
  
  return(seurat_obj)
}

#' @title Get QC Suggestions
#' @description Get species-specific QC parameter suggestions
#' @param species Species identification string
#' @return A list of recommended QC parameters
#' @keywords internal
get_qc_suggestions <- function(species) {
  # Default suggestions
  default_min_feature <- 500
  default_max_feature <- 5000
  default_max_mt <- 5
  
  # Species-specific suggestions
  if (!is.null(species)) {
    if (species == "human") {
      default_min_feature <- 500
      default_max_feature <- 6000
      default_max_mt <- 10
    } else if (species == "mouse") {
      default_min_feature <- 500
      default_max_feature <- 5000
      default_max_mt <- 5
    } else if (species == "zebrafish") {
      default_min_feature <- 300
      default_max_feature <- 3000
      default_max_mt <- 8
    } else if (species == "fly") {
      default_min_feature <- 300
      default_max_feature <- 3000
      default_max_mt <- 10
    }
  }
  
  return(list(
    min_feature = default_min_feature,
    max_feature = default_max_feature,
    max_mt = default_max_mt
  ))
}

#' @title Get Species from Seurat
#' @description Extract species information from a Seurat object
#' @param seurat_obj A Seurat object
#' @return A list containing the species information
#' @keywords internal
get_species_from_seurat <- function(seurat_obj) {
  # Try to get species information
  species <- NULL
  tryCatch({
    if ("species" %in% colnames(seurat_obj@meta.data)) {
      species <- seurat_obj$species[1]
    }
  }, error = function(e) {
    # Species info not available
  })
  
  return(list(
    species = species
  ))
}

#' @title Build QC Status Message
#' @description Creates a status message with information about filtering
#' @param seurat_obj Seurat object being processed
#' @param filtered_samples Vector of sample IDs to include
#' @param filtered_conditions Vector of condition values to include
#' @param condition_column Name of the condition column
#' @param min_feature Minimum feature count threshold
#' @param max_feature Maximum feature count threshold
#' @param max_mt Maximum mitochondrial percentage threshold
#' @return A list with message components
#' @keywords internal
build_qc_status_message <- function(seurat_obj, filtered_samples, filtered_conditions, 
                                    condition_column, min_feature, max_feature, max_mt) {
  num_cells <- ncol(seurat_obj)
  
  # Get sample information
  sample_info <- ""
  if (!is.null(filtered_samples) && length(filtered_samples) > 0) {
    sample_info <- paste0(" from ", length(filtered_samples), " samples")
  }
  
  # Get condition information
  condition_info <- ""
  if (!is.null(condition_column) && !is.null(filtered_conditions) && 
      length(filtered_conditions) > 0 && 
      condition_column %in% colnames(seurat_obj@meta.data)) {
    condition_info <- paste0(" in ", length(filtered_conditions), " conditions")
  }
  
  # Filter summary
  filter_summary <- paste0("QC thresholds - Min features: ", min_feature, 
                           ", Max features: ", max_feature, 
                           ", Max MT%: ", max_mt)
  
  # Initial message
  initial_message <- paste0(
    "Processing QC with ", num_cells, " total cells", 
    sample_info, condition_info, "; ", filter_summary
  )
  
  return(list(
    initial_message = initial_message,
    sample_info = sample_info,
    condition_info = condition_info,
    filter_summary = filter_summary
  ))
}