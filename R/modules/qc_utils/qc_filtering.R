#' @title Filter Seurat for QC Visualization
#' @description Applies filtering to a Seurat object for QC visualization
#' @param seurat_obj Seurat object to filter
#' @param filtered_samples Vector of sample IDs to include
#' @param filtered_conditions Vector of condition values to include
#' @param condition_column Name of condition column
#' @param sample_labels Named list of sample labels
#' @return Filtered Seurat object
#' @keywords internal
filter_seurat_for_qc <- function(seurat_obj, filtered_samples, filtered_conditions, 
                                 condition_column, sample_labels = NULL) {
  # Create a filtered Seurat object
  filtered_seurat <- seurat_obj
  
  # Apply sample filtering if available
  if (!is.null(filtered_samples) && length(filtered_samples) > 0) {
    # Filter to show only cells from active samples
    cells_to_keep <- filtered_seurat$sample %in% filtered_samples
    if (any(cells_to_keep)) {
      filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
    }
  }
  
  # Apply condition filtering if available
  if (!is.null(condition_column) && !is.null(filtered_conditions) && 
      length(filtered_conditions) > 0 && 
      condition_column %in% colnames(filtered_seurat@meta.data)) {
    
    # Filter to show only cells from active conditions
    cells_to_keep <- filtered_seurat@meta.data[[condition_column]] %in% filtered_conditions
    if (any(cells_to_keep)) {
      filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
    }
  }
  
  # Apply sample limit to prevent overcrowding (get first 5 samples)
  filtered_seurat <- getSamplesToPlot(filtered_seurat)
  
  # For large datasets, downsample to improve plotting performance
  if (ncol(filtered_seurat) > 10000) {
    set.seed(42)  # For reproducibility
    cells_to_keep <- sample(colnames(filtered_seurat), min(10000, ncol(filtered_seurat)))
    filtered_seurat <- subset(filtered_seurat, cells = cells_to_keep)
  }
  
  # If we have sample labels, apply them to the plot data
  if (!is.null(sample_labels)) {
    # Create temporary column for sample labels if needed
    if (!"sample_label" %in% colnames(filtered_seurat@meta.data)) {
      filtered_seurat$sample_label <- filtered_seurat$sample
    }
    
    # For each sample in the dataset, update its display name if we have a label
    current_samples <- unique(filtered_seurat$sample)
    for (sample in current_samples) {
      if (sample %in% names(sample_labels)) {
        # Replace sample IDs with their labels
        filtered_seurat$sample_label[filtered_seurat$sample == sample] <- sample_labels[[sample]]
      }
    }
  }
  
  return(filtered_seurat)
}

#' @title Filter Seurat for Processing
#' @description Applies filtering to a Seurat object for QC processing
#' @param seurat_obj Seurat object to filter
#' @param filtered_samples Vector of sample IDs to include
#' @param filtered_conditions Vector of condition values to include
#' @param condition_column Name of condition column
#' @return Filtered Seurat object
#' @keywords internal
filter_seurat_for_processing <- function(seurat_obj, filtered_samples, filtered_conditions, condition_column) {
  filtered_seurat <- seurat_obj
  
  # Apply sample filtering if available
  if (!is.null(filtered_samples) && length(filtered_samples) > 0) {
    # Filter to show only cells from active samples
    cells_to_keep <- filtered_seurat$sample %in% filtered_samples
    if (any(cells_to_keep)) {
      filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
    } else {
      stop("No cells match the active sample selection. Cannot process data.")
    }
  }
  
  # Apply condition filtering if available
  if (!is.null(condition_column) && !is.null(filtered_conditions) && 
      length(filtered_conditions) > 0 && 
      condition_column %in% colnames(filtered_seurat@meta.data)) {
    
    # Filter to show only cells from active conditions
    cells_to_keep <- filtered_seurat@meta.data[[condition_column]] %in% filtered_conditions
    if (any(cells_to_keep)) {
      filtered_seurat <- subset(filtered_seurat, cells = colnames(filtered_seurat)[cells_to_keep])
    } else {
      stop("No cells match the active condition selection. Cannot process data.")
    }
  }
  
  return(filtered_seurat)
}

#' @title Process QC Filtering
#' @description Filters cells based on QC metrics and runs preprocessing steps (normalization, scaling, PCA).
#' @param seurat_obj The Seurat object to process
#' @param min_feature Minimum feature count threshold
#' @param max_feature Maximum feature count threshold
#' @param max_mt Maximum mitochondrial percentage threshold
#' @return A processed Seurat object
#' @export
processQCFiltering <- function(seurat_obj, min_feature, max_feature, max_mt) {
  withProgress(message = 'Processing data', value=0, {
    # Filter cells
    incProgress(0.05, detail = "Filtering cells")
    seurat <- subset(seurat_obj, subset = nFeature_RNA > min_feature &
                       nFeature_RNA < max_feature &
                       percent.mt < max_mt)
    
    # Normalize data
    incProgress(0.2, detail = "Normalizing data")
    seurat <- JoinLayers(seurat)
    seurat <- NormalizeData(seurat)
    
    # Find variable features
    incProgress(0.2, detail = "Finding variable features")
    seurat <- FindVariableFeatures(seurat)
    
    # Scale data
    incProgress(0.2, detail = "Scaling data")
    seurat <- ScaleData(seurat)
    
    # Run PCA
    incProgress(0.2, detail = "Running PCA")
    seurat <- RunPCA(seurat, npcs = 50)
    
    return(seurat)
  })
}