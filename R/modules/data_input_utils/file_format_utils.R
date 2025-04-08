#' @title Read H5 File
#' @description Reads a 10X Genomics HDF5 file and creates a sparse matrix
#' @param file_path Path to the H5 file
#' @return A list containing a sparse matrix of counts
#' @export
Read_10X_H5 <- function(file_path) {
  # Check if hdf5r is available
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required to read H5 files. Please install it.")
  }
  
  # Use Read10X_h5 from Seurat to read 10X h5 file
  message("Reading H5 file:", file_path)
  counts <- Seurat::Read10X_h5(file_path)
  
  # If it's a list (multimodal data), use just the Gene Expression
  if (is.list(counts) && "Gene Expression" %in% names(counts)) {
    counts <- counts[["Gene Expression"]]
  }
  
  return(list(counts))
}

#' @title Detect File Format
#' @description Detects the format of a file based on its extension
#' @param file_name Name of the file
#' @return A string indicating the file format ("txt.gz", "h5", or "unknown")
#' @export
Detect_File_Format <- function(file_name) {
  if (grepl("\\.txt\\.gz$", file_name)) {
    return("txt.gz")
  } else if (grepl("\\.h5$", file_name)) {
    return("h5")
  } else {
    return("unknown")
  }
}

#' @title Read Data File
#' @description Reads a data file based on its format
#' @param data_dir Directory containing the data file
#' @param file_name Name of the file
#' @return A list containing a sparse matrix of counts
#' @export
Read_Data_File <- function(data_dir, file_name) {
  file_path <- file.path(data_dir, file_name)
  file_format <- Detect_File_Format(file_name)
  
  if (file_format == "txt.gz") {
    return(Read_GEO_Delim(data_dir, file_name))
  } else if (file_format == "h5") {
    return(Read_10X_H5(file_path))
  } else {
    stop(paste("Unsupported file format for file:", file_name))
  }
}

#' @title Calculate Mitochondrial Percentage for Seurat Object
#' @description Identifies mitochondrial genes and calculates percentage for a Seurat object
#' @param seurat_obj Seurat object
#' @return The Seurat object with percent.mt column added
#' @export
Calculate_MT_Percent <- function(seurat_obj) {
  # Get gene names
  gene_names <- rownames(seurat_obj)
  
  # Check for presence of common MT gene patterns
  has_human_mt <- any(grepl("^MT-", gene_names))
  has_mouse_mt <- any(grepl("^mt-", gene_names, ignore.case = TRUE))
  has_ensembl_mt <- any(grepl("^ENSMUSG00000064", gene_names))
  has_generic_mt <- any(grepl("^[Mm][Tt]-|^[Mm][Tt]:|^[Mm]ito", gene_names))
  
  print(paste("Has human MT genes:", has_human_mt))
  print(paste("Has mouse MT genes:", has_mouse_mt))
  print(paste("Has Ensembl MT genes:", has_ensembl_mt))
  print(paste("Has generic MT pattern:", has_generic_mt))
  
  # Try different patterns
  if (has_human_mt) {
    print("Using human MT- pattern")
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  } else if (has_mouse_mt) {
    # Instead of using ignore.case parameter, use a case-insensitive regex with both options
    print("Using mouse mt- pattern")
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^(mt-|MT-)")
  } else if (has_ensembl_mt) {
    print("Using Ensembl mitochondrial pattern")
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^ENSMUSG00000064")
  } else if (has_generic_mt) {
    print("Using generic mitochondrial pattern")
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^[Mm][Tt]-|^[Mm][Tt]:|^[Mm]ito")
  } else {
    # Look for any MT-related feature
    mt_patterns <- c("MT", "Mt", "mt", "Mito", "mito")
    found_pattern <- FALSE
    
    for (pattern in mt_patterns) {
      # For case-insensitive search, we can use grepl directly
      mt_indices <- grep(pattern, gene_names, ignore.case = TRUE)
      if (length(mt_indices) > 0) {
        matching_features <- gene_names[mt_indices]
        print(paste("Found", length(matching_features), "mitochondrial genes using pattern:", pattern))
        print(paste("Example matches:", paste(head(matching_features, 3), collapse=", ")))
        seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = matching_features)
        found_pattern <- TRUE
        break
      }
    }
    
    if (!found_pattern) {
      print("No mitochondrial genes identified. Setting percent.mt to 0.")
      seurat_obj[["percent.mt"]] <- 0
    }
  }
  
  # Check the range of percent.mt values
  print(paste("Range of percent.mt:", min(seurat_obj$percent.mt), "to", max(seurat_obj$percent.mt)))
  
  return(seurat_obj)
}