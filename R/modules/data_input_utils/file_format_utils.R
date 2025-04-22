#' @title Read 10X H5 File
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
  print(paste("Reading H5 file:", file_path))
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
#' @param file_name Name of the file or GSM ID
#' @return A list containing a sparse matrix of counts
#' @export
Read_Data_File <- function(data_dir, file_name) {
  # Determine file format type
  format_type <- Detect_File_Format_Type(data_dir, file_name)
  print(paste("Detected file format:", format_type, "for", file_name))
  
  if (format_type == "h5") {
    # For H5 files
    h5_files <- list.files(data_dir, pattern = paste0("^", file_name, ".*\\.h5$"), full.names = TRUE)
    if (length(h5_files) > 0) {
      return(Read_10X_H5(h5_files[1]))
    } else {
      stop(paste("No H5 files found for", file_name))
    }
  } else if (format_type == "mtx") {
    # For 10X MTX format
    mtx_files <- Detect_10X_MTX_Format(data_dir, file_name)
    if (!is.null(mtx_files)) {
      return(Read_10X_MTX_Format(mtx_files))
    } else {
      stop(paste("MTX format files not found correctly for", file_name))
    }
  } else if (format_type == "txt.gz") {
    # For GEO text format
    return(Read_GEO_Delim(data_dir, file_name))
  } else {
    stop(paste("Unsupported file format for file:", file_name))
  }
}

#' @title Detect 10X MTX Format
#' @description Detects if a set of files represents a 10X Genomics MTX format dataset
#' @param dir_path Directory to check
#' @param gsm_id GSM ID to look for
#' @return A list of file paths if format is found, NULL otherwise
#' @export
Detect_10X_MTX_Format <- function(dir_path, gsm_id) {
  # Check if directory exists
  if (!dir.exists(dir_path)) {
    print(paste("Directory does not exist:", dir_path))
    return(NULL)
  }
  
  # List all files in the directory
  all_files <- list.files(dir_path, full.names = FALSE)
  print(paste("Checking for 10X MTX format files for GSM ID:", gsm_id))
  
  # Look for the triplet of files needed for 10X MTX format
  base_pattern <- paste0("^", gsm_id, ".*")
  
  # Check for the three required files
  barcode_files <- list.files(dir_path, pattern = paste0(base_pattern, "barcodes\\.tsv\\.gz$"), full.names = TRUE)
  features_files <- list.files(dir_path, pattern = paste0(base_pattern, "genes\\.tsv\\.gz$|", base_pattern, "features\\.tsv\\.gz$"), full.names = TRUE)
  matrix_files <- list.files(dir_path, pattern = paste0(base_pattern, "matrix\\.mtx\\.gz$"), full.names = TRUE)
  
  print(paste("Found barcode files:", paste(barcode_files, collapse=", ")))
  print(paste("Found feature/gene files:", paste(features_files, collapse=", ")))
  print(paste("Found matrix files:", paste(matrix_files, collapse=", ")))
  
  # If we found all three files, return their paths
  if (length(barcode_files) > 0 && length(features_files) > 0 && length(matrix_files) > 0) {
    return(list(
      barcodes = barcode_files[1],
      features = features_files[1],
      matrix = matrix_files[1]
    ))
  }
  
  # If we didn't find all files, return NULL
  print("Did not find all three required files for 10X MTX format")
  return(NULL)
}

#' @title Read 10X MTX Format
#' @description Reads 10X Genomics MTX format files into a sparse matrix
#' @param file_paths List of file paths for barcodes, features/genes, and matrix
#' @return A list containing a sparse matrix of counts
#' @export
Read_10X_MTX_Format <- function(file_paths) {
  print(paste("Reading 10X MTX format files from:", 
              paste(unlist(file_paths), collapse=", ")))
  
  # Create a temporary directory
  temp_dir <- file.path(tempdir(), "temp_10x_data")
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)  # Remove old temp dir if it exists
  }
  dir.create(temp_dir, recursive = TRUE)
  print(paste("Created temporary directory:", temp_dir))
  
  # Copy files with standardized names
  file.copy(file_paths$barcodes, file.path(temp_dir, "barcodes.tsv.gz"))
  
  # Decide whether to use "features.tsv.gz" or "genes.tsv.gz" based on original name
  features_file <- basename(file_paths$features)
  if (grepl("features", features_file)) {
    file.copy(file_paths$features, file.path(temp_dir, "features.tsv.gz"))
  } else {
    file.copy(file_paths$features, file.path(temp_dir, "genes.tsv.gz"))
  }
  
  file.copy(file_paths$matrix, file.path(temp_dir, "matrix.mtx.gz"))
  
  # List files in temp dir to confirm
  temp_files <- list.files(temp_dir)
  print(paste("Files in temporary directory:", paste(temp_files, collapse=", ")))
  
  # Try to read the data
  tryCatch({
    counts <- Seurat::Read10X(data.dir = temp_dir)
    print("Successfully read 10X MTX data")
    print(paste("Matrix dimensions:", nrow(counts), "x", ncol(counts)))
    return(list(counts))
  }, error = function(e) {
    print(paste("Error reading 10X MTX data:", e$message))
    print("Attempting to use Seurat::ReadMtx directly...")
    
    # Try using ReadMtx directly
    mtx_file <- file.path(temp_dir, "matrix.mtx.gz")
    features_file <- if ("features.tsv.gz" %in% temp_files) {
      file.path(temp_dir, "features.tsv.gz")
    } else {
      file.path(temp_dir, "genes.tsv.gz")
    }
    barcode_file <- file.path(temp_dir, "barcodes.tsv.gz")
    
    counts <- Seurat::ReadMtx(
      mtx = mtx_file,
      features = features_file,
      cells = barcode_file
    )
    
    print("Successfully read 10X MTX data using ReadMtx")
    print(paste("Matrix dimensions:", nrow(counts), "x", ncol(counts)))
    return(list(counts))
  })
}

#' @title Detect File Format Type
#' @description Determines what type of scRNA-seq file format is being used
#' @param dir_path Directory path
#' @param file_name Filename or GSM ID to examine
#' @return A character string describing the format type: "h5", "mtx", "txt.gz", or "unknown"
#' @export
Detect_File_Format_Type <- function(dir_path, file_name) {
  print(paste("Detecting file format for:", file_name, "in directory:", dir_path))
  
  # Check if it's a direct file that exists
  full_path <- file.path(dir_path, file_name)
  if (file.exists(full_path)) {
    if (grepl("\\.h5$", file_name)) {
      return("h5")
    } else if (grepl("\\.txt\\.gz$", file_name)) {
      return("txt.gz")
    }
  }
  
  # For GSM IDs, we need to check if any matching files exist
  # First check for 10X MTX format - this should be first since it's more specific
  mtx_files <- Detect_10X_MTX_Format(dir_path, file_name)
  if (!is.null(mtx_files)) {
    print("10X MTX format detected")
    return("mtx")
  }
  
  # If not found as MTX, check for H5 files
  print("Checking for H5 files")
  h5_files <- list.files(dir_path, pattern = paste0("^", file_name, ".*\\.h5$"), full.names = FALSE)
  if (length(h5_files) > 0) {
    print(paste("H5 file found:", h5_files[1]))
    return("h5")
  }
  
  # Finally check for txt.gz files
  print("Checking for txt.gz files")
  txt_files <- list.files(dir_path, pattern = paste0("^", file_name, ".*\\.txt\\.gz$"), full.names = FALSE)
  if (length(txt_files) > 0) {
    print(paste("txt.gz file found:", txt_files[1]))
    return("txt.gz")
  }
  
  # If we couldn't determine the format
  print("No supported format detected")
  return("unknown")
}

#' @title Read GEO Delimited File
#' @description Reads a GEO data file in delimited format
#' @param data_dir Directory containing the file
#' @param file_suffix Filename or suffix to read
#' @return A list containing a sparse matrix of gene counts
#' @export
Read_GEO_Delim <- function(data_dir, file_suffix) {
  # This keeps your original implementation
  # Determine the actual file path
  if (file.exists(file.path(data_dir, file_suffix))) {
    # Direct file name provided
    file_path <- file.path(data_dir, file_suffix)
  } else {
    # GSM ID provided, find the matching file
    matching_files <- list.files(data_dir, pattern = paste0("^", file_suffix, ".*\\.txt\\.gz$"), full.names = TRUE)
    if (length(matching_files) == 0) {
      stop(paste("No matching .txt.gz files found for", file_suffix))
    }
    file_path <- matching_files[1]
  }
  
  print(paste("Reading text file:", file_path))
  
  # Read the data
  tryCatch({
    # Try to read with tab as delimiter
    data_matrix <- read.delim(file_path, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
    
    # Convert to sparse matrix for memory efficiency
    sparse_matrix <- Matrix::Matrix(as.matrix(data_matrix), sparse=TRUE)
    
    # Return as a list with a single element (to match other read functions)
    print(paste("Successfully read text file with dimensions:", nrow(sparse_matrix), "x", ncol(sparse_matrix)))
    return(list(sparse_matrix))
  }, error = function(e) {
    print(paste("Error reading file:", e$message))
    stop(paste("Failed to read file:", e$message))
  })
}

#' @title Calculate Mitochondrial Percentage
#' @description Calculates percentage of mitochondrial gene expression in a Seurat object
#' @param seurat_obj Seurat object to analyze
#' @param species Species to determine mitochondrial gene patterns (auto, human, mouse, etc.)
#' @return Seurat object with added percent.mt column
#' @export
Calculate_MT_Percent <- function(seurat_obj, species = NULL) {
  # Get gene names
  gene_names <- rownames(seurat_obj)
  
  # Detect species if not provided
  if (is.null(species)) {
    # Try to detect from metadata or project name
    if ("organism" %in% colnames(seurat_obj@meta.data)) {
      detected_species <- tolower(seurat_obj$organism[1])
      print(paste("Detected species from metadata:", detected_species))
      species <- detected_species
    } else {
      # Default to auto-detection
      species <- "auto"
    }
  }
  
  # Define species-specific MT patterns
  mt_patterns <- list(
    human = list(
      pattern = "^MT-",
      name = "Human (Homo sapiens)"
    ),
    mouse = list(
      patterns = c("^mt-", "^MT-", "^ENSMUSG00000064"),  # Include both patterns for mouse
      name = "Mouse (Mus musculus)"
    ),
    rat = list(
      pattern = "^(Mt-|MT-|RNO)",
      name = "Rat (Rattus norvegicus)"
    ),
    zebrafish = list(
      pattern = "^(mt-|MT-)",
      name = "Zebrafish (Danio rerio)"
    ),
    fly = list(
      pattern = "^mt:|^mt-",
      name = "Fruit Fly (Drosophila melanogaster)"
    ),
    worm = list(
      pattern = "^MTCE",
      name = "C. elegans"
    ),
    auto = list(
      patterns = c("^MT-", "^mt-", "^Mt-", "^mt:", "^MTCE", "^ENSMUSG00000064"),
      name = "Auto-detection"
    )
  )
  
  print(paste("Using pattern for species:", ifelse(is.null(species), "auto", species)))
  
  # Initialize flag
  mt_found <- FALSE
  
  # Special handling for mouse to ensure both patterns are tried
  if (species == "mouse") {
    # Try each mouse pattern separately
    for (pattern in mt_patterns$mouse$patterns) {
      print(paste("Trying mouse pattern:", pattern))
      if (any(grepl(pattern, gene_names))) {
        print(paste("Found matches with pattern:", pattern))
        tryCatch({
          seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = pattern)
          
          # Verify we got non-zero values
          if (all(seurat_obj$percent.mt == 0)) {
            print("Warning: All MT percentages are zero. Trying next pattern.")
          } else {
            print(paste("Successfully calculated MT percentage with pattern:", pattern))
            mt_found <- TRUE
            break
          }
        }, error = function(e) {
          print(paste("Error with pattern", pattern, ":", e$message))
        })
      }
    }
  } else if (species == "auto") {
    # Handle auto-detection mode - try all patterns
    for (pattern in mt_patterns$auto$patterns) {
      print(paste("Trying pattern:", pattern))
      if (any(grepl(pattern, gene_names))) {
        print(paste("Found matches with pattern:", pattern))
        tryCatch({
          seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = pattern)
          
          # Check if we got non-zero values
          if (all(seurat_obj$percent.mt == 0)) {
            print("Warning: All MT percentages are zero. Trying next pattern.")
          } else {
            print(paste("Successfully calculated MT percentage with pattern:", pattern))
            mt_found <- TRUE
            break
          }
        }, error = function(e) {
          print(paste("Error with pattern", pattern, ":", e$message))
        })
      }
    }
  } else if (species %in% names(mt_patterns)) {
    # Use species-specific pattern(s)
    if (is.list(mt_patterns[[species]]) && "patterns" %in% names(mt_patterns[[species]])) {
      # Multiple patterns
      for (pattern in mt_patterns[[species]]$patterns) {
        print(paste("Trying pattern:", pattern))
        tryCatch({
          seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = pattern)
          
          # Check if we got non-zero values
          if (all(seurat_obj$percent.mt == 0)) {
            print("Warning: All MT percentages are zero. Trying next pattern.")
          } else {
            print(paste("Successfully calculated MT percentage with pattern:", pattern))
            mt_found <- TRUE
            break
          }
        }, error = function(e) {
          print(paste("Error with pattern", pattern, ":", e$message))
        })
      }
    } else {
      # Single pattern
      pattern <- mt_patterns[[species]]$pattern
      print(paste("Using species-specific pattern:", pattern))
      
      tryCatch({
        seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = pattern)
        if (!all(seurat_obj$percent.mt == 0)) {
          print(paste("Successfully calculated MT percentage for", mt_patterns[[species]]$name))
          mt_found <- TRUE
        } else {
          print("Warning: All MT percentages are zero.")
        }
      }, error = function(e) {
        print(paste("Error calculating MT percentage for", species, ":", e$message))
      })
    }
  }
  
  # Fallback to general pattern matching if specific patterns didn't work
  if (!mt_found) {
    # Look for any MT-related feature
    generic_patterns <- c("MT", "Mt", "mt", "Mito", "mito")
    
    for (pattern in generic_patterns) {
      # For case-insensitive search
      mt_indices <- grep(pattern, gene_names, ignore.case = TRUE)
      if (length(mt_indices) > 0) {
        matching_features <- gene_names[mt_indices]
        print(paste("Found", length(matching_features), "potential mitochondrial genes using pattern:", pattern))
        print(paste("Example matches:", paste(head(matching_features, 3), collapse=", ")))
        
        tryCatch({
          seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = matching_features)
          if (!all(seurat_obj$percent.mt == 0)) {
            mt_found <- TRUE
            break
          } else {
            print("Warning: All MT percentages are zero with this pattern.")
          }
        }, error = function(e) {
          print(paste("Error calculating MT percentage with features:", e$message))
        })
      }
    }
  }
  
  # If still no MT genes found, set to 0
  if (!mt_found) {
    print("No mitochondrial genes identified. Setting percent.mt to 0.")
    seurat_obj[["percent.mt"]] <- 0
  }
  
  # Ensure no NA or infinite values
  if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
    # Check for NA or infinite values
    na_count <- sum(is.na(seurat_obj$percent.mt))
    inf_count <- sum(is.infinite(seurat_obj$percent.mt))
    
    if (na_count > 0) {
      print(paste("Found", na_count, "NA values in percent.mt. Replacing with 0."))
      seurat_obj$percent.mt[is.na(seurat_obj$percent.mt)] <- 0
    }
    
    if (inf_count > 0) {
      print(paste("Found", inf_count, "infinite values in percent.mt. Replacing with 0."))
      seurat_obj$percent.mt[is.infinite(seurat_obj$percent.mt)] <- 0
    }
    
    # Check the range of percent.mt values
    print(paste("Range of percent.mt:", min(seurat_obj$percent.mt), "to", max(seurat_obj$percent.mt)))
  } else {
    print("Warning: percent.mt column was not created successfully")
  }
  
  return(seurat_obj)
}