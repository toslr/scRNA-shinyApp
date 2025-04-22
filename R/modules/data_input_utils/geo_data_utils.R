#' @title Read GEO Delimited File
#' @description Reads a GEO data file in delimited format
#' @param data_dir Directory containing the file
#' @param file_suffix Filename or suffix to read
#' @return A list containing a sparse matrix of gene counts
#' @export
Read_GEO_Delim <- function(data_dir, file_suffix) {
  # Find the actual file path
  if (grepl("\\.txt\\.gz$", file_suffix)) {
    file_path <- file.path(data_dir, file_suffix)
  } else {
    # Look for a file matching the GSM pattern
    files <- list.files(data_dir, pattern = paste0("^", file_suffix, ".*\\.txt\\.gz$"), full.names = TRUE)
    if (length(files) == 0) {
      stop(paste("No matching .txt.gz files found for", file_suffix))
    }
    file_path <- files[1]
  }
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Read the first few lines to determine format
  con <- gzfile(file_path, "r")
  header <- readLines(con, n = 5)
  close(con)
  
  # Determine separator and skip based on header
  if (any(grepl("^\\s*ENSMUSG|^\\s*ENSG", header))) {
    # Standard GEO format with gene IDs in first column
    separator <- "\t"
    skip_lines <- 0
  } else {
    # Try to detect format
    if (any(grepl("\t", header))) {
      separator <- "\t"
    } else if (any(grepl(",", header))) {
      separator <- ","
    } else {
      separator <- " "
    }
    
    # Check for header lines to skip
    skip_lines <- sum(grepl("^#|^\"#", header))
  }
  
  # Read the data
  data <- tryCatch({
    df <- read.delim(file_path, sep = separator, skip = skip_lines, 
                     header = TRUE, row.names = 1, check.names = FALSE)
    
    # Convert to sparse matrix
    sparse_mat <- as(as.matrix(df), "dgCMatrix")
    list(sparse_mat)
  }, error = function(e) {
    stop(paste("Error reading file:", e$message))
  })
  
  return(data)
}

#' @title Find Matching GEO Files
#' @description Finds files that match the selected GEO Sample accessions (GSM).
#' @param all_files Character vector of all available files in the directory
#' @param selected_samples Character vector of GSM IDs to be matched
#' @return Character vector of matching file names
#' @export
findMatchingFiles <- function(all_files, selected_samples) {
  selected_files <- character(0)
  file_to_gsm <- list()  # Using assignment to make this accessible
  
  for(gsm in selected_samples) {
    pattern <- paste0("^", gsm)
    matches <- grep(pattern, all_files, value=TRUE)
    selected_files <- c(selected_files, matches)
    file_to_gsm[matches] <- gsm
  }
  
  if(length(selected_files) == 0) {
    stop("No matching files found for selected GSMs")
  }
  
  return(selected_files)
}

#' @title Add Metadata to Seurat Object
#' @description Adds GEO metadata to a Seurat object based on the sample accession.
#' @param seurat_obj The Seurat object to add metadata to
#' @param gsm The GEO Sample accession (GSM) ID
#' @param metadata_function Function to retrieve metadata for the samples
#' @return The Seurat object with added metadata
#' @export
addMetadata <- function(seurat_obj, gsm, metadata_function) {
  if (!is.null(metadata_function)) {
    metadata_data <- metadata_function()
    if (!is.null(metadata_data)) {
      sample_meta <- metadata_data[metadata_data$geo_accession == gsm, ]
      
      for(col in setdiff(colnames(sample_meta), "geo_accession")) {
        seurat_obj[[col]] <- sample_meta[[col]][1]
      }
    }
  }
  
  return(seurat_obj)
}