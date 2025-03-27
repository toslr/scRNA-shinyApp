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

#' @title Create Seurat Objects from GEO Files
#' @description Creates Seurat objects from GEO data files, adding metadata and gene mappings.
#' @param dir_path Path to the directory containing data files
#' @param selected_files Character vector of selected file names to process
#' @param selected_samples Character vector of GSM IDs corresponding to the files
#' @param gene_mapping Named vector for mapping Ensembl IDs to gene names
#' @param metadata_function Function to retrieve metadata for the samples
#' @return List of Seurat objects, one for each sample
#' @export
createSeuratObjects <- function(dir_path, selected_files, selected_samples, gene_mapping, metadata_function) {
  seurat_objects <- list()
  count = 0
  
  for(file in selected_files) {
    count = count + 1
    incProgress(0.9/(length(selected_files)+1), 
                detail = paste("Reading file", count, "of", length(selected_files)))
    
    # Read the file
    data <- Read_GEO_Delim(data_dir = dir_path, 
                           file_suffix = file)
    
    if(length(data) > 0 && !is.null(data[[1]])) {
      # Create Seurat object for this sample
      gsm <- file_to_gsm[[file]]
      seurat <- CreateSeuratObject(counts = data[[1]], project = gsm)
      seurat$sample <- gsm
      
      # Add GEO metadata if available
      addMetadata(seurat, gsm, metadata_function)
      
      # Add percent.mt
      seurat[["percent.mt"]] <- PercentageFeatureSet(seurat,
                                                     pattern = "^ENSMUSG00000064")
      
      # Store gene mapping
      seurat@misc$gene_mapping <- gene_mapping
      
      # Store in list
      seurat_objects[[gsm]] <- seurat
    }
  }
  
  return(seurat_objects)
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