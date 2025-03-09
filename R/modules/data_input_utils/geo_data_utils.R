# R/utils/geo_data_utils.R

findMatchingFiles <- function(all_files, selected_samples) {
  selected_files <- character(0)
  file_to_gsm <<- list()  # Using global assignment to make this accessible
  
  for(gsm in selected_samples) {
    pattern <- paste0("^", gsm)
    matches <- grep(pattern, all_files, value=TRUE)
    print(paste("Looking for files matching GSM:", gsm))
    print(paste("Found matches:", paste(matches, collapse=", ")))
    selected_files <- c(selected_files, matches)
    file_to_gsm[matches] <- gsm
  }
  
  if(length(selected_files) == 0) {
    stop("No matching files found for selected GSMs")
  }
  
  return(selected_files)
}

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
      print(paste("Creating Seurat object for GSM:", gsm))
      
      seurat <- CreateSeuratObject(counts = data[[1]], project = gsm)
      print("Created Seurat object")
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