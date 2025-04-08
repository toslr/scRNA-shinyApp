#' @title Get Gene Mapping for Species
#' @description Loads or creates gene ID mapping for a specific species
#' @param species Species name (e.g., "human", "mouse")
#' @param cache_dir Directory to cache gene mapping files
#' @return A named vector mapping Ensembl IDs to gene symbols
#' @export
get_gene_mapping <- function(species = "mouse", cache_dir = "gene_mappings") {
  # Ensure cache directory exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  # Define file path for cached mapping
  cache_file <- file.path(cache_dir, paste0(species, "_gene_mapping.rds"))
  
  # Check if mapping already exists
  if (file.exists(cache_file)) {
    message("Loading cached gene mapping for ", species)
    return(readRDS(cache_file))
  }
  
  # If not cached, create mapping
  message("Creating gene mapping for ", species, "...")
  
  # Install biomaRt if needed
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("Installing biomaRt package...")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("biomaRt")
  }
  
  # Define species-specific Ensembl dataset
  dataset_map <- list(
    human = "hsapiens_gene_ensembl",
    mouse = "mmusculus_gene_ensembl",
    rat = "rnorvegicus_gene_ensembl",
    zebrafish = "drerio_gene_ensembl",
    fly = "dmelanogaster_gene_ensembl",
    worm = "celegans_gene_ensembl"
  )
  
  # Default to mouse if species not recognized
  if (!species %in% names(dataset_map)) {
    warning("Species '", species, "' not recognized. Using mouse as default.")
    species <- "mouse"
  }
  
  # Get Ensembl dataset
  dataset <- dataset_map[[species]]
  
  tryCatch({
    # Connect to Ensembl
    mart <- biomaRt::useMart("ensembl", dataset = dataset)
    
    # Get gene ID to symbol mapping
    gene_data <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      mart = mart
    )
    
    # Convert to named vector
    gene_mapping <- setNames(gene_data$external_gene_name, gene_data$ensembl_gene_id)
    
    # Cache the result
    saveRDS(gene_mapping, cache_file)
    
    message("Created gene mapping with ", length(gene_mapping), " genes")
    return(gene_mapping)
  }, error = function(e) {
    warning("Error creating gene mapping: ", e$message)
    # Return empty mapping if failed
    return(setNames(character(0), character(0)))
  })
}

#' @title Detect Species from Data
#' @description Attempts to detect species from gene IDs or metadata
#' @param seurat_obj Seurat object or gene IDs vector
#' @return Character string with detected species name
#' @export
detect_species <- function(seurat_obj) {
  # Get gene IDs to analyze
  gene_ids <- if (class(seurat_obj)[1] == "Seurat") {
    rownames(seurat_obj)
  } else {
    seurat_obj  # Assume it's already a vector of gene IDs
  }
  
  # Sample some genes for pattern matching (use up to 1000 genes)
  sample_size <- min(1000, length(gene_ids))
  sampled_genes <- sample(gene_ids, sample_size)
  
  # Define patterns to check
  patterns <- list(
    human = list(
      pattern = "^ENSG0|^HGNC:|^MT-",
      name = "human"
    ),
    mouse = list(
      pattern = "^ENSMUSG0|^ENSMUS|^mt-",
      name = "mouse"
    ),
    rat = list(
      pattern = "^ENSRNOG0|^RGD:|^Mt-",
      name = "rat"
    ),
    zebrafish = list(
      pattern = "^ENSDARG0|^ZDB-",
      name = "zebrafish"
    ),
    fly = list(
      pattern = "^FBgn|^CG\\d+",
      name = "fly"
    ),
    worm = list(
      pattern = "^WBGene|^MTCE",
      name = "worm"
    )
  )
  
  # Check species-specific patterns
  matches <- sapply(patterns, function(p) {
    sum(grepl(p$pattern, sampled_genes))
  })
  
  # Get the species with the most matches
  best_match <- names(which.max(matches))
  best_count <- max(matches)
  
  # If there's a clear winner and it has a significant number of matches
  if (best_count > 10 && best_count > 0.1 * sample_size) {
    message("Detected species: ", patterns[[best_match]]$name, " (",
            best_count, " matching genes out of ", sample_size, " sampled)")
    return(patterns[[best_match]]$name)
  } else {
    message("Could not confidently detect species. Using 'mouse' as default.")
    return("mouse")
  }
}