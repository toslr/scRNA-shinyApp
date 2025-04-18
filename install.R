# install.R

# Function to install if package is missing
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    message(sprintf("Installing %s...", package))
    install.packages(package, dependencies = TRUE, ask = FALSE)
    library(package, character.only = TRUE)
  } else {
    message(sprintf("%s is already installed.", package))
  }
}

# Function to install if Bioconductor package is missing
install_bioc_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    message(sprintf("Installing %s from Bioconductor...", package))
    BiocManager::install(package, ask = FALSE, update = FALSE)
    library(package, character.only = TRUE)
  } else {
    message(sprintf("%s is already installed.", package))
  }
}

# Function to install if GitHub package is missing
install_github_if_missing <- function(package, repo) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    message(sprintf("Installing %s from GitHub (%s)...", package, repo))
    devtools::install_github(repo, dependencies = TRUE, upgrade = "never")
    library(package, character.only = TRUE)
  } else {
    message(sprintf("%s is already installed.", package))
  }
}

# List of required CRAN packages
required_packages <- c(
  # Core Shiny packages
  "shiny",
  "shinyjs",
  "shinyFiles",
  "DT",
  
  # Data manipulation
  "dplyr",
  "Matrix",
  "tibble",
  "tidyr",
  
  # Visualization
  "ggplot2",
  "patchwork",
  "pheatmap",
  "plotly",
  "htmlwidgets",
  "RColorBrewer",
  "fontawesome",
  "ggpubr",
  
  # Development tools
  "devtools",
  "remotes",
  
  # Single-cell analysis
  "Seurat",
  "hdf5r",
  
  # Additional utilities
  "jsonlite",
  "R.utils",
  "future",
  "future.apply"
)

# List of required Bioconductor packages
required_bioc_packages <- c(
  "GEOquery",
  "biomaRt",
  "SingleCellExperiment",
  "AnnotationDbi"
  #"org.Hs.eg.db",    # Human annotations
  #"org.Mm.eg.db",    # Mouse annotations
  #"org.Rn.eg.db",    # Rat annotations
  #"org.Dr.eg.db",    # Zebrafish annotations
  #"org.Dm.eg.db",    # Drosophila annotations
  #"org.Ce.eg.db"     # C. elegans annotations
)

# GitHub packages
github_packages <- list(
  "scCustomize" = "samuel-marsh/scCustomize",
  "presto" = "immunogenomics/presto"
)

# Create a log file
log_file <- "install_log.txt"
cat("Installation Log\n", file = log_file, append = FALSE)
cat(paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n"), file = log_file, append = TRUE)

# Function to log messages
log_message <- function(msg) {
  message(msg)
  cat(paste(format(Sys.time(), "%H:%M:%S"), "-", msg, "\n"), file = log_file, append = TRUE)
}

# Install BiocManager if missing
if (!require("BiocManager", quietly = TRUE)) {
  log_message("Installing BiocManager...")
  install.packages("BiocManager", dependencies = TRUE, ask = FALSE)
}
options(BiocManager.check_repositories = FALSE)

# Set up parallel processing for faster installations
log_message("Setting up parallel processing for installations...")
if (requireNamespace("parallel", quietly = TRUE)) {
  ncores <- min(4, parallel::detectCores())
  options(Ncpus = ncores)
  log_message(paste("Using", ncores, "cores for parallel installation"))
}

# Install each CRAN package if missing
log_message("Checking and installing CRAN packages...")
for (package in required_packages) {
  tryCatch({
    install_if_missing(package)
  }, error = function(e) {
    log_message(paste("Error installing", package, ":", e$message))
  })
}

# Install each Bioconductor package if missing
log_message("Checking and installing Bioconductor packages...")
for (package in required_bioc_packages) {
  tryCatch({
    install_bioc_if_missing(package)
  }, error = function(e) {
    log_message(paste("Error installing", package, ":", e$message))
  })
}

# More detailed error handling for GitHub packages
log_message("Checking and installing GitHub packages...")
for (package_name in names(github_packages)) {
  repo <- github_packages[[package_name]]
  tryCatch({
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      log_message(paste("Installing", package_name, "from GitHub:", repo))
      
      # Ensure remotes is available
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes", dependencies = TRUE)
        library(remotes)
      }
      
      # Force installation with no interactive prompts
      remotes::install_github(
        repo,
        dependencies = TRUE,
        upgrade = "never",
        quiet = FALSE,  # Set to FALSE to see more output
        force = TRUE,
        INSTALL_opts = c("--no-multiarch", "--no-test-load")
      )
      
      # Check if installation succeeded
      if (require(package_name, character.only = TRUE, quietly = TRUE)) {
        log_message(paste("Successfully installed", package_name))
      } else {
        log_message(paste("WARNING: Failed to load", package_name, "after installation"))
      }
    } else {
      log_message(paste(package_name, "is already installed"))
    }
  }, error = function(e) {
    log_message(paste("ERROR installing", package_name, ":", e$message))
  })
}


log_message("All required packages have been installed!")
cat(paste("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = log_file, append = TRUE)

# Print summary
message("\n=== Installation Summary ===")
message("Installation log saved to: ", normalizePath(log_file))
message("To start the application, run: shiny::runApp()")
