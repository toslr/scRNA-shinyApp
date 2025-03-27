# install.R

# Function to install if package is missing
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    message(sprintf("Installing %s...", package))
    install.packages(package)
    library(package, character.only = TRUE)
  } else {
    message(sprintf("%s is already installed.", package))
  }
}

# Function to install if Bioconductor package is missing
install_bioc_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    message(sprintf("Installing %s from Bioconductor...", package))
    BiocManager::install(package)
    library(package, character.only = TRUE)
  } else {
    message(sprintf("%s is already installed.", package))
  }
}

# List of required CRAN packages
required_packages <- c(
  "shiny",
  "patchwork",
  "shinyFiles",
  "Seurat",
  "dplyr",
  "Matrix",
  "ggplot2",
  "shinyjs",
  "DT",
  "devtools",
  "pheatmap",
  "plotly",
  "htmlwidgets",
  "RColorBrewer",
  "fontawesome"
)

# List of required Bioconductor packages
required_bioc_packages <- c(
  "GEOquery"
)

# Install BiocManager if missing
if (!require("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager...")
  install.packages("BiocManager")
}

# Install each CRAN package if missing
message("Checking and installing CRAN packages...")
for (package in required_packages) {
  install_if_missing(package)
}

# Install each Bioconductor package if missing
message("Checking and installing Bioconductor packages...")
for (package in required_bioc_packages) {
  install_bioc_if_missing(package)
}

# Install GitHub packages
message("Checking and installing GitHub packages...")
if (!require("presto", quietly = TRUE)) {
  message("Installing presto from GitHub...")
  devtools::install_github("immunogenomics/presto")
}

# Install any custom packages needed for scCustomize
if (!require("scCustomize", quietly = TRUE)) {
  message("Installing scCustomize from GitHub...")
  devtools::install_github("samuel-marsh/scCustomize")
}

message("All required packages have been installed!")