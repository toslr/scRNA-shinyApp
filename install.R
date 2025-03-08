# Function to install if package is missing
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# List of required packages
required_packages <- c(
  "shiny",
  "patchwork",
  "shinyFiles",
  "Seurat",
  "dplyr",
  "Matrix",
  "ggplot2",
  "scCustomize",
  "shinyjs",
  "DT",
  "devtools",
  "GEOquery",
  "pheatmap",
  "plotly",
  "htmlwidgets",
  "RColorBrewer"
)

# Install each package if missing
for (package in required_packages) {
  message(sprintf("Checking %s...", package))
  install_if_missing(package)
}

# Bioconductor packages if any
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("presto")) {
  devtools::install_github("immunogenomics/presto")
}

message("All required packages have been installed!")