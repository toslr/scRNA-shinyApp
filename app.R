# app.R

# Load required libraries
library(shiny)
library(patchwork)
library(shinyFiles)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scCustomize)
library(shinyjs)
library(DT)
library(GEOquery)
library(pheatmap)

# Source de_module
source("R/modules/de_analysis_module/de_analysis_module.R")
source("R/modules/de_analysis_module/cluster_utils.R")
source("R/modules/de_analysis_module/de_computation.R")
source("R/modules/de_analysis_module/visualization.R")
source("R/modules/de_analysis_module/ui_components.R")

#source other modules
source("R/modules/metadata_module.R")
source("R/modules/data_input_module.R")
source("R/modules/qc_module.R")
source("R/modules/dimension_reduction_module.R")
#source("R/modules/save_load_module.R")

# Source server files
source("R/server/observers.R")
source("R/server/navigation.R")
source("R/server/sections.R")
source("R/server.R")
source("R/ui.R")

shinyApp(ui = buildUI(), server = buildServer())