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
library(plotly)
library(htmlwidgets)
library(RColorBrewer)
library(fontawesome)

# Source metadata module
source("R/modules/metadata_module.R")

# Source data input
source("R/modules/data_input_module.R")
source("R/modules/data_input_utils/geo_data_utils.R")

# Source qc module
source("R/modules/qc_module.R")

# Source dimred
source("R/modules/dimension_reduction_module.R")
source("R/modules/dimension_reduction_utils/dimred_utils.R")

# Source de_module
source("R/modules/de_analysis_module.R")
source("R/modules/de_analysis_utils/cluster_utils.R")
source("R/modules/de_analysis_utils/de_computation.R")
source("R/modules/de_analysis_utils/visualization.R")
source("R/modules/de_analysis_utils/ui_components.R")

# Source save/load module
source("R/modules/save_load_module.R")

# Source management modules
source("R/modules/cluster_management_module.R")
source("R/modules/sample_management_module.R")
source("R/modules/condition_management_module.R")

# Source server
source("R/server/observers.R")
source("R/server/navigation.R")
source("R/server/sections.R")
source("R/server.R")
source("R/ui.R")

shinyApp(ui = buildUI(), server = buildServer())