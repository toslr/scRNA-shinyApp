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

# Source all module files
source("R/ui.R")
source("R/server.R")
source("R/server/observers.R")
source("R/server/navigation.R")
source("R/server/sections.R")
source("R/modules/data_input_module.R")
source("R/modules/qc_module.R")
source("R/modules/dimension_reduction_module.R")
source("R/modules/de_analysis_module.R")

shinyApp(ui = buildUI(), server = buildServer())