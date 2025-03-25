# shinyscRNA: Single-Cell RNA Analysis Tool 

A comprehensive Shiny application for analyzing single-cell RNA sequencing data with an intuitive modular workflow.

## Table of Contents
- [Installation](#installation)
- [User Guide](#user-guide)
  - [Metadata Import](#metadata-import)
  - [Data Input](#data-input)
  - [Quality Control](#quality-control)
  - [Dimension Reduction](#dimension-reduction)
  - [Differential Expression](#differential-expression)
  - [Saving and Loading](#saving-and-loading)
- [Code Architecture](#code-architecture)
  - [App Structure](#app-structure)
  - [Module System](#module-system)
  - [Code Documentation](#code-documentation)

## Installation

### Prerequisites
- R (version 4.0.0 or higher)
- RStudio (recommended for development)

### Installation Steps

1. Clone the repository:
```bash
git clone https://github.com/your-username/shinyscRNA.git
cd shinyscRNA
```

2. Install required R packages:
```r
source("install.R")
```

This will install all dependencies required by the application, including:
- shiny
- Seurat
- patchwork
- dplyr
- ggplot2
- plotly
- DT
- pheatmap
- GEOquery
- RColorBrewer
- and other needed packages

3. Launch the application:
```r
shiny::runApp()
```

## User Guide

The shinyscRNA app provides a step-by-step workflow for analyzing single-cell RNA sequencing data. Follow these steps to complete your analysis:

### Metadata Import

1. **GEO Series Input**:
   - Enter a GEO Series ID (e.g., "GSE123456") in the input field.
   - Click "Fetch GEO Metadata".
   - The application will load sample metadata from the Gene Expression Omnibus.
   - After loading, you'll see a table showing all samples in the dataset.
   - Use the checkboxes to select which samples to include in your analysis.

### Data Input

1. **Data Selection**:
   - After selecting samples from the metadata step, click "Select Data Directory".
   - Navigate to the folder containing your expression data files.
   - The app expects files matching the GEO IDs selected in the previous step.
   - Click "Read Data" to import the expression matrices into the application.

### Quality Control

1. **QC Metrics Inspection**:
   - After data input, the QC section will display violin plots showing:
     - Number of features (genes) per cell
     - Number of UMI counts per cell
     - Percentage of mitochondrial reads
   - These plots help identify low-quality cells.

2. **Filter Parameters**:
   - Adjust these parameters based on your QC plots:
     - Minimum Features: Cells with fewer features than this will be filtered out.
     - Maximum Features: Cells with more features than this will be filtered out (removes potential doublets).
     - Maximum MT %: Cells with higher mitochondrial percentage will be filtered out (typically dying cells).
   - Click "Filter and run PCA" to apply these filters and normalize the data.

### Dimension Reduction

1. **PCA Elbow Plot**:
   - After QC filtering, an elbow plot will show variance explained by each principal component.
   - A suggested number of dimensions will be highlighted based on the elbow point.
   - You can adjust this number if needed.

2. **UMAP Visualization**:
   - Click "Confirm Dimensions" to proceed with the selected number of PCs.
   - Enter a clustering resolution (0.4-0.8 is common; higher values create more clusters).
   - Click "Run Clustering" to perform graph-based clustering.
   - The UMAP plots will show cells colored by:
     - Sample origin (left panel default)
     - Cluster assignment (right panel default)
   - Toggle between 2D and 3D visualizations using the buttons.
   - Search for specific genes to visualize their expression on the UMAP.
   - Download plots using the "Save Plot" buttons.

3. **Cluster Management**:
   - In the sidebar, you can rename clusters based on marker genes.
   - Enable/disable clusters for downstream analysis.
   - Click "Save Labels" to store your cluster annotations.

### Differential Expression

1. **DE Analysis Options**:
   - **One vs All Analysis**: Compare gene expression in one cluster against all others.
   - **One vs One Analysis**: Compare gene expression between two specific clusters.
   - **General Cluster Map**: Create a heatmap showing top marker genes for each cluster.

2. **Results Visualization**:
   - **Volcano Plot**: Shows significantly up/down-regulated genes.
   - **Gene Table**: Displays differential expression statistics for all genes.
   - **Heatmap**: Shows expression patterns of top differentially expressed genes.
   - All visualizations can be saved for inclusion in reports.

### Saving and Loading

1. **Save Analysis**:
   - Click "Save Analysis" in the sidebar.
   - Enter a name for your analysis.
   - This saves the entire state of your analysis, including:
     - Loaded data
     - QC parameters
     - Clustering results
     - Differential expression results

2. **Load Analysis**:
   - Click "Load Analysis" in the sidebar.
   - Select a previously saved analysis from the dropdown.
   - The application will restore all states from the saved session.

## Code Architecture

### App Structure

The application follows a modular design pattern with these main components:

```
shinyscRNA/
├── app.R                  # Main application entry point
├── R/
│   ├── ui.R               # User interface definition
│   ├── server.R           # Server function definition
│   ├── modules/           # Functional modules
│   │   ├── data_input_module.R
│   │   ├── metadata_module.R
│   │   ├── qc_module.R
│   │   ├── dimension_reduction_module.R
│   │   ├── de_analysis_module/
│   │   │   ├── de_analysis_module.R
│   │   │   ├── cluster_utils.R
│   │   │   ├── de_computation.R
│   │   │   ├── visualization.R
│   │   │   └── ui_components.R
│   │   ├── cluster_management_module.R
│   │   └── save_load_module.R
│   ├── server/
│   │   ├── navigation.R
│   │   ├── observers.R
│   │   └── sections.R
│   └── modules/*_utils/   # Utility functions for modules
└── www/                   # Web assets
    ├── script.js
    └── styles.css
```

### Module System

Each major functionality is encapsulated in a module following the Shiny module pattern:

1. **Metadata Module**:
   - Imports sample metadata from GEO
   - Manages sample selection

2. **Data Input Module**:
   - Handles file selection and reading
   - Creates initial Seurat objects

3. **QC Module**:
   - Generates QC metrics visualizations
   - Applies filtering and normalization

4. **Dimension Reduction Module**:
   - Performs PCA and determines optimal dimensions
   - Runs UMAP and clustering
   - Visualizes dimensionality reduction results

5. **Cluster Management Module**:
   - Allows renaming and activation/deactivation of clusters

6. **DE Analysis Module**:
   - Performs differential expression analysis between clusters
   - Generates volcano plots, tables, and heatmaps

7. **Save/Load Module**:
   - Handles saving and loading analysis states

### Code Documentation

The code uses function-level documentation with descriptive names. For comprehensive documentation:

1. **Adding Roxygen-style Documentation**:
   - Consider adding Roxygen2 documentation to functions:
   ```r
   #' Process Seurat QC Filtering
   #' 
   #' This function applies quality control filtering to a Seurat object based on specified parameters.
   #' 
   #' @param seurat_obj A Seurat object to filter
   #' @param min_feature Minimum number of features (genes) per cell
   #' @param max_feature Maximum number of features per cell
   #' @param max_mt Maximum percentage of mitochondrial reads
   #' @return A filtered and normalized Seurat object with PCA run
   #' 
   processQCFiltering <- function(seurat_obj, min_feature, max_feature, max_mt) {
      # Function implementation
   }
   ```

2. **Module Design**:
   - Each module follows the standard Shiny module pattern:
     - `moduleNameUI(id)` function for UI components
     - `moduleNameServer(id, ...)` function for server logic
     - Returns reactive expressions for inter-module communication

3. **Extending the Application**:
   - Add new modules by following the existing pattern
   - Update `app.R` to include new module source files
   - Update the UI and server functions to integrate new modules

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the [MIT License](LICENSE).

## Acknowledgements

- [Seurat](https://satijalab.org/seurat/) - The core single-cell analysis framework
- [Shiny](https://shiny.rstudio.com/) - The web application framework
- Stanford University, Zuchero Lab