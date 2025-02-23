# scRNA Analysis App

## Installation

1. Make sure you have R installed on your system
2. Clone this repository
3. Run the installation script:
   ```R
   source("install.R")
   ```

## Running the App

After installation, you can run the app using:
```R
shiny::runApp()
```


## App walkthrough
### What is the usage?
This app is designed to analyze preprocessed single cell RNA sequencing data. The app is designed to be user friendly and allow for easy analysis of single cell data, contains multiple renderings and allows some flexibility in the analysis.  
If you have a `_series_matrix.txt.gz` file and a folder containing count matrices in `.txt.gz` format, this app is for you!

### What is not the usage?
If you have other types of data, including raw scRNA, or other RNA sequencing data, this app does not support your type of data.

### How to install the app?

Currently supported versions include Windows and MacOS.

1. Clone or download the repository. If you downloaded it, unzip the file.

2. In a terminal, navigate to the directory of the app and launch the install.R script. This will install all necessary packages for the app to run.

### How to use the app?
1. Launch the application in a browser by clicking on the desktop icon or by running the command `shiny::runApp()` in the R console.

2. In the sidebar, input your GEO Series ID and click on Fetch. The app will download the clinical metadata from your study and display it in a table. You can then select which samples to include in your following analysis. The navigation part in the sidebar gets updated at each step and allows you to easily navigate between the different parts of your analysis.

3. In the sidebar, click on "Select Data Directory" and browse to the folder containing your count matrices. Hit "Select" and then "Read Data". This is the most computationally intense part, depending on your data size and computer this can take a few minutes. A Seurat object containing your data is created.

4. In the main panel, violin plots show quality metrics for your samples. The controls below allow you to filter out outliers (cells with too much or too little expression,...). Click Filter and run PCA to go to the next step. As with other step, you can go back later to this stage and restart the analysis from here.

5. In the dimension reduction section, a PCA is run to reduce the dimensionality of the dataset. If you're not familiar with PCA, I recommend ******** for reference. Long-story short, you want to maximize the components that discriminate your data, i.e. those with the highest standard deviation. The elbow plot shows this metrics, and the control panel below allows you to select how many dimensions to keep. There is no good answer to that, generally we aim for the elbow of the curve. You can hit Run UMAP.

6. The UMAP plot shows your data in a 2D space, considering the 2 prevailing dimensions.You can adjust the threshold below for clustering resolution. The higher the resolution, the more clusters there will be. This step uses Louvain algorithm to find an optimal number of clusters. Hit Run Clustering. If the results are not satisfying, adjust the resolution and hit Run Clustering again. 

7. In the Differential Expression section, you are given the possibility to select clusters you want to work with and rename them depending on their charateristics. Remember to hit the Update All Cluster Labels every time you make a change to the clusters. Running One vs ALl Analysis allows you to identify most expressed genes in a specific cluster, and for instance determine what it corresponds to. Running the One vs One analysis allows you to compare one cluster to another in terms of differential expression. Once the DE is run, you have the possibility to show a heatmap with the x most expressed genes from your cluster, and compare this expression across the other clusters.

8. As of now, the saving functionality is still under development. 


### 
## Requirements
- R version 4.0.0 or higher
- All packages listed in install.R

## Contact
Tom Soulaire
Stanford University