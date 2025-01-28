library(Seurat)
library(data.table)

file_path = '/Users/tom/Desktop/Stanford/RA/GEO182846/GSE182846_series_matrix.txt'
file_path = '/Users/tom/Downloads/GSE100178_series_matrix.txt'

lines <- readLines(file_path)
sample_geos <- grep("!Sample_geo_accession", lines, value = TRUE)
sample_chars <- grep("!Sample_characteristics_ch1", lines, value = TRUE)
print(sample_chars)
sample_titles <- grep("!Sample_title", lines, value = TRUE)

geo_accessions <- strsplit(sample_geos, "\t")[[1]][-1]
treatments <- strsplit(sample_chars, "\t")[[1]][-1]  # Remove the label column
treatments <- gsub("treatment: ", "", treatments)  # Clean up treatment labels

metadata <- data.frame(
  sample_id = geo_accessions,
  treatment = treatments,
  row.names = geo_accessions
)

library(GEOquery)
gset <- getGEO("GSE100178", GSEMatrix = TRUE)
gset <- getGEO("GSE182846", GSEMatrix = TRUE)
gset <- gset[[1]]
pheno <- phenoData(gset)
phe <- pData(pheno)
