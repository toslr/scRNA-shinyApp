library(Seurat)
library(dplyr)
library(Matrix)
library(scCustomize)

GEO_data <- Read_GEO_Delim(data_dir = '/Users/tom/Desktop/Stanford/RA/shinyscRNA', 
                           file_suffix = '61.txt.gz')
seurat <- CreateSeuratObject(counts = GEO_data, project = "DS1")
gene_names <- rownames(GEO_data[[1]])

# Load required libraries
library(biomaRt)

# Connect to mouse ensembl database
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

# Convert gene IDs
gene_conversion <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = gene_names,
  mart = mart
)

# Merge with original data
converted_genes <- merge(
  data.frame(ensembl_gene_id = gene_names), 
  gene_conversion, 
  by = "ensembl_gene_id", 
  all.x = TRUE
)

converted_genes$external_gene_name[is.na(converted_genes$external_gene_name)] <- 
  converted_genes$ensembl_gene_id[is.na(converted_genes$external_gene_name)]


write.csv(converted_genes, file = "gene_conversion_results.csv", row.names = FALSE)
