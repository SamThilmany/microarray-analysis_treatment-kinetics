# ####################################
# Create an up-to-date annotation file
# ####################################

# Install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

require(biomaRt)

# Get annotation table
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'agilent_sureprint_g3_ge_8x60k_v2',
    'wikigene_description',
    'ensembl_gene_id',
    'entrezgene_id',
    'gene_biotype',
    'external_gene_name'
  )
)

# Export annotation table as tsv file
today <- format(Sys.Date(), format="%Y-%m")

write.table(
  annotLookup,
  paste0('Human_agilent_sureprint_g3_ge_8x60k_v2_', gsub("-", "_", as.character(today)), '.tsv'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
)