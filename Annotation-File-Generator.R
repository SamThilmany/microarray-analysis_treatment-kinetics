# ####################################
# Create an up-to-date annotation file
# ####################################

# Install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt", update = TRUE, ask = FALSE, checkBuilt = TRUE)

require(biomaRt)

# Get annotation table
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

annotLookup <- getBM(
  mart = ensembl,
  attributes = c(
    'agilent_sureprint_g3_ge_8x60k_v2',
    'ensembl_gene_id',
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