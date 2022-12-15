# Reference: https://support.bioconductor.org/p/9147231/

# #####################
# General configuration
# #####################

# Set the base directory
baseDir <- getwd()

# Define which experiment(s) should be analyzed
experiments = c(1, 2)

# Define which min logFC values should be analyzed
logFC_values <- c(2,4)

# Define which p values should be considered
p_values <- c(0.05)

# Define the degree for the polynomial fit
polyDegree <- 3

# Define name of negative control treatment
treat_NC <- "NC"

# Define the graphics dimensions
graphics_dimensions <- c(1200, 700)

# Create a folder for graphics
graphicsDir <- paste0(baseDir, '/graphics')

if (!dir.exists(graphicsDir)) {
  dir.create(graphicsDir)
}

# Create a folder for results
resultsDir <- paste0(baseDir, '/results')

if (!dir.exists(resultsDir)) {
  dir.create(resultsDir)
}

# Targets file
targetsFile <- 'Targets.tsv'

# Prevent scientific notation
options(scipen = 99)

# Install the Bioconductor Manager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor Packages
BiocManager::install("limma", update = TRUE, ask = FALSE, force = TRUE, checkBuilt = TRUE)
require(limma)

# Load R packages
require(statmod)
require(stringr)
require(gplots)
require(tidyr)
require(dplyr)
require(VennDiagram)



# ##################################
# Read or create the annotation file
# ##################################

# Create a folder for annotation files
annotationDir <- paste0(baseDir, "/annotation-files")

if (!dir.exists(annotationDir)) {
  dir.create(annotationDir)
}

setwd(annotationDir)

tryCatch(
  expr = {
    currentMonth <- format(Sys.Date(), format="%Y-%m")
    
    annotationFile <- paste0(
      'Human_agilent_sureprint_g3_ge_8x60k_v2_', 
      gsub("-", "_", as.character(currentMonth)), 
      '.tsv'
    )
    
    if (!file.exists(annotationFile)) {
      source(paste0(baseDir, '/', 'Annotation-File-Generator.R'), local = TRUE)
    }
  },
  finally = {
    setwd(baseDir)
  }
)



# ################
# Run the Analysis
# ################

dim_eset <- list()

for (experiment in experiments) {
  graphicsDirExp <- paste0(graphicsDir, '/exp_', experiment)
  if (dir.exists(graphicsDirExp)) {
    unlink(graphicsDirExp, recursive = TRUE)
  }
  dir.create(graphicsDirExp)
  dir.create(paste0(graphicsDirExp, '/data-processing'))

  resultsDirExp <- paste0(resultsDir, '/exp_', experiment)
  if (dir.exists(resultsDirExp)) {
    unlink(resultsDirExp, recursive = TRUE)
  }
  dir.create(resultsDirExp)
  
  # Capture the output to a log file
  sink(file = paste0(resultsDirExp, '/log.txt'), append = TRUE, type = c('output', 'message'), split = TRUE)
  
  source(paste0(baseDir, '/', '_analysis-script.R'), local = TRUE)
  
  # Remove tmp variables
  rm(list = ls(pattern = '^tmp_'))
  
  # End the output capture
  sink(file = NULL, split = FALSE)
  
  # Close all connections
  closeAllConnections()
}