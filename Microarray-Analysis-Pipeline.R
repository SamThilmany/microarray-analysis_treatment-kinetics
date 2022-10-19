# #####################
# General configuration
# #####################

baseDir <- getwd()
print(baseDir)
targetsFile <- 'targets.tsv'

options(scipen = 99) # prevent scientific notation

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma"))
library(limma)

require(statmod)
require(stringr)
require(gplots)
require(tidyr)



# ##################################
# Read or create the annotation file
# ##################################

# Create a folder for annotation files
dir_annotation <- paste0(baseDir, "/annotation-files")

if (!dir.exists(dir_annotation)) {
  dir.create(dir_annotation)
}

setwd(dir_annotation)

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
# Read in the data
# ################

# Targets
targetinfo <- readTargets(targetsFile, row.names = 'Name')
targetinfo[order(targetinfo$Name),]

# Select a subset for a specific experiment
targetinfo = subset(targetinfo, Experiment == 1)

# Converts the raw data to an EListRaw object
wtAgilent.GFilter <- function(qta) { qta[,"gIsPosAndSignif"] }
eset <- read.maimages(
  targetinfo,
  source = 'agilent.median',
  green.only = TRUE,
  path = "data",
  other.columns = 'gIsWellAboveBG',
  wt.fun = wtAgilent.GFilter
)

colnames(eset) <- row.names(targetinfo)
eset <- eset[, order(colnames(eset))]

# Add the spot type
spotTypes <- readSpotTypes(file = 'SpotTypes.tsv')
eset$genes$Status <- controlStatus(spotTypes, eset)



# ###################
# Annotate the probes
# ###################

annotLookup <- read.csv(
  paste0(dir_annotation, '/', annotationFile),
  header = TRUE,
  sep = '\t',
  stringsAsFactors = FALSE
)

colnames(annotLookup)[1] <- 'AgilentID'

annotLookup <- annotLookup[which(annotLookup$AgilentID %in% eset$genes$ProbeName),]
annotLookup <- annotLookup[match(eset$genes$ProbeName, annotLookup$AgilentID),]
table(eset$genes$ProbeName == annotLookup$AgilentID) # check that annotations are aligned

eset$genes$AgilentID <- annotLookup$AgilentID
eset$genes$wikigene_description <- annotLookup$wikigene_description
eset$genes$ensembl_gene_id <- annotLookup$ensembl_gene_id
eset$genes$entrezgene <- annotLookup$entrezgene
eset$genes$gene_biotype <- annotLookup$gene_biotype
eset$genes$external_gene_name <- annotLookup$external_gene_name



# #############################
# Perform background correction
# #############################

eset <- backgroundCorrect(eset, method = 'normexp')
plotDensities(eset, legend = FALSE, main = 'Density plot with background-corrected data')



# ##############
# Normalize data
# ##############

eset <- normalizeBetweenArrays(eset, method = 'quantile')
plotDensities(eset, legend = FALSE, main = 'Density plot with normalized, background-corrected data')



# #######
# MA-Plot
# #######

# plotMDS(eset)



# ##########################
# Filter out control probes, 
# those with no symbol, 
# and those that fail
# ##########################

Control <- eset$genes$ControlType != 0
NoSymbol <- is.na(eset$genes$external_gene_name)
IsExpr <- rowSums(eset$other$gIsWellAboveBG > 0) == 20

eset <- eset[!Control & !NoSymbol & IsExpr, ]



# ##########################################
# Remove annotation columns no longer needed
# ##########################################

eset$genes <- eset$genes[,c(
  'ProbeName','wikigene_description','ensembl_gene_id','entrezgene','gene_biotype','external_gene_name'
)]



# #####################
# Array quality weights
# #####################

array.weights <- arrayWeights(eset)
barplot(array.weights, xlab = "Array", ylab = "Weight", main = "Array weights", col = "white", las = 2)
abline(h = 1, lwd = 1, lty = 2)



# #################
# Expression values
# #################

boxplot(log2(eset$E), range = 0, ylab = "log2 intensity")


# ###############################
# Create a linear fit of the data
# ###############################

# array_batch <- factor(targetinfo$Array_Batch)
# treatment <- factor(targetinfo$Treatment)
# time <- factor(targetinfo$Time)
# 
# design <- model.matrix(~ array_batch + treatment + time)
# 
# fit <- lmFit(eset, design = design, weights = array.weights)
# fit <- eBayes(fit)



# ###########################
# Create a Student's Q-Q Plot
# ###########################

# qqt(fit$t, df = fit$df.prior + fit$df.residual, pch = 16, cex = 0.2)
# abline(0,1)



# ################
# Create a boxplot
# ################

# cols <- eset$targets$CellType
# cols[cols == names$SH_SY5Y] <- "steelblue1"
# cols[cols == names$SK_N_BE_2] <- "steelblue2"
# cols[cols == names$Kelly] <- "steelblue3"
# cols[cols == names$ReNcell_VM] <- "steelblue4"
# cols[cols == names$Adult_Brain] <- "springgreen1"
# cols[cols == names$Hippocampus] <- "springgreen2"
# cols[cols == names$Insula] <- "springgreen3"
# 
# boxplot(eset$E ~ col(eset$E), names = rownames(eset$targets), col = cols, xlab = "Cell Type", ylab = "E-values", las = 2, mar = c(15, 4, 4, 2) + 0.1)



# ##############################
# Pairwise comparison of FC data
# ##############################

# S3 object containing the difference strings
# differences <- list(Adult_Brain = '', Hippocampus = '', Insula = '')
# class(differences) <- "Difference strings for FC comparison of different brain regions and cell lines"
# 
# differences$Adult_Brain <- list(SH_SY5Y = '', SK_N_BE_2 = '', Kelly = '', ReNcell_VM = '')
# class(differences$Adult_Brain) <- "Difference strings for the adult brain"
# 
# differences$Hippocampus <- list(SH_SY5Y = '', SK_N_BE_2 = '', Kelly = '', ReNcell_VM = '')
# class(differences$Hippocampus) <- "Difference strings for the hippocampus"
# 
# differences$Insula <- list(SH_SY5Y = '', SK_N_BE_2 = '', Kelly = '', ReNcell_VM = '')
# class(differences$Insula) <- "Difference strings for the insula"

# S3 object containing the contrast matrices
# contrast_matrices <- list(SH_SY5Y = '', SK_N_BE_2 = '', Kelly = '', ReNcell_VM = '')
# class(contrast_matrices) <- "Constrast matrices for FC comparison different brain regions and cell lines"

# S3 object containing the fits
# fits <- list(SH_SY5Y = '', SK_N_BE_2 = '', Kelly = '', ReNcell_VM = '')
# class(fits) <- "Fits for FC comparison different brain regions and cell lines"

# S3 object containing the results
# results <- list(SH_SY5Y = '', SK_N_BE_2 = '', Kelly = '', ReNcell_VM = '')
# class(results) <- "Results for FC comparison between different brain regions and cell lines"


# Adult brain vs cell lines
# differences$Adult_Brain$SH_SY5Y <- paste0(names$Adult_Brain,"-",names$SH_SY5Y)
# differences$Adult_Brain$SK_N_BE_2 <- paste0(names$Adult_Brain,"-",names$SK_N_BE_2)
# differences$Adult_Brain$Kelly <- paste0(names$Adult_Brain,"-",names$Kelly)
# differences$Adult_Brain$ReNcell_VM <- paste0(names$Adult_Brain,"-",names$ReNcell_VM)

# Hippocampus vs cell lines
# differences$Hippocampus$SH_SY5Y <- paste0(names$Hippocampus,"-",names$SH_SY5Y)
# differences$Hippocampus$SK_N_BE_2 <- paste0(names$Hippocampus,"-",names$SK_N_BE_2)
# differences$Hippocampus$Kelly <- paste0(names$Hippocampus,"-",names$Kelly)
# differences$Hippocampus$ReNcell_VM <- paste0(names$Hippocampus,"-",names$ReNcell_VM)

# Insula vs cell lines
# differences$Insula$SH_SY5Y <- paste0(names$Insula,"-",names$SH_SY5Y)
# differences$Insula$SK_N_BE_2 <- paste0(names$Insula,"-",names$SK_N_BE_2)
# differences$Insula$Kelly <- paste0(names$Insula,"-",names$Kelly)
# differences$Insula$ReNcell_VM <- paste0(names$Insula,"-",names$ReNcell_VM)


# Tests for the SH-SY5Y cell line
# contrast_matrices$SH_SY5Y <- makeContrasts(
#   differences$Adult_Brain$SH_SY5Y,
#   differences$Hippocampus$SH_SY5Y,
#   differences$Insula$SH_SY5Y,
#   levels = design
# )
# 
# fits$SH_SY5Y <- contrasts.fit(fit, contrast_matrices$SH_SY5Y)
# fits$SH_SY5Y <- eBayes(fits$SH_SY5Y, trend = TRUE, robust = TRUE)
# results$SH_SY5Y <- decideTests(fits$SH_SY5Y, method = "separate", adjust.method = "BH", p.value = 0.01, lfc = 2)

# Tests for the SK-N-BE(2) cell line
# contrast_matrices$SK_N_BE_2 <- makeContrasts(
#   differences$Adult_Brain$SK_N_BE_2,
#   differences$Hippocampus$SK_N_BE_2,
#   differences$Insula$SK_N_BE_2,
#   levels = design
# )
# 
# fits$SK_N_BE_2 <- contrasts.fit(fit, contrast_matrices$SK_N_BE_2)
# fits$SK_N_BE_2 <- eBayes(fits$SK_N_BE_2, trend = TRUE, robust = TRUE)
# results$SK_N_BE_2 <- decideTests(fits$SK_N_BE_2, method = "separate", adjust.method = "BH", p.value = 0.01, lfc = 2)

# Tests for the Kelly cell line
# contrast_matrices$Kelly <- makeContrasts(
#   differences$Adult_Brain$Kelly,
#   differences$Hippocampus$Kelly,
#   differences$Insula$Kelly,
#   levels = design
# )
# 
# fits$Kelly <- contrasts.fit(fit, contrast_matrices$Kelly)
# fits$Kelly <- eBayes(fits$Kelly, trend = TRUE, robust = TRUE)
# results$Kelly <- decideTests(fits$Kelly, method = "separate", adjust.method = "BH", p.value = 0.01, lfc = 2)

# Tests for the ReNcell VM cell line
# contrast_matrices$ReNcell_VM <- makeContrasts(
#   differences$Adult_Brain$ReNcell_VM,
#   differences$Hippocampus$ReNcell_VM,
#   differences$Insula$ReNcell_VM,
#   levels = design
# )
# 
# fits$ReNcell_VM <- contrasts.fit(fit, contrast_matrices$ReNcell_VM)
# fits$ReNcell_VM <- eBayes(fits$ReNcell_VM, trend = TRUE, robust = TRUE)
# results$ReNcell_VM <- decideTests(fits$ReNcell_VM, method = "separate", adjust.method = "BH", p.value = 0.01, lfc = 2)



# ############################################
# Replace replicate arrays with the mean value
# ############################################

# eset_ave <- avearrays(eset, ID = eset$targets$CellType, weights = eset$weights)
# eset_ave <- eset_ave[rowSums(is.na(eset_ave$E)) == 0,]
# colnames(eset_ave$E) <- c('Adult Brain', 'Hippocampus', 'Insula', 'Kelly', 'ReNcell VM', 'SH-SY5Y', 'SK-N-BE(2)')



# ################
# Create a heatmap
# ################

# coolmap(eset_ave, cluster.by="de pattern", show.dendrogram = "column", margins = c(7, 1), srtCol=45, labRow='')
# 
# coolmap(eset_ave, cluster.by="expression level", show.dendrogram = "column", col = "redblue", margins = c(7, 1), srtCol=45, labRow='')



# #################
# Create a PCA plot
# #################

# plotMDS(eset_ave, labels = NULL, pch = c(rep(15, 3), rep(16, 4)), col = 1:7)
# legend(-3.5, 4, legend = colnames(eset_ave$E), col = 1:7, pch = c(rep(15, 3), rep(16, 4)), bty = 'n', cex = 0.9)
