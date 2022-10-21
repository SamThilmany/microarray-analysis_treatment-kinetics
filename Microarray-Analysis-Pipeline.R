# #####################
# General configuration
# #####################

baseDir <- getwd()
print(baseDir)

graphicsDir <- paste0(baseDir, "/graphics")
resultsDir <- paste0(baseDir, "/results")

treatment = "GTx"
experiment = 1

targetsFile <- 'Targets.tsv'

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
# Read in the data
# ################

# Targets
targetinfo <- readTargets(targetsFile, row.names = 'Name')
targetinfo[order(targetinfo$Name),]

# Select a subset for a specific experiment
targetinfo <- subset(targetinfo, Experiment == experiment)
targetinfo <- subset(targetinfo, Treatment == treatment | Treatment == "NC")

# Converts the raw data to an EListRaw object
wtAgilent.GFilter <- function(qta) { qta[,"gIsPosAndSignif"] }
eset <- read.maimages(
  targetinfo,
  source = 'agilent.median',
  green.only = TRUE,
  path = "data",
  names = targetinfo$Name,
  other.columns = 'gIsWellAboveBG',
  wt.fun = wtAgilent.GFilter
)

# Add the spot type
spotTypes <- readSpotTypes(file = 'SpotTypes.tsv')
eset$genes$Status <- controlStatus(spotTypes, eset)



# ###################
# Annotate the probes
# ###################

annotLookup <- read.csv(
  paste0(annotationDir, '/', annotationFile),
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

png(file = paste0(graphicsDir, '/', treatment, '_density-plot_bg-corrected.png'), width = 600, height = 350)
plotDensities(eset, legend = FALSE, main = 'Density plot after background correction')
dev.off()



# ##############
# Normalize data
# ##############

eset <- normalizeBetweenArrays(eset, method = 'quantile')

png(file = paste0(graphicsDir, '/', treatment, '_density-plot_normalized.png'), width = 600, height = 350)
plotDensities(eset, legend = FALSE, main = 'Density plot after normalization')
dev.off()



# ##########################
# Filter out control probes, 
# those with no symbol, 
# and those that fail
# ##########################

Control <- eset$genes$ControlType != 0
NoSymbol <- is.na(eset$genes$external_gene_name)
IsExpr <- rowSums(eset$other$gIsWellAboveBG > 0) == length(colnames(eset))

eset <- eset[!Control & !NoSymbol & IsExpr, ]



# ##########################################
# Remove annotation columns no longer needed
# ##########################################

eset$genes <- eset$genes[,c(
  'ProbeName','wikigene_description','ensembl_gene_id','entrezgene','gene_biotype','external_gene_name'
)]



# #################
# Expression values
# #################

png(file = paste0(graphicsDir, '/', treatment, '_expression-values.png'), width = 600, height = 350)
boxplot(log2(eset$E), main = "Expression values", ylab = "log2 intensity")
dev.off()



# ############
# Batch effect
# ############

png(file = paste0(graphicsDir, '/', treatment, '_batch-effect.png'), width = 600, height = 350)
plotMDS(eset, labels = eset$targets$Array_Batch)
dev.off()



# ################
# Create a heatmap
# ################

png(file = paste0(graphicsDir, '/', treatment, '_de-pattern.png'), width = 600, height = 350)
coolmap(eset$E, cluster.by = "de pattern", show.dendrogram = "column", margins = c(7, 1), srtCol = 45, labRow = '', main = "Clustered by DE")
dev.off()

png(file = paste0(graphicsDir, '/', treatment, '_expression-levels.png'), width = 600, height = 350)
coolmap(eset$E, cluster.by="expression level", show.dendrogram = "column", col = "redblue", margins = c(7, 1), srtCol=45, labRow='', main = "Clustered by EL")
dev.off()



# #####################
# Array quality weights
# #####################

array.weights <- arrayWeights(eset)

png(file = paste0(graphicsDir, '/', treatment, '_array-weights.png'), width = 600, height = 350)
barplot(array.weights, xlab = "Array", ylab = "Weight", main = "Array weights", col = "white", las = 2)
abline(h = 1, lwd = 1, lty = 2)
dev.off()



# ###############################
# Create a linear fit of the data
# ###############################

time_factor <- factor(targetinfo$Time)

design <- model.matrix(~ 0 + time_factor)

fit <- lmFit(eset, design = design, weights = array.weights)



# ##############
# Make contrasts
# ##############

contrasts <- makeContrasts(
  time_factor2-time_factor0,
  time_factor5-time_factor0,
  time_factor7-time_factor0,
  time_factor10-time_factor0,
  levels = design
)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)



# #############
# Venn diagrams
# #############

png(file = paste0(graphicsDir, '/', treatment, '_vennDiagram.png'), width = 600, height = 350)
vennDiagram(fit, include = c("up", "down"), names = c("D2-D0", "D5-D0", "D7-D0", "D10-D0"), counts.col = c("red", "blue"), main = "Up- and down-regulated genes")
dev.off()



# ############
# Save results
# ############

# Create a folder for result files
if (!dir.exists(resultsDir)) {
  dir.create(resultsDir)
}

treatment_deg_D2 <- topTable(fit, coef = 1, n = Inf)
treatment_deg_D5 <- topTable(fit, coef = 2, n = Inf)
treatment_deg_D7 <- topTable(fit, coef = 3, n = Inf)
treatment_deg_D10 <- topTable(fit, coef = 4, n = Inf)

treatment_deg_D2_top20 <- topTable(fit, coef = 1, n = 20)
treatment_deg_D5_top20 <- topTable(fit, coef = 2, n = 20)
treatment_deg_D7_top20 <- topTable(fit, coef = 3, n = 20)
treatment_deg_D10_top20 <- topTable(fit, coef = 4, n = 20)

setwd(resultsDir)

tryCatch(
  expr = {
    write.csv(treatment_deg_D2, file = paste0("treatment", treatment, "_deg_D2.csv"), row.names = FALSE)
    write.csv(treatment_deg_D5, file = paste0("treatment", treatment, "_deg_D5.csv"), row.names = FALSE)
    write.csv(treatment_deg_D7, file = paste0("treatment", treatment, "_deg_D7.csv"), row.names = FALSE)
    write.csv(treatment_deg_D10, file = paste0("treatment", treatment, "_deg_D10.csv"), row.names = FALSE)
    
    write.csv(treatment_deg_D2_top20, file = paste0("treatment", treatment, "_deg_D2_top20.csv"), row.names = FALSE)
    write.csv(treatment_deg_D5_top20, file = paste0("treatment", treatment, "_deg_D5_top20.csv"), row.names = FALSE)
    write.csv(treatment_deg_D7_top20, file = paste0("treatment", treatment, "_deg_D7_top20.csv"), row.names = FALSE)
    write.csv(treatment_deg_D10_top20, file = paste0("treatment", treatment, "_deg_D10_top20.csv"), row.names = FALSE)
  },
  finally = {
    setwd(baseDir)
  }
)

png(file = paste0(graphicsDir, '/', treatment, '_volcanoplot.png'), width = 600, height = 350)
volcanoplot(fit, coef = 1)
dev.off()

