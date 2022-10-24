# Reference: https://support.bioconductor.org/p/9147231/

# #####################
# General configuration
# #####################

# Set the base directory
baseDir <- getwd()

# Define which experiment should be analyzed
experiment = 2

# Create a folder for graphics
graphicsDir <- paste0(baseDir, '/graphics')
graphicsDirExp <- paste0(baseDir, '/graphics/exp', experiment)

if (!dir.exists(graphicsDir)) {
  dir.create(graphicsDir)
}

if (!dir.exists(graphicsDirExp)) {
  dir.create(graphicsDirExp)
}

# Create a folder for results
resultsDir <- paste0(baseDir, '/results')
resultsDirExp <- paste0(baseDir, '/results/exp', experiment)

if (!dir.exists(resultsDir)) {
  dir.create(resultsDir)
}

if (!dir.exists(resultsDirExp)) {
  dir.create(resultsDirExp)
}

# Targets file
targetsFile <- 'Targets.tsv'

dim_eset <- list()

# Prevent scientific notation
options(scipen = 99)

# Install the Bioconductor Manager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor Packages
BiocManager::install(c("limma"))
library(limma)

# Load R packages
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

dim_eset[['raw']] <- dim(eset)
cat(paste0('The raw data has been loaded. The dataset includes ', dim_eset$raw[1], ' genes.'), fill = TRUE)

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

cat('The data has been annotated.', fill = TRUE)



# #############################
# Perform background correction
# #############################

eset <- backgroundCorrect(eset, method = 'normexp')

png(file = paste0(graphicsDirExp, '/density-plot_bg-corrected.png'), width = 600, height = 350)
plotDensities(eset, legend = FALSE, main = 'Density plot after background correction')
dev.off()

cat('Background correction was executed.', fill = TRUE)



# ##########################
# Filter out control probes, 
# those with no symbol, 
# and those that fail
# ##########################

Control <- eset$genes$ControlType != 0
NoSymbol <- is.na(eset$genes$external_gene_name) | eset$genes$external_gene_name == ''
IsExpr <- rowSums(eset$other$gIsWellAboveBG > 0) == length(colnames(eset))

dim_eset[['controlSamples']] <- dim(eset[Control, ])
dim_eset[['notExpressedSamples']] <- dim(eset[!IsExpr, ])
dim_eset[['samplesWoSymbol']] <- dim(eset[NoSymbol, ])
dim_eset[['samplesToFilter']] <- dim(eset[Control | NoSymbol | !IsExpr, ])

eset <- eset[!Control & !NoSymbol & IsExpr, ]

dim_eset[['filtered']] <- dim(eset)
cat(paste0('The data has been filtered. In total, ', dim_eset$samplesToFilter[1], ' samples were omitted (', dim_eset$controlSamples[1], ' control samples, ', dim_eset$notExpressedSamples[1], ' samples that were not significantly above the background, and ', dim_eset$samplesWoSymbol[1], ' samples that have no gene name). The dataset now includes: ', dim_eset$filtered[1], ' genes.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/density-plot_filtered.png'), width = 600, height = 350)
plotDensities(eset, legend = FALSE, main = 'Density plot after filtering')
dev.off()

png(file = paste0(graphicsDirExp, '/expression-values_filtered.png'), width = 600, height = 350)
boxplot(log2(eset$E), main = "Expression values after filtering", ylab = "log2 intensity")
dev.off()



# ##########################################
# Remove annotation columns no longer needed
# ##########################################

eset$genes <- eset$genes[, c(
  'ProbeName', 'AgilentID', 'SystematicName', 'Status'
)]



# ##############
# Normalize data
# ##############

eset <- normalizeBetweenArrays(eset, method = 'quantile')

png(file = paste0(graphicsDirExp, '/density-plot_normalized.png'), width = 600, height = 350)
plotDensities(eset, legend = FALSE, main = 'Density plot after normalization')
dev.off()

png(file = paste0(graphicsDirExp, '/expression-values_normalized.png'), width = 600, height = 350)
boxplot(log2(eset$E), main = "Expression values after normalization", ylab = "log2 intensity")
dev.off()

cat('The data has been normalized.', fill = TRUE)

png(file = paste0(graphicsDirExp, '/mds.png'), width = 600, height = 350)
plotMDS(eset, labels = substring(eset$targets$Name, 1, nchar(eset$targets$Name) - 3))
dev.off()



# ############
# Batch effect
# ############

png(file = paste0(graphicsDirExp, '/batch-effect.png'), width = 600, height = 350)
plotMDS(eset, labels = substring(eset$targets$Array_Batch, nchar(eset$targets$Array_Batch) - 4 + 1))
dev.off()

cat('An MDS plot for the visualization of a batch effect has been generated.', fill = TRUE)



# ################
# Create a heatmap
# ################

png(file = paste0(graphicsDirExp, '/expression-levels.png'), width = 600, height = 350)
coolmap(eset$E, cluster.by="expression level", show.dendrogram = "column", col = "redblue", margins = c(7, 1), srtCol=45, labRow='', main = "Clustered by EL")
dev.off()

cat('The data has been clustered according to its expression levels.', fill = TRUE)

png(file = paste0(graphicsDirExp, '/de-pattern.png'), width = 600, height = 350)
coolmap(eset$E, cluster.by = "de pattern", show.dendrogram = "column", margins = c(7, 1), srtCol = 45, labRow = '', main = "Clustered by DE")
dev.off()

cat('The data has been clustered according to its DE pattern.', fill = TRUE)



# #####################
# Array quality weights
# #####################

array.weights <- arrayWeights(eset)

png(file = paste0(graphicsDirExp, '/array-weights.png'), width = 600, height = 350)
barplot(array.weights, xlab = "Array", ylab = "Weight", main = "Array weights", col = "white", las = 2)
abline(h = 1, lwd = 1, lty = 2)
dev.off()

cat('The array weights have been calculated.', fill = TRUE)



# ###############################
# Create a linear fit of the data
# ###############################

X <- poly(targetinfo$Time, degree = 3)
Baseline.Lin  <- X[,1]
Baseline.Quad <- X[,2]
Baseline.Cubic <- X[,3]

treat_EE.Lin  <- (targetinfo$Treatment=="EE" & targetinfo$Time > 0) * X[,1]
treat_EE.Quad <- (targetinfo$Treatment=="EE" & targetinfo$Time > 0) * X[,2]
treat_EE.Cubic <- (targetinfo$Treatment=="EE" & targetinfo$Time > 0) * X[,3]

treat_LNG.Lin  <- (targetinfo$Treatment=="LNG" & targetinfo$Time > 0) * X[,1]
treat_LNG.Quad <- (targetinfo$Treatment=="LNG" & targetinfo$Time > 0) * X[,2]
treat_LNG.Cubic <- (targetinfo$Treatment=="LNG" & targetinfo$Time > 0) * X[,3]

treat_GTx.Lin  <- (targetinfo$Treatment=="GTx" & targetinfo$Time > 0) * X[,1]
treat_GTx.Quad <- (targetinfo$Treatment=="GTx" & targetinfo$Time > 0) * X[,2]
treat_GTx.Cubic <- (targetinfo$Treatment=="GTx" & targetinfo$Time > 0) * X[,3]

design <- model.matrix( ~ Baseline.Lin + Baseline.Quad + Baseline.Cubic
                        + treat_EE.Lin + treat_EE.Quad + treat_EE.Cubic
                        + treat_LNG.Lin + treat_LNG.Quad + treat_LNG.Cubic
                        + treat_GTx.Lin + treat_GTx.Quad + treat_GTx.Cubic )

fit <- lmFit(eset, design, weights = array.weights)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

cat('The data has been fit to a polynomial regression.', fill = TRUE)



# ###################
# Create volcanoplots
# ###################

png(file = paste0(graphicsDirExp, '/volcanoplot_EE_linear.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_EE.Lin', highlight = 10, names = fit$genes$external_gene_name, hl.col = 'red', main = 'Volcanoplot for EE (linear regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()

png(file = paste0(graphicsDirExp, '/volcanoplot_EE_quadratic.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_EE.Quad', highlight = 10, names = fit$genes$external_gene_name, hl.col = 'red', main = 'Volcanoplot for EE (quadratic regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()

png(file = paste0(graphicsDirExp, '/volcanoplot_EE_cubic.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_EE.Cubic', highlight = 10, names = fit$genes$external_gene_name, hl.col = 'red', main = 'Volcanoplot for EE (cubic regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()


png(file = paste0(graphicsDirExp, '/volcanoplot_LNG_linear.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_LNG.Lin', highlight = 10, names = fit$genes$external_gene_name, main = 'Volcanoplot for LNG (linear regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()

png(file = paste0(graphicsDirExp, '/volcanoplot_LNG_quadratic.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_LNG.Quad', highlight = 10, names = fit$genes$external_gene_name, main = 'Volcanoplot for LNG (quadratic regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()

png(file = paste0(graphicsDirExp, '/volcanoplot_LNG_cubic.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_LNG.Cubic', highlight = 10, names = fit$genes$external_gene_name, main = 'Volcanoplot for LNG (cubic regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()


png(file = paste0(graphicsDirExp, '/volcanoplot_GTx_linear.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_GTx.Lin', highlight = 10, names = fit$genes$external_gene_name, main = 'Volcanoplot for GTx (linear regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()

png(file = paste0(graphicsDirExp, '/volcanoplot_GTx_quadratic.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_GTx.Quad', highlight = 10, names = fit$genes$external_gene_name, main = 'Volcanoplot for GTx (quadratic regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()

png(file = paste0(graphicsDirExp, '/volcanoplot_GTx_cubic.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_GTx.Cubic', highlight = 10, names = fit$genes$external_gene_name, main = 'Volcanoplot for GTx (cubic regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()



# ############
# Save results
# ############

treatment_EE_top10_linear <- topTable(fit, sort.by = 'p', coef = 'treat_EE.Lin')
treatment_EE_top10_quadratic <- topTable(fit, sort.by = 'p', 'treat_EE.Quad')
treatment_EE_top10_cubic <- topTable(fit, sort.by = 'p', coef = 'treat_EE.Cubic')

treatment_LNG_top10_linear <- topTable(fit, sort.by = 'p', coef = 'treat_LNG.Lin')
treatment_LNG_top10_quadratic <- topTable(fit, sort.by = 'p', coef = 'treat_LNG.Quad')
treatment_LNG_top10_cubic <- topTable(fit, sort.by = 'p', coef = 'treat_LNG.Cubic')

treatment_GTx_top10_linear <- topTable(fit, sort.by = 'p', coef = 'treat_GTx.Lin')
treatment_GTx_top10_quadratic <- topTable(fit, sort.by = 'p', coef = 'treat_GTx.Quad')
treatment_GTx_top10_cubic <- topTable(fit, sort.by = 'p', coef = 'treat_GTx.Cubic')


treatment_EE_linear <- topTable(fit, sort.by = 'p', coef = 'treat_EE.Lin', number = Inf)
treatment_EE_quadratic <- topTable(fit, sort.by = 'p', coef = 'treat_EE.Quad', number = Inf)
treatment_EE_cubic <- topTable(fit, sort.by = 'p', coef = 'treat_EE.Cubic', number = Inf)

treatment_LNG_linear <- topTable(fit, sort.by = 'p', coef = 'treat_LNG.Lin', number = Inf)
treatment_LNG_quadratic <- topTable(fit, sort.by = 'p', coef = 'treat_LNG.Quad', number = Inf)
treatment_LNG_cubic <- topTable(fit, sort.by = 'p', coef = 'treat_LNG.Cubic', number = Inf)

treatment_GTx_linear <- topTable(fit, sort.by = 'p', coef = 'treat_GTx.Lin', number = Inf)
treatment_GTx_quadratic <- topTable(fit, sort.by = 'p', coef = 'treat_GTx.Quad', number = Inf)
treatment_GTx_cubic <- topTable(fit, sort.by = 'p', coef = 'treat_GTx.Cubic', number = Inf)

setwd(resultsDirExp)

tryCatch(
  expr = {
    write.csv(treatment_EE_top10_linear, file = paste0('treatment_EE_top10_linear.csv'), row.names = FALSE)
    write.csv(treatment_EE_top10_quadratic, file = paste0('treatment_EE_top10_quadratic.csv'), row.names = FALSE)
    write.csv(treatment_EE_top10_cubic, file = paste0('treatment_EE_top10_cubic.csv'), row.names = FALSE)
    
    write.csv(treatment_LNG_top10_linear, file = paste0('treatment_LNG_top10_linear.csv'), row.names = FALSE)
    write.csv(treatment_LNG_top10_quadratic, file = paste0('treatment_LNG_top10_quadratic.csv'), row.names = FALSE)
    write.csv(treatment_LNG_top10_cubic, file = paste0('treatment_LNG_top10_cubic.csv'), row.names = FALSE)
    
    write.csv(treatment_GTx_top10_linear, file = paste0('treatment_GTx_top10_linear.csv'), row.names = FALSE)
    write.csv(treatment_GTx_top10_quadratic, file = paste0('treatment_GTx_top10_quadratic.csv'), row.names = FALSE)
    write.csv(treatment_GTx_top10_cubic, file = paste0('treatment_GTx_top10_cubic.csv'), row.names = FALSE)
    
    
    write.csv(treatment_EE_linear, file = paste0('treatment_EE_linear.csv'), row.names = FALSE)
    write.csv(treatment_EE_quadratic, file = paste0('treatment_EE_quadratic.csv'), row.names = FALSE)
    write.csv(treatment_EE_cubic, file = paste0('treatment_EE_cubic.csv'), row.names = FALSE)
    
    write.csv(treatment_LNG_linear, file = paste0('treatment_LNG_linear.csv'), row.names = FALSE)
    write.csv(treatment_LNG_quadratic, file = paste0('treatment_LNG_quadratic.csv'), row.names = FALSE)
    write.csv(treatment_LNG_cubic, file = paste0('treatment_LNG_cubic.csv'), row.names = FALSE)
    
    write.csv(treatment_GTx_linear, file = paste0('treatment_GTx_linear.csv'), row.names = FALSE)
    write.csv(treatment_GTx_quadratic, file = paste0('treatment_GTx_quadratic.csv'), row.names = FALSE)
    write.csv(treatment_GTx_cubic, file = paste0('treatment_GTx_cubic.csv'), row.names = FALSE)
  },
  finally = {
    setwd(baseDir)
  }
)

cat('The results have been saved to `.csv` files. The pipeline has come to a successful end.', fill = TRUE)
