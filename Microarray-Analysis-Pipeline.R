# Reference: https://support.bioconductor.org/p/9147231/

# #####################
# General configuration
# #####################

# Set the base directory
baseDir <- getwd()

# Define which experiment should be analyzed
experiment = 1

# Create a folder for graphics
graphicsDir <- paste0(baseDir, '/graphics')
graphicsDirExp <- paste0(baseDir, '/graphics/exp', experiment)

if (!dir.exists(graphicsDir)) {
  dir.create(graphicsDir)
}

if (dir.exists(graphicsDirExp)) {
  unlink(graphicsDirExp, recursive = TRUE)
}
dir.create(graphicsDirExp)

# Create a folder for results
resultsDir <- paste0(baseDir, '/results')
resultsDirExp <- paste0(baseDir, '/results/exp', experiment)

if (!dir.exists(resultsDir)) {
  dir.create(resultsDir)
}

if (dir.exists(resultsDirExp)) {
  unlink(resultsDirExp, recursive = TRUE)
}
dir.create(resultsDirExp)

# Targets file
targetsFile <- 'Targets.tsv'

dim_eset <- list()

# Prevent scientific notation
options(scipen = 99)

# Install the Bioconductor Manager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor Packages
BiocManager::install("limma", update = TRUE, ask = FALSE, checkBuilt = TRUE)
require(limma)

# Load R packages
require(statmod)
require(stringr)
require(gplots)
require(tidyr)

# Capture the output to a log file
sink(file = paste0(resultsDirExp, '/log.txt'), append = TRUE, type = c('output', 'message'), split = TRUE)



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

dim_eset$raw <- dim(eset)
cat(paste0(Sys.time(), ': ', 'The raw data has been loaded. The dataset includes ', dim_eset$raw[1], ' probes'), fill = TRUE)

# Add the spot type
spotTypes <- readSpotTypes(file = 'SpotTypes.tsv')
eset$genes$Status <- controlStatus(spotTypes, eset)



# ###################
# Annotate the probes
# ###################

annotation <- read.csv(
  paste0(annotationDir, '/', annotationFile),
  header = TRUE,
  sep = '\t',
  stringsAsFactors = FALSE
)
colnames(annotation)[1] <- 'AgilentID'
annotation <- annotation[!(is.na(annotation$AgilentID) | annotation$AgilentID == ''),]

annotation <- annotation[which(annotation$AgilentID %in% eset$genes$ProbeName),]
annotation <- annotation[match(eset$genes$ProbeName, annotation$AgilentID),]
table(eset$genes$ProbeName == annotation$AgilentID)

eset$genes$AgilentID <- annotation$AgilentID
eset$genes$ensembl_gene_id <- annotation$ensembl_gene_id
eset$genes$external_gene_name <- annotation$external_gene_name

dim_eset$annotation <- dim(eset)
cat(paste0(Sys.time(), ': ', 'The data has been annotated. This step also removed the control probes, resulting in a reduced number of ', dim(dim_eset$annotation)[1], ' probes/genes.'), fill = TRUE)



# #############################
# Perform background correction
# #############################

eset <- backgroundCorrect(eset, method = 'normexp')
cat(paste0(Sys.time(), ': Background correction was executed.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/density-plot_bg-corrected.png'), width = 600, height = 350)
plotDensities(eset, legend = FALSE, main = 'Density plot after background correction')
dev.off()
cat(paste0(Sys.time(), ': The density plot with the background-corrected data was created'), fill = TRUE)



# ##########################
# Filter out samples with no 
# symbol, and those that are
# not well over the BG
# ##########################

Control <- eset$genes$ControlType != 0
NoSymbol <- is.na(eset$genes$external_gene_name) | eset$genes$external_gene_name == ''
IsExpr <- rowSums(eset$other$gIsWellAboveBG > 0) == length(colnames(eset))

dim_eset$controlSamples <- dim(eset[!Control, ])
dim_eset$notExpressedSamples <- dim(eset[!IsExpr, ])
dim_eset$samplesWoSymbol <- dim(eset[NoSymbol, ])
dim_eset$samplesToFilter <- dim(eset[Control | NoSymbol | !IsExpr, ])

eset <- eset[!Control & !NoSymbol & IsExpr, ]

dim_eset$filtered <- dim(eset)
cat(paste0(Sys.time(), ': ', 'The data has been filtered. In total, ', dim_eset$samplesToFilter[1], ' samples were omitted (', dim_eset$controlSamples[1], ' control samples, ', dim_eset$notExpressedSamples[1], ' samples that were not significantly above the background, and ', dim_eset$samplesWoSymbol[1], ' samples that have no gene name). The dataset now includes: ', dim_eset$filtered[1], ' genes.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/density-plot_filtered.png'), width = 600, height = 350)
plotDensities(eset, legend = FALSE, main = 'Density plot after filtering')
dev.off()
cat(paste0(Sys.time(), ': The density plot of the filtered data was created'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/expression-values_filtered.png'), width = 600, height = 350)
boxplot(log2(eset$E), main = "Expression values after filtering", ylab = "log2 intensity")
dev.off()
cat(paste0(Sys.time(), ': The boxplot of the filtered values was created'), fill = TRUE)



# ##########################################
# Remove annotation columns no longer needed
# ##########################################

eset$genes <- eset$genes[, c(
  'ProbeName', 'ensembl_gene_id', 'external_gene_name'
)]



# ##############
# Normalize data
# ##############

eset <- normalizeBetweenArrays(eset, method = 'quantile')
cat(paste0(Sys.time(), ': The expression values were normalized so that they have similar distributions across the arrays.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/density-plot_normalized.png'), width = 600, height = 350)
plotDensities(eset, legend = FALSE, main = 'Density plot after normalization')
dev.off()
cat(paste0(Sys.time(), ': The density plot of the normalized data was created.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/expression-values_normalized.png'), width = 600, height = 350)
boxplot(log2(eset$E), main = "Expression values after normalization", ylab = "log2 intensity")
dev.off()
cat(paste0(Sys.time(), ': The boxplot of the normalized data was created.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/mds.png'), width = 600, height = 350)
plotMDS(eset, labels = substring(eset$targets$Name, 1, nchar(eset$targets$Name) - 3))
dev.off()
cat(paste0(Sys.time(), ': The multidimensional scaling plot visualizing the distance between gene expression profiles of the different samples was created.'), fill = TRUE)



# ############
# Batch effect
# ############

png(file = paste0(graphicsDirExp, '/batch-effect.png'), width = 600, height = 350)
plotMDS(eset, labels = substring(eset$targets$Array_Batch, nchar(eset$targets$Array_Batch) - 4 + 1))
dev.off()
cat(paste0(Sys.time(), ': The multidimensional scaling plot visualizing the distance between gene expression profiles of the different array batches was created.'), fill = TRUE)



# ################
# Create a heatmap
# ################

png(file = paste0(graphicsDirExp, '/expression-levels.png'), width = 600, height = 350)
coolmap(eset$E, cluster.by="expression level", show.dendrogram = "column", col = "redblue", margins = c(7, 1), srtCol=45, labRow='', main = "Clustered by EL")
dev.off()
cat(paste0(Sys.time(), ': The data has been clustered according to its expression levels.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/de-pattern.png'), width = 600, height = 350)
coolmap(eset$E, cluster.by = "de pattern", show.dendrogram = "column", margins = c(7, 1), srtCol = 45, labRow = '', main = "Clustered by DE")
dev.off()
cat(paste0(Sys.time(), ': The data has been clustered according to its DE pattern.'), fill = TRUE)



# #####################
# Array quality weights
# #####################

array.weights <- arrayWeights(eset)
cat(paste0(Sys.time(), ': The array weights have been calculated.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/array-weights.png'), width = 600, height = 350)
barplot(array.weights, xlab = "Array", ylab = "Weight", main = "Array weights", col = "white", las = 2)
abline(h = 1, lwd = 1, lty = 2)
dev.off()
cat(paste0(Sys.time(), ': The barplot visualizing the array weights has was created.'), fill = TRUE)



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
fit$coefficients <- as.data.frame(fit$coefficients)

cat(paste0(Sys.time(), ': The data has been fit to a polynomial regression.'), fill = TRUE)



# ###################
# Create volcanoplots
# ###################

png(file = paste0(graphicsDirExp, '/volcanoplot_EE_linear.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_EE.Lin', highlight = 10, names = fit$genes$external_gene_name, hl.col = 'red', main = 'Volcanoplot for EE (linear regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()
cat(paste0(Sys.time(), ': A volcano plot showing the most significant DEGs of the linear fit for EE was created.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/volcanoplot_LNG_linear.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_LNG.Lin', highlight = 10, names = fit$genes$external_gene_name, main = 'Volcanoplot for LNG (linear regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()
cat(paste0(Sys.time(), ': A volcano plot showing the most significant DEGs of the linear fit for LNG was created.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/volcanoplot_GTx_linear.png'), width = 600, height = 350)
volcanoplot(fit, coef = 'treat_GTx.Lin', highlight = 10, names = fit$genes$external_gene_name, main = 'Volcanoplot for GTx (linear regression)')
abline(h = -log10(0.05), lwd = 1, lty = 2)
abline(v = -2, lwd = 1, lty = 2)
abline(v = 2, lwd = 1, lty = 2)
dev.off()
cat(paste0(Sys.time(), ': A volcano plot showing the most significant DEGs of the linear fit for GTx was created.'), fill = TRUE)



# ############
# Save results
# ############

treatment_EE_linear <- topTable(fit, coef = 'treat_EE.Lin', sort.by = 'p', number = Inf)
cat('\n')
cat(paste0(Sys.time(), ': ', 'Minimal p value for EE: ', round(min(treatment_EE_linear$P.Value), digits = 4)), fill = TRUE)
cat(paste0(Sys.time(), ': ', 'Maximal absolute log2FC for EE: ', round(max(abs(treatment_EE_linear$logFC)), digits = 4)), fill = TRUE)
cat('\n')

treatment_LNG_linear <- topTable(fit, coef = 'treat_LNG.Lin', sort.by = 'p', number = Inf)
cat(paste0(Sys.time(), ': ', 'Minimal p value for LNG: ', round(min(treatment_LNG_linear$P.Value), digits = 4)), fill = TRUE)
cat(paste0(Sys.time(), ': ', 'Maximal absolute log2FC for LNG: ', round(max(abs(treatment_LNG_linear$logFC)), digits = 4)), fill = TRUE)
cat('\n')

treatment_GTx_linear <- topTable(fit, coef = 'treat_GTx.Lin', sort.by = 'p', number = Inf)
cat(paste0(Sys.time(), ': ', 'Minimal p value for GTx: ', round(min(treatment_GTx_linear$P.Value), digits = 4)), fill = TRUE)
cat(paste0(Sys.time(), ': ', 'Maximal absolute log2FC for GTx: ', round(max(abs(treatment_GTx_linear$logFC)), digits = 4)), fill = TRUE)
cat('\n')


treatment_EE_linear_filtered_log2FC <- treatment_EE_linear[abs(treatment_EE_linear$logFC) >= 1,]
treatment_EE_linear_filtered_log2FC <- treatment_EE_linear_filtered_log2FC[order(abs(treatment_EE_linear_filtered_log2FC$logFC), decreasing = TRUE),]
cat(paste0(Sys.time(), ': ', 'After filtering the data of EE for an absolute log2FC value of at least 1, ', dim(treatment_EE_linear_filtered_log2FC)[1], ' genes were left. The best p value of these remaining genes was ', round(max(treatment_EE_linear_filtered_log2FC$P.Value), digits = 4), '.'), fill = TRUE)

treatment_EE_linear_filtered_pValue <- treatment_EE_linear[treatment_EE_linear$P.Value <= 0.05,]
treatment_EE_linear_filtered_pValue <- treatment_EE_linear_filtered_pValue[order(treatment_EE_linear_filtered_pValue$P.Value),]
cat(paste0(Sys.time(), ': ', 'After filtering the data of EE for a p value <= 0.05, ', dim(treatment_EE_linear_filtered_pValue)[1], ' genes were left. The best absolute log2FC of these remaining genes was ', round(max(abs(treatment_EE_linear_filtered_pValue$logFC)), digits = 4), '.'), fill = TRUE)

treatment_EE_linear_filtered_pValue_log2FC <- treatment_EE_linear_filtered_pValue[abs(treatment_EE_linear_filtered_pValue$logFC) >= 1,]
treatment_EE_linear_filtered_pValue_log2FC <- treatment_EE_linear_filtered_pValue_log2FC[order(treatment_EE_linear_filtered_pValue_log2FC$logFC),]
cat(paste0(Sys.time(), ': ', 'After filtering the data of EE for a p value <= 0.05 and an absolute log2FC of at least 1, ', dim(treatment_EE_linear_filtered_pValue_log2FC)[1], ' genes were left.'), fill = TRUE)
cat('\n')


treatment_LNG_linear_filtered_log2FC <- treatment_LNG_linear[abs(treatment_LNG_linear$logFC) >= 1,]
treatment_LNG_linear_filtered_log2FC <- treatment_LNG_linear_filtered_log2FC[order(abs(treatment_LNG_linear_filtered_log2FC$logFC), decreasing = TRUE),]
cat(paste0(Sys.time(), ': ', 'After filtering the data of LNG for an absolute log2FC value of at least 1, ', dim(treatment_LNG_linear_filtered_log2FC)[1], ' genes were left. The best p value of these remaining genes was ', round(max(treatment_LNG_linear_filtered_log2FC$P.Value), digits = 4), '.'), fill = TRUE)

treatment_LNG_linear_filtered_pValue <- treatment_LNG_linear[treatment_LNG_linear$P.Value <= 0.05,]
treatment_LNG_linear_filtered_pValue <- treatment_LNG_linear_filtered_pValue[order(treatment_LNG_linear_filtered_pValue$P.Value),]
cat(paste0(Sys.time(), ': ', 'After filtering the data of LNG for a p value <= 0.05, ', dim(treatment_LNG_linear_filtered_pValue)[1], ' genes were left. The best absolute log2FC of these remaining genes was ', round(max(abs(treatment_LNG_linear_filtered_pValue$logFC)), digits = 4), '.'), fill = TRUE)

treatment_LNG_linear_filtered_pValue_log2FC <- treatment_LNG_linear_filtered_pValue[abs(treatment_LNG_linear_filtered_pValue$logFC) >= 1,]
treatment_LNG_linear_filtered_pValue_log2FC <- treatment_LNG_linear_filtered_pValue_log2FC[order(treatment_LNG_linear_filtered_pValue_log2FC$logFC),]
cat(paste0(Sys.time(), ': ', 'After filtering the data of LNG for a p value <= 0.05 and an absolute log2FC of at least 1, ', dim(treatment_LNG_linear_filtered_pValue_log2FC)[1], ' genes were left.'), fill = TRUE)
cat('\n')


treatment_GTx_linear_filtered_log2FC <- treatment_GTx_linear[abs(treatment_GTx_linear$logFC) >= 1,]
treatment_GTx_linear_filtered_log2FC <- treatment_GTx_linear_filtered_log2FC[order(abs(treatment_GTx_linear_filtered_log2FC$logFC), decreasing = TRUE),]
cat(paste0(Sys.time(), ': ', 'After filtering the data of GTx for an absolute log2FC value of at least 1, ', dim(treatment_GTx_linear_filtered_log2FC)[1], ' genes were left. The best p value of these remaining genes was ', round(max(treatment_GTx_linear_filtered_log2FC$P.Value), digits = 4), '.'), fill = TRUE)

treatment_GTx_linear_filtered_pValue <- treatment_GTx_linear[treatment_GTx_linear$P.Value <= 0.05,]
treatment_GTx_linear_filtered_pValue <- treatment_GTx_linear_filtered_pValue[order(treatment_GTx_linear_filtered_pValue$P.Value),]
cat(paste0(Sys.time(), ': ', 'After filtering the data of GTx for a p value <= 0.05, ', dim(treatment_GTx_linear_filtered_pValue)[1], ' genes were left. The best absolute log2FC of these remaining genes was ', round(max(abs(treatment_GTx_linear_filtered_pValue$logFC)), digits = 4), '.'), fill = TRUE)

treatment_GTx_linear_filtered_pValue_log2FC <- treatment_GTx_linear_filtered_pValue[abs(treatment_GTx_linear_filtered_pValue$logFC) >= 1,]
treatment_GTx_linear_filtered_pValue_log2FC <- treatment_GTx_linear_filtered_pValue_log2FC[order(treatment_GTx_linear_filtered_pValue_log2FC$logFC),]
cat(paste0(Sys.time(), ': ', 'After filtering the data of GTx for a p value <= 0.05 and an absolute log2FC of at least 1, ', dim(treatment_GTx_linear_filtered_pValue_log2FC)[1], ' genes were left.'), fill = TRUE)
cat('\n')


setwd(resultsDirExp)

tryCatch(
  expr = {
    write.csv(treatment_EE_linear_filtered_log2FC, file = paste0('treatment_EE_linear_filtered_log2FC.csv'), row.names = FALSE)
    write.csv(treatment_EE_linear_filtered_pValue, file = paste0('treatment_EE_linear_filtered_pValue.csv'), row.names = FALSE)
    write.csv(treatment_EE_linear_filtered_pValue_log2FC, file = paste0('treatment_EE_linear_filtered_pValue_log2FC.csv'), row.names = FALSE)
    
    write.csv(treatment_LNG_linear_filtered_log2FC, file = paste0('treatment_LNG_linear_filtered_log2FC.csv'), row.names = FALSE)
    write.csv(treatment_LNG_linear_filtered_pValue, file = paste0('treatment_LNG_linear_filtered_pValue.csv'), row.names = FALSE)
    write.csv(treatment_LNG_linear_filtered_pValue_log2FC, file = paste0('treatment_LNG_linear_filtered_pValue_log2FC.csv'), row.names = FALSE)
    
    write.csv(treatment_GTx_linear_filtered_log2FC, file = paste0('treatment_GTx_linear_filtered_log2FC.csv'), row.names = FALSE)
    write.csv(treatment_GTx_linear_filtered_pValue, file = paste0('treatment_GTx_linear_filtered_pValue.csv'), row.names = FALSE)
    write.csv(treatment_GTx_linear_filtered_pValue_log2FC, file = paste0('treatment_GTx_linear_filtered_pValue_log2FC.csv'), row.names = FALSE)
  },
  finally = {
    setwd(baseDir)
  }
)

cat(paste0(Sys.time(), ': The results have been saved to `.csv` files. The pipeline has come to a successful end.'), fill = TRUE)

# End the output capture
sink(file = NULL, split = FALSE)
