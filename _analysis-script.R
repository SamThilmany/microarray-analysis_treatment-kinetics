# ################
# Read in the data
# ################

# Targets
targetinfo <- readTargets(targetsFile, row.names = 'Name')
targetinfo[order(targetinfo$Name),]

# Extract the treatments
treatments <- unique(targetinfo$Treatment)[unique(targetinfo$Treatment) != "NC"]

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
cat(paste0(Sys.time(), ': ', 'The data has been annotated. This step also removed the control probes, resulting in a reduced number of ', dim_eset$annotation[1], ' probes/genes.'), fill = TRUE)



# #############################
# Perform background correction
# #############################

eset <- backgroundCorrect(eset, method = 'normexp')
cat(paste0(Sys.time(), ': Background correction was executed.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/data-processing/density-plot_bg-corrected.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
plotDensities(eset, legend = FALSE, main = 'Density plot after background correction')
dev.off()
cat(paste0(Sys.time(), ': The density plot with the background-corrected data was created'), fill = TRUE)



# ##########################
# Filter out samples with no 
# symbol, and those that are
# not well over the BG
# ##########################

Control <- eset$genes$ControlType != 0
NoSymbol <- is.na(eset$genes$external_gene_name) | eset$genes$external_gene_name == '' | is.na(eset$genes$ensembl_gene_id) |eset$genes$ensembl_gene_id == ''
IsExpr <- rowSums(eset$other$gIsWellAboveBG > 0) == length(colnames(eset))

dim_eset$controlSamples <- dim(eset[!Control, ])
dim_eset$notExpressedSamples <- dim(eset[!IsExpr, ])
dim_eset$samplesWoSymbol <- dim(eset[NoSymbol, ])
dim_eset$samplesToFilter <- dim(eset[Control | NoSymbol | !IsExpr, ])

eset <- eset[!Control & !NoSymbol & IsExpr, ]

dim_eset$filtered <- dim(eset)
cat(paste0(Sys.time(), ': ', 'The data has been filtered. In total, ', dim_eset$samplesToFilter[1], ' samples were omitted (', dim_eset$controlSamples[1], ' control samples, ', dim_eset$notExpressedSamples[1], ' samples that were not significantly above the background, and ', dim_eset$samplesWoSymbol[1], ' samples that have no gene name). The dataset now includes: ', dim_eset$filtered[1], ' genes.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/data-processing/density-plot_filtered.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
plotDensities(eset, legend = FALSE, main = 'Density plot after filtering')
dev.off()
cat(paste0(Sys.time(), ': The density plot of the filtered data was created'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/data-processing/expression-values_filtered.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
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

png(file = paste0(graphicsDirExp, '/data-processing/density-plot_normalized.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
plotDensities(eset, legend = FALSE, main = 'Density plot after normalization')
dev.off()
cat(paste0(Sys.time(), ': The density plot of the normalized data was created.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/data-processing/expression-values_normalized.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
boxplot(log2(eset$E), main = "Expression values after normalization", ylab = "log2 intensity")
dev.off()
cat(paste0(Sys.time(), ': The boxplot of the normalized data was created.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/data-processing/mds.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
plotMDS(eset, labels = substring(eset$targets$Name, 1, nchar(eset$targets$Name) - 3))
dev.off()
cat(paste0(Sys.time(), ': The multidimensional scaling plot visualizing the distance between gene expression profiles of the different samples was created.'), fill = TRUE)



# ############
# Batch effect
# ############

png(file = paste0(graphicsDirExp, '/data-processing/batch-effect.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
plotMDS(eset, labels = substring(eset$targets$Array_Batch, nchar(eset$targets$Array_Batch) - 4 + 1))
dev.off()
cat(paste0(Sys.time(), ': The multidimensional scaling plot visualizing the distance between gene expression profiles of the different array batches was created.'), fill = TRUE)



# ################
# Create a heatmap
# ################

png(file = paste0(graphicsDirExp, '/data-processing/expression-levels.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
coolmap(eset$E, cluster.by="expression level", show.dendrogram = "column", col = "redblue", margins = c(7, 1), srtCol=45, labRow='', main = "Clustered by EL")
dev.off()
cat(paste0(Sys.time(), ': The data has been clustered according to its expression levels.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/data-processing/de-pattern.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
coolmap(eset$E, cluster.by = "de pattern", show.dendrogram = "column", margins = c(7, 1), srtCol = 45, labRow = '', main = "Clustered by DE")
dev.off()
cat(paste0(Sys.time(), ': The data has been clustered according to its DE pattern.'), fill = TRUE)



# #####################
# Array quality weights
# #####################

array.weights <- arrayWeights(eset)
cat(paste0(Sys.time(), ': The array weights have been calculated.'), fill = TRUE)

png(file = paste0(graphicsDirExp, '/data-processing/array-weights.png'), width = graphics_dimensions[1], height = graphics_dimensions[2])
barplot(array.weights, xlab = "Array", ylab = "Weight", main = "Array weights", col = "white", las = 2)
abline(h = 1, lwd = 1, lty = 2)
dev.off()
cat(paste0(Sys.time(), ': The barplot visualizing the array weights has was created.'), fill = TRUE)



# ###################################
# Create a polynomial fit of the data
# ###################################

orthogonalPolynomials <- poly(targetinfo$Time, degree = polyDegree)

model_factors <- list()

for (x in seq(1, polyDegree)) {
  model_factors[[paste0('Baseline.Degree', x)]] <- orthogonalPolynomials[,x]
}

for (tmp_treatment in treatments) {
  for (x in seq(1, polyDegree)) {
    model_factors[[paste0('treat_', tmp_treatment, '.Degree', x)]] <- (targetinfo$Treatment == tmp_treatment & targetinfo$Time > 0) * orthogonalPolynomials[,x]
  }
}

modelMatrixFormula <- 'model.matrix(~'

for (i in seq_along(model_factors)) {
  if (i == 1) {
    modelMatrixFormula <- paste0(modelMatrixFormula, 'model_factors$', names(model_factors[i]))
  } else if (i < length(model_factors)) {
    modelMatrixFormula <- paste0(modelMatrixFormula, '+', 'model_factors$', names(model_factors[i]))
  } else {
    modelMatrixFormula <- paste0(modelMatrixFormula, '+', 'model_factors$', names(model_factors[i]), ')')
  }
}

design <- eval(parse(text = modelMatrixFormula))

fit <- lmFit(eset, design, weights = array.weights)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
fit$coefficients <- as.data.frame(fit$coefficients)

cat('\n')
cat(paste0(Sys.time(), ': The data has been fit to a polynomial regression using the model matrix "', modelMatrixFormula, '".'), fill = TRUE)
cat('\n')



# ###################
# Create volcanoplots
# ###################

for (tmp_treatment in treatments) {
  png(file = paste0(graphicsDirExp, paste0('/volcanoplot_', tmp_treatment, '.png')), width = graphics_dimensions[1], height = graphics_dimensions[2])
  volcanoplot(fit, coef = paste0('model_factors$treat_', tmp_treatment, '.Degree1'), highlight = 10, names = fit$genes$external_gene_name, hl.col = 'red', main = paste0('Volcanoplot for ', tmp_treatment))
  abline(h = -log10(0.05), lwd = 1, lty = 2)
  abline(v = -2, lwd = 1, lty = 2)
  abline(v = 2, lwd = 1, lty = 2)
  dev.off()
  cat(paste0(Sys.time(), ': A volcano plot showing the most significant DEGs of the linear fit for ', tmp_treatment, ' was created.'), fill = TRUE)
}
cat('\n')



# #################
# Create Top Tables
# #################

topTables <- list()

for (tmp_treatment in treatments) {
  tmp_table <- topTable(fit, coef = paste0('model_factors$treat_', tmp_treatment, '.Degree1'), sort.by = 'p', number = Inf)
  cat(paste0(Sys.time(), ': ', 'Minimal p value for ', tmp_treatment, ': ', round(min(tmp_table$P.Value), digits = 4)), fill = TRUE)
  cat(paste0(Sys.time(), ': ', 'Maximal absolute logFC for ', tmp_treatment, ': ', round(max(abs(tmp_table$logFC)), digits = 4)), fill = TRUE)
  
  topTables$noFilter[[paste0('treatment_', tmp_treatment)]] <- tmp_table
}
cat('\n')



# ##################################
# Filter the top tables by the logFC
# ##################################

for (tmp_treatment in treatments) {
  for (tmp_logFC_value in logFC_values) {
    tmp_noFilters_topTable <- topTables$noFilter[[paste0('treatment_', tmp_treatment)]]
    tmp_table <- tmp_noFilters_topTable[abs(tmp_noFilters_topTable$logFC) >= tmp_logFC_value,]
    tmp_table <- tmp_table[order(abs(tmp_table$logFC), decreasing = TRUE),]
    cat(paste0(Sys.time(), ': ', 'After filtering the data of ', tmp_treatment, ' for an absolute logFC value of at least ', tmp_logFC_value, ', ', dim(tmp_table)[1], ' genes were left. The best p value of these remaining genes was ', round(max(tmp_table$P.Value), digits = 4), '.'), fill = TRUE)
    
    topTables[[paste0('logFC_GTE_', tmp_logFC_value)]][[paste0('treatment_', tmp_treatment)]] <- tmp_table
  } 
}
cat('\n')



# ###################################
# Filter the top table by the p value
# ###################################

for (tmp_treatment in treatments) {
  for (tmp_p_value in p_values) {
    tmp_noFilters_topTable <- topTables$noFilter[[paste0('treatment_', tmp_treatment)]]
    tmp_table <- tmp_noFilters_topTable[abs(tmp_noFilters_topTable$P.Value) <= tmp_p_value,]
    tmp_table <- tmp_table[order(abs(tmp_table$P.Value), decreasing = TRUE),]
    cat(paste0(Sys.time(), ': ', 'After filtering the data of ', tmp_treatment, ' for an absolute p value value of at most ', tmp_p_value, ', ', dim(tmp_table)[1], ' genes were left. The best absolute logFC value of these remaining genes was ', round(max(tmp_table$logFC), digits = 4), '.'), fill = TRUE)
    
    topTables[[paste0('pvalue_LTE_', tmp_p_value)]][[paste0('treatment_', tmp_treatment)]] <- tmp_table
  } 
}
cat('\n')



# ##################################################
# Filter the top table by the p value and the log2FC
# ##################################################

for (tmp_treatment in treatments) {
  for (tmp_p_value in p_values) {
    for (tmp_logFC_value in logFC_values) {
      tmp_pValueFiltered_topTable <- topTables[[paste0('pvalue_LTE_', tmp_p_value)]][[paste0('treatment_', tmp_treatment)]]
      tmp_table <- tmp_pValueFiltered_topTable[abs(tmp_pValueFiltered_topTable$logFC) >= tmp_logFC_value,]
      tmp_table <- tmp_table[order(abs(tmp_table$logFC), decreasing = TRUE),]
      cat(paste0(Sys.time(), ': ', 'After filtering the data of ', tmp_treatment, ' for an absolute p value value of at most ', tmp_p_value, ', and a logFC of at least ', tmp_logFC_value, ', ', dim(tmp_table)[1], ' genes were left. The best absolute logFC value of these remaining genes was ', round(max(tmp_table$logFC), digits = 4), '.'), fill = TRUE)
      
      topTables[[paste0('pvalue_LTE_', tmp_p_value, '_logFC_GTE_', tmp_logFC_value)]][[paste0('treatment_', tmp_treatment)]] <- tmp_table
    }
  } 
}
cat('\n')



# ####################
# Create Venn diagrams
# ####################

topTables_group_names <- names(topTables)

for (tmp_group_name in topTables_group_names) {
  if (tmp_group_name == 'noFilter') {
    next
  }
  
  tmp_groupGraphicsDir <- paste0(graphicsDirExp, '/', tmp_group_name)
  
  if (!dir.exists(tmp_groupGraphicsDir)) {
    dir.create(tmp_groupGraphicsDir)
  }
  
  tmp_list <- list()
  for (tmp_list_name in names(topTables[[tmp_group_name]])) {
    tmp_list[[tmp_list_name]] <- topTables[[tmp_group_name]][[tmp_list_name]]$ensembl_gene_id
  }
  
  venn.diagram(
    x = tmp_list,
    category.names = treatments,
    force.unique = TRUE,
    main = paste0('DEG per treatment'),
    fontfamily = "sans",
    main.fontfamily = "sans",
    main.fontface = "bold",
    main.cex = 0.3,
    sub.fontfamily = "sans",
    cat.fontfamily = "sans",
    cex = 0.3,
    cat.cex = 0.3,
    cat.dist = 0.05,
    lwd = 1,
    filename = paste0(tmp_groupGraphicsDir, '/venndiagram.png'),
    output = TRUE,
    imagetype="png",
    width = graphics_dimensions[1],
    height = graphics_dimensions[2],
    resolution = 300,
    margin = 0.1,
    disable.logging = TRUE
  )
}



# ##########################
# Find all overlapping genes
# ##########################

overlapping_genes <- list()
topTables_group_names <- names(topTables)

for (tmp_group_name in topTables_group_names) {
  if (tmp_group_name == 'noFilter') {
    next
  }
  
  
  # Find the overlapping genes between different combinations of treatments
  
  tmp_names <- names(topTables[[tmp_group_name]])
  tmp_combinations <- combn(tmp_names, 2)
  
  for (tmp_combination in seq(1, ncol(tmp_combinations))) {
    overlapping_genes[[tmp_group_name]][[toString(tmp_combinations[,tmp_combination])]] <- intersect(topTables[[tmp_group_name]][[tmp_combinations[,tmp_combination][1]]]$ensembl_gene_id, topTables[[tmp_group_name]][[tmp_combinations[,tmp_combination][2]]]$ensembl_gene_id)
    overlapping_genes[[tmp_group_name]][[toString(tmp_combinations[,tmp_combination])]] <- unique(overlapping_genes[[tmp_group_name]][[toString(tmp_combinations[,tmp_combination])]])

    cat(paste0(Sys.time(), ': ', length(overlapping_genes[[tmp_group_name]][[toString(tmp_combinations[,tmp_combination])]]), ' genes overlap between ', tmp_combinations[,tmp_combination][1], ' and ', tmp_combinations[,tmp_combination][2], '.'), fill = TRUE)
  }
  cat('\n')
  
  
  # Find all the genes that overlap for at least one combination of treatments
  
  overlapping_genes[[tmp_group_name]]$all <- c()

  for (tmp_overlapping_genes_combination in overlapping_genes[[tmp_group_name]]) {
    if (identical(tmp_overlapping_genes_combination, 'all') || identical(tmp_overlapping_genes_combination, 'shared')) {
      next
    }
    
    overlapping_genes[[tmp_group_name]]$all <- c(overlapping_genes[[tmp_group_name]]$all, tmp_overlapping_genes_combination)
  }
  overlapping_genes[[tmp_group_name]]$all <- unique(overlapping_genes[[tmp_group_name]]$all)

  cat(paste0(Sys.time(), ': ', length(overlapping_genes[[tmp_group_name]]$all), ' genes overlap in total for the following filter set: ', tmp_group_name, '.'), fill = TRUE)
  cat('\n')
  
  
  # Find the genes that overlap between all the combinations of treatments
  
  overlapping_genes[[tmp_group_name]]$shared <- c()
  tmp_iterator <- 1
  
  for (tmp_overlapping_genes_combination in names(overlapping_genes[[tmp_group_name]])) {
    if (identical(tmp_overlapping_genes_combination, 'all') || identical(tmp_overlapping_genes_combination, 'shared')) {
      next
    }

    if (tmp_iterator == 1) {
      tmp_expression <- paste0("overlapping_genes[['", tmp_group_name, "']][['", tmp_overlapping_genes_combination, "']]")
    }

    if (tmp_iterator > 1) {
      tmp_expression <- paste0(tmp_expression, " %>% intersect(overlapping_genes[['", tmp_group_name, "']][['", tmp_overlapping_genes_combination, "']])")
    }

    tmp_iterator <- tmp_iterator + 1
  }

  overlapping_genes[[tmp_group_name]]$shared <- eval(parse(text = tmp_expression))
  overlapping_genes[[tmp_group_name]]$shared <- unique(overlapping_genes[[tmp_group_name]]$shared)

  cat(paste0(Sys.time(), ': ', length(overlapping_genes[[tmp_group_name]]$shared), ' genes overlap in all treatments for the following filter set: ', tmp_group_name, '.'), fill = TRUE)
  cat('\n')
}



# ####################################################
# Add a column to the topTables data for unique genes,
# i.e., genes that do not overlap for any treatment
# ####################################################

topTables_group_names <- names(topTables)

for (tmp_group_name in topTables_group_names) {
  if (tmp_group_name == 'noFilter') {
    next
  }
  
  tmp_group_dfs <- topTables[[tmp_group_name]]
  for (tmp_group_df_name in names(tmp_group_dfs)) {
    tmp_df <- tmp_group_dfs[[tmp_group_df_name]]
    tmp_df$unique <- rep(NaN, each = nrow(tmp_df))
    
    for(i in 1:nrow(tmp_df)) {
      tmp_df[i,]$unique <- toString(as.logical(!(tmp_df[i,]$ensembl_gene_id %in% overlapping_genes[[tmp_group_name]]$all)))
    }
    
    topTables[[tmp_group_name]][[tmp_group_df_name]] <- tmp_df
    cat(paste0(Sys.time(), ': ', sum(as.logical(tmp_df$unique)), ' genes were unique for the filter set ', tmp_group_name, ' and the treatment ', tmp_group_df_name, '.'), fill = TRUE)
  }
  cat('\n')
}
cat('\n')



# ############################
# Save all the tables to files
# ############################

# Save the topTables to .csv files

for (tmp_group_name in names(topTables)) {
  tmp_groupResultsDir <- paste0(resultsDirExp, '/', tmp_group_name)
  
  if (!dir.exists(tmp_groupResultsDir)) {
    dir.create(tmp_groupResultsDir)
  }
  
  for (tmp_df_name in names(topTables[[tmp_group_name]])) {
    write.csv(topTables[[tmp_group_name]][[tmp_df_name]], file = paste0(tmp_groupResultsDir, '/topTable_', tmp_df_name, '.csv'), row.names = FALSE)
    cat(paste0(Sys.time(), ': The file ', paste0(tmp_groupResultsDir, '/topTable_', tmp_df_name, '.csv'), ' was created.'), fill = TRUE)
  }
  cat('\n')
}
cat('\n')


# Save the overlapping genes to .csv files

for (tmp_group_name in names(overlapping_genes)) {
  tmp_groupResultsDir <- paste0(resultsDirExp, '/', tmp_group_name)
  
  if (!dir.exists(tmp_groupResultsDir)) {
    dir.create(tmp_groupResultsDir)
  }
  
  for (tmp_list_name in names(overlapping_genes[[tmp_group_name]])) {
    write.csv(overlapping_genes[[tmp_group_name]][[tmp_list_name]], file = paste0(tmp_groupResultsDir, '/overlappingGenes_', tmp_list_name, '.csv'), row.names = FALSE)
    cat(paste0(Sys.time(), ': The file ', paste0(tmp_groupResultsDir, '/overlappingGenes_', tmp_list_name, '.csv'), ' was created.'), fill = TRUE)
  }
  cat('\n')
}
cat('\n')