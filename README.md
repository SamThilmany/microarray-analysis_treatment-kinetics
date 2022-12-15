# Microarray Analysis of Treatment Kinetics

This R script is used to analyze microarray data acquired by an Agilent SureScan Microarray Scanner. Internally the package `limma` by Gordon K. Smyth *et al.* [[1]](https://doi.org/10.18129/B9.bioc.limma "Limma Package") is used to read and analyze the data.

## Aim of the Experiment

This script is used to analyze the raw readouts of two biological experiments. In both experiments, the cell culture media were supplemented with three different pharmaceuticals, resulting in four samples per time point; three treatments, and one negative control. The data should show whether treatment with different pharmaceuticals affects gene expression and, if so, how this influence changes over time. Furthermore, we are interested in the specific genes that changed over time and how these genes relate to each treatment, *i.e.*, whether a particular gene is differentially expressed for only one treatment or multiple treatments and whether a gene that is differentially expressed for multiple treatments is always up- or down-regulated.

## Raw Data Generation

Cells were harvested at different time points, and RNA was isolated using a modified protocol of Qiagen's RNeasy Protect Cell Mini Kit Protocol.[[2]](https://www.qiagen.com/us/resources/resourcedetail?id=bea04757-b25e-4eb1-86b4-3ef1cb4f94b0&lang=en "Qiagen Protocol") Further sample preparation and microarray-based gene expression analysis was performed according to Agilent's protocol.[[3]](https://www.agilent.com/cs/library/usermanuals/Public/G4140-90040_GeneExpression_OneColor_6.9.pdf "Microarray Protocol") The microarray used was a SurePrint G3 Human Gene Expression v3 8x60K (P/N G4851C). Each glass slide carries eight high-definition 60K arrays containing cDNA for 26,803 unique Entrez genes and 30,606 unique lncRNAs, and 3000 replicates.[[4]](https://agilent.com/store/de_DE/Prod-G4851C/G4851C "SurePrint G3 Human Gene Expression v3 8x60K Microarray Kit")

## Raw Data Analysis

The Agilent SureScan generates a QC Report and a `.txt` file with various data and metadata for each sample. Eight samples could be analyzed per glass slide; the files of these eight samples were stored in a folder named after the serial number of the glass slide. (See `data/*`)

After importing the raw data via `limma`'s `read.maimages` function [[5]](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/read.maimages "read.maimages"), the data was annotated with the ensembl gene IDs. Background correction was done via the `backgroundCorrect` function from `limma` using the `normexp` method.[[6]](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/backgroundCorrect "backgroundCorrect") After background correction, the data were filtered; here, the data that could not be annotated with data from ensembl and whose expression level was not significantly above the background were excluded. The remaining data were normalized with the `limma` function `normalizeBetweenArrays` using the `quantile` method to ensure a similar distribution of expression levels between the different arrays.[[7]](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/normalizeBetweenArrays "normalizeBetweenArrays")

A polynomial trend was allowed for the baseline and the individual treatments to examine the effect of each treatment over the treatment period. These polynomial trends were used to create a design matrix fitted to the data using the `lmFit` function from `limma`.[[9]](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/lmFit "lmFit") Professor Gordon K. Smyth of the Walter and Eliza Hall Institute of Medical Research in Melbourne, creator of the `limma` package, advised this procedure.[[8]](https://support.bioconductor.org/p/9147231/ "Advice from Professor Smyth")

The fitted data were then statistically analyzed using empirical Bayes statistics for differential expression. The `eBayes` function from `limma` was used for this purpose.[[10]](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/ebayes "eBayes")

The generated data were filtered according to different criteria (logFC, p-value, and p-value with logFC). For all criteria, lists of genes that were differentially expressed in multiple or all treatments were generated; genes that were differentially expressed in only one treatment were also flagged.

### Want to know more?

For more precise information on the evaluation, you may look at the code and refer to the individual functions in the manual of `limma` [[1]](https://doi.org/10.18129/B9.bioc.limma "Limma Package") or the respective R package.

## About this Project

This project is part of Teresa Hardy's Master's thesis, which was conducted under the supervision of Sam Thilmany at the Federal Institute for Drugs and Medical Devices, Bonn, Germany.

The concept for the data analysis was a joint effort of Teresa Hardy and Sam Thilmany; Sam Thilmany did the programming.
