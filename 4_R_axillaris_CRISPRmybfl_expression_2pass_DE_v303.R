rm(list=ls())
library(BiocManager)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library("VennDiagram")
library("plotrix")
library("Rmisc")
library("flashClust")
library("MBCluster.Seq")
library("dplyr")
library("sfsmisc")
library("DESeq2")

# The following code was adapted from Stephen Turner's pipeline https://gist.github.com/stephenturner/f60c1934405c127f09a6

#### samples are from stage 4 limb tissue, two independent lines of mutant/axN control ####
#Name	ID	seqname
#KML 3-11-1 axN	L3_axN_1	RNASeq1
#KML 3-12-1 axN	L3_axN_2	RNASeq2
#KML 3-15-3 axN	L3_axN_3	RNASeq3
#KML 3-17-1 axN	L3_axN_4	RNASeq4
#KML 3-6-1 mutant	L3_mut_1	RNASeq5
#KML 3-9-1 mutant	L3_mut_2	RNASeq6
# KML 3-20-2 mutant	L3_mut_3	RNASeq7
# KML 3-24-1 mutant	L3_mut_4	RNASeq8
# KML 10-3-1 axN	L10_axN_1	RNASeq9
# KML 10-6-2 axN	L10_axN_2	RNASeq10
# KML 10-8-1 axN	L10_axN_3	RNASeq11
# KML 10-12-1 axN	L10_axN_4	RNASeq12
# KML 10-2-1 mutant	L10_mut_1	RNASeq13
# KML 10-9-1 mutant	L10_mut_2	RNASeq14
# KML 10-16-1 mutant	L10_mut_3	RNASeq15
# KML 10-20-1 mutant	L10_mut_4	RNASeq16


setwd("/path to folder")

# Import & pre-process ----------------------------------------------------

# Import data from featureCounts
## Previously ran at command line
# Import countdata from featureCounts
countdata <- read.table("countdata_file.txt", header=TRUE, skip=1, row.names=1)
names(countdata)

#store this data for later merges
countdata.info <- countdata[ ,1:5]

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
# Remove other stuff from filenames
colnames(countdata) <- gsub("\\_Aligned.sortedByCoord.out$", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Make sure the order of samples is correct
head(countdata)

# Assign condition (reps of axN and mutant, could also try just axN and mutant to alter the DE analysis depending on your question)
(condition <- factor(c(rep("axN", 4), rep("crispr_mut", 4), rep("axN", 4), rep("crispr_mut", 4))))

# Analysis with DESeq2 ----------------------------------------------------

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline, which estimates size factors of each sample, estimates dispersion, and negative binomial GLM fitting and Wald statistics
dds <- DESeq(dds)

## which transformation method is best? #https://www.bioconductor.org/help/workflows/rnaseqGene/#the-deseqdataset-object-sample-information-and-the-design-formula
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))
#VST
vsd <- vst(dds)
head(assay(vsd))
hist(assay(vsd))

#### Sample distances ####
#A useful first step in an RNA-seq analysis is often to assess overall similarity between samples: Which samples are similar to each other, which are different? Does this fit to the expectation from the experiment’s design?
#We use the R function dist to calculate the Euclidean distance between samples. To ensure we have a roughly equal contribution from all genes, we use it on the rlog-transformed data. We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.
sampleDists <- dist(t(assay(rld)))
head(sampleDists)

# Principal components analysis
## Could do with built-in DESeq2 function:
DESeq2::plotPCA(rld, intgroup="condition")

#### DE expression analysis ####
dds.st4 <- DESeq(dds)

#overall model with LRT method #https://www.biostars.org/p/115685/
#does gene expression change due to species condition, in general?
#ddsLRT.st4 <- DESeq(dds.st4, test="LRT", reduced= ~ 1) #somehow the general 
#resLRT.st4 <- results(ddsLRT.st4, tidy = TRUE) #tidy puts it into a DF automatically

#contrasts, add the alpha up here, also use tidy to make DF
res.st4.ax.mut <- results(dds.st4, contrast=c("condition","axN","crispr_mut"), alpha = 0.01) #DESeqResults object, retains your comparison information
  res.st4.ax.mut <- results(dds.st4, contrast=c("condition","axN","crispr_mut"), alpha = 0.01, tidy = TRUE)
  res.st4.ax.mut <- na.omit(res.st4.ax.mut)
