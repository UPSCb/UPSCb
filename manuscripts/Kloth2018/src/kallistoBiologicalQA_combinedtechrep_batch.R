#' ---
#' title: "Arabidopsis aphid project - biological QA - combined technical replicates including batch factor"
#' author: "Nicolas Delhomme & Karen Kloth"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/arabidopsis/balbrectsen/arabidopsis-aphids")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/arabidopsis/balbrectsen/arabidopsis-aphids")
#' ```
#' 

#' Librarius
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))

#' Define a palette
pal <- brewer.pal(8,"Dark2")

#' Save the default margin parameters
mar=par("mar")

#' source a few helper scripts
source("~/Git/UPSCb/src/R/featureSelection.R")
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Create a sample file
#' tab delimited
#' SciLifeID  ID Genotype Replicate Time
#' P8905_117_S67_L007 yourID WT 1 T0
#' ...
#' Upload it to ~/Git/UPSCb/projects/arabidopsis-pib/doc/sample.txt
samples <- read.delim("/mnt/picea/projects/arabidopsis/balbrectsen/arabidopsis-aphids/analysis/kallisto/Sample_names.txt")

#' Read the tx2gene translation
tx2gene <- read.table("/mnt/picea/storage/reference/Arabidopsis-thaliana/ARAPORT11/kallisto/tx2gene.txt",
                      header=TRUE)

#' Get the data
kallisto <- list.files("kallisto", 
                    recursive = TRUE, 
                    pattern = "abundance.tsv",
                    full.names = TRUE)

#' name them
names(kallisto) <- sub("_sortmerna.*","",
                    sapply(strsplit(kallisto, "/"), .subset, 2))

#' match the samples to the files
kallisto <- kallisto[match(samples$SciLifeID,names(kallisto))]

#' load the data
tx <- suppressMessages(tximport(files = kallisto, type = "kallisto", txOut = TRUE))
kt <- round(tx$counts)

#'  And summarize it at the gene level
counts <- round(summarizeToGene(txi=tx,tx2gene=tx2gene)$counts)

#' Combine the technical replicates together
counts <- do.call(cbind,lapply(split.data.frame(t(counts),samples$SampleID),colSums))
#dimnames(scounts)

#' Create a directory
dir.create("analysis/kallisto",recursive = TRUE, showWarnings = FALSE)
write.csv(counts,file="analysis/kallisto/raw-unormalised-gene-expression_data_combtechrep_batch.csv")
length(rowSums(counts)) #32833 genes
length(which(rowSums(counts) == 0)) #7295 genes without counts
geneswithcount <- which(rowSums(counts) != 0) #genes with expression
write.csv(geneswithcount,file="analysis/kallisto/genes_with_counts.csv")

#' The cumulative transcript coverage is deep, about 1000X
plot(density(log10(rowMeans(counts))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by genotype. The samples are extremely similar
#' 
plot.multidensity(lapply(1:ncol(counts),function(k){log10(counts)[,k]}),
                  col=pal[as.integer(samples$Genotype)],
#                  col=sample(pal,size = ncol(counts),replace = TRUE),
                  legend.x="topright",
                  legend=levels(samples$Genotype),
                  legend.col=pal[1:nlevels(samples$Genotype)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' select only one sample per biological replicate
samples2 <- samples[grep("L001", samples$SciLifeID),]

#' Create the DESeq object
samples2$Treatment <- relevel(samples2$Treatment,"c")
samples2$Batch <- factor(samples2$Batch)
dds <- DESeqDataSetFromMatrix(counts,samples2[order(samples2$SampleID),],~Batch+Treatment*Genotype)

#' Attach sample information
colnames(dds)
samples <- read.delim("/mnt/picea/projects/arabidopsis/balbrectsen/arabidopsis-aphids/analysis/kallisto/Sample_names.txt")
samples <- samples[grep("L001", samples$SciLifeID),]
samples <- samples[order(samples$SampleID),]
samples$lab <- paste(samples$SampleID, samples$Genotype, samples$Treatment, sep = "_")
colnames(dds) <- samples$lab
#' Columns in order of treatment (c-4-8)
dds <- dds[,c(1,8,10,7,11,18,2,3,12,
              4,9,15,5,13,16,6,14,17)]

save(dds,file="analysis/kallisto/DESeq-object_combtechrep_batch.rda")
dds
colData(dds)
stopifnot(all(colnames(counts) == rownames(colData(dds))))

#' Estimating the size factor
#' 
#' There is little difference in the size factor.
dds <- estimateSizeFactors(dds)
boxplot(colData(dds)$sizeFactor,main="Library size factor",ylab="proportion")
abline(h=1,lty=2)

#' ## Variance Stabilising Transformation
#' 
#' Since there is almost no difference in library size, a VST is
#' perfectly applicable.

#' This is the situation prior to normalisation
meanSdPlot(log2(counts(dds)[rowSums(counts(dds)) > 0,]))

#' Normalisation - blind, we give no prior as we want to assess quality
vst <- varianceStabilizingTransformation(dds,blind=FALSE)

#' Extract the normalised counts
vsd <- assay(vst)

#' Look at the VST fit. It looks good, almost flat, around 0.25 sd on average
meanSdPlot(vsd[rowSums(counts(dds)) > 0,])

#' The VST introduces an offset
range(vsd)

#' Which we remove, so that 0 means no expression
vsd <- vsd - min(vsd)
write.csv(vsd, "analysis/kallisto/raw-normalized-gene-expression_data_combtechrep_batch.csv")

#' # Quality Assessment
#' ## Principal Component Analysis
pc <- prcomp(t(vsd))
percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
#' 
samples2$PCHgeno <- ifelse(samples2$Genotype == "Col-0", 19, 17)

scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples2$Genotype)],
              #color=sample(pal,ncol(counts),replace=TRUE),
              pch=19)
legend("topleft",pch=19,
       col=pal[1:nlevels(samples2$Genotype)],
       legend=levels(samples2$Genotype))

scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[samples2$Treatment],
              pch=19)
legend("topleft",pch=19,
       col=pal[1:nlevels(as.factor(samples2$Treatment))],
       legend=levels(as.factor(samples2$Treatment)))

scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[samples2$Treatment],
              pch=samples2$PCHgeno)
legend("topleft",pch=c(rep(19,4),17),
       col=c(pal[1:nlevels(as.factor(samples2$Treatment))],"gray","gray"),
       legend=c(levels(as.factor(samples2$Treatment)),"Col-0","pae9-2"))

scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples2$Batch)],
              #color=sample(pal,ncol(counts),replace=TRUE),
              pch=19)
legend("topleft",pch=19,
       col=pal[1:nlevels(as.factor(samples2$Batch))],
       legend=levels(as.factor(samples2$Batch)))

par(mar=mar)

#' General: some variation between biological replicates
#' Most variation explained by treatment, not genotype. Is to be expected due to overruling aphid influence.
#' c and 4 hpi are relatively similar, which is to be expected after only 4 h.
#' Genotypes are separated by PC2 and differ mostly in treatment c.
#' After 8 h aphids, differences between genotypes diminish

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     pch=19,
     col=pal[as.integer(samples2$Treatment)])
legend("topleft",pch=19,
       col=pal[1:nlevels(as.factor(samples2$Treatment))],
       legend=levels(as.factor(samples2$Treatment)))
       

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     pch=19,
     col=pal[as.integer(samples2$Genotype)])
legend("topleft",pch=19,
       col=pal[1:nlevels(as.factor(samples2$Genotype))],
       legend=levels(as.factor(samples2$Genotype)))

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     pch=samples2$PCHgeno,
     col=pal[as.integer(samples2$Treatment)])
legend("topleft",pch=c(rep(19,4),17),
       col=c(pal[1:nlevels(as.factor(samples2$Treatment))],"gray","gray"),
       legend=c(levels(as.factor(samples2$Treatment)),"Col-0","pae9-2"))

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     pch=19,
     col=pal[as.integer(samples2$Batch)])
legend("topleft",pch=19,
       col=pal[1:nlevels(as.factor(samples2$Batch))],
       legend=levels(as.factor(samples2$Batch)))

#' ## Select the genes that are expressed 
sels <- sapply(1:10,function(i){
  featureSelect(vsd,conditions = factor(paste0(samples2$Genotype,samples2$Treatment)),
              exp=i)})


plot(colSums(sels),type="l",xlab="vst cutoff",
     main="number of genes selected at cutoff",ylab="number of genes")

sel <- sels[,2]

#' Hierarchical clustering of the data
plot(hclust(dist(t(vsd[sel,]))),labels=paste(samples2$SampleID,samples2$Genotype,
                                             samples2$Treatment,samples2$Batch))


#' Create a heatmap
hpal <- colorRampPalette(c("blue","white","red"))(100)

heatmap.2(vsd[sel,],trace="none",col=hpal)

#' z scale for gene counts (difference between mean expression and sample expression)
s.vst <- t(scale(t(vsd)))

library(hyperSpec)

heatmap.2(s.vst[sel,],distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",col=hpal,labRow = FALSE,
          labCol = paste(samples2$Genotype,samples2$Treatment,samples2$Batch,sep="-"))









