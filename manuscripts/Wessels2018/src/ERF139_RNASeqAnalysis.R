##############################################################################################
############                                 ########################################################
############  Analysis for DRGs ERF139SRDX   #######################################
############                                 ########################################################
##############################################################################################


count_table <- read.csv("C:/Users/Carolin/Box Sync/1_Projects/ERFs/Bioinformatic part/ERF139/HTSEq/HTSeq_ERF139SRDX.txt", sep="\t", row.names = 1)

count_table <- as.matrix(count_table)
storage.mode(count_table) = "integer"
head (count_table)
meta<-data.frame(row.names=colnames(count_table),condition=c(rep("T89",3), rep("ERF139SRDX",4)))

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count_table,
  colData = data.frame(condition=meta$condition),
  design = ~ condition)

dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count_table)
library(pander)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#Check how many genes are never expressed
sel <- rowSums(count_table) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count_table),digits=1),
        sum(sel),
        nrow(count_table))
[1] "21.7% percent (8990) of 41335 genes are not expressed"

#Display the per-gene mean expression i.e. the mean raw count of every gene across samples is calculated and displayed on a log10 scale.
plot(density(log10(count_table[,1])))
lines(density(log10(count_table[,2])),col="#4d4d4d")
lines(density(log10(count_table[,3])),col="#01665e")
lines(density(log10(count_table[,4])),col="#35978f")
lines(density(log10(count_table[,5])),col="#80cdc1")
lines(density(log10(count_table[,6])),col="#8c510a")
lines(density(log10(count_table[,7])),col="#bf812d")
abline(v=log10(10))

#The cumulative coverage in our samples is good
#sequencing depth that we achieved is ok 


library("pander")
library("vsn")

vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count_table)
vst <- vst - min(vst)
write.table(vst, file="20180514_vst_ERF139SRDX.txt", row.names=T, col.names=T, sep="\t", quote=F)
#Validate the VST
meanSdPlot(vst[rowSums(count_table)>0,])
#Visualize the corrected mean - sd relationship. 
#It is fairly linear, meaning we can assume homoscedasticity. 
#First perform a Principal Component Analysis (PCA) 
pc <- prcomp(t(vst))
summary(pc)

percent <- round(summary(pc)$importance[2,]*100)

#Plot the PCA first two dimensions
pal <- brewer.pal(8,"Dark2")
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(meta$condition)],
     main="Principal Component Analysis",sub="variance stabilized counts")


# do boxplot with PC1 and PC2 to see differences between samples
data_for_plot <- cbind(pc$x, meta)
data_for_plot <- as.data.frame(data_for_plot)

boxplot <- ggplot(data_for_plot, aes(x = condition, y= PC1)) #change y for RDA2
z <- boxplot + geom_boxplot() + scale_x_discrete("")
z

library(corrplot)
cor = cor(count_table)
corrplot.mixed(cor(count_table, method = "pearson"))
col1 <- colorRampPalette(c("navy", "white", "#b2182b"))
corrplot(cor, order="alphabet", method="circle", col=col1(10), title="Correlation TW_NW_ACC_Mock", tl.pos="lt", type="upper", tl.col="black", tl.cex=0.6, tl.srt=45, addCoef.col="black", pch.cex=0.1, insig="blank")

library(ape)
data.dist = as.dist(1 - cor)
out = hclust(data.dist, method="complete")
plot(out)
plot(as.phylo(out), type = "fan")

##Differential expression and Histogram on log2 and p-Values

sum(rowSums(is.na(count_table)) > 0)
cdsFull<-newCountDataSet(count_table,as.character(meta$condition))
cdsFull <- estimateSizeFactors(cdsFull)
sizeFactors(cdsFull)
cdsFullBlind <- estimateDispersions(cdsFull, method="blind")
cdsFull <- estimateDispersions(cdsFull)
vsdFull <- varianceStabilizingTransformation(cdsFullBlind)
res = nbinomTest(cdsFull, "T89", "ERF139SRDX")
write.table(res, file="res_T89_ERF139SRDX_noselesction.txt", row.names=T, col.names=T, sep="\t", quote=F)

plotMA(res)
plot(res$baseMean,res$log2FoldChange,log="x",pch=ifelse(res$padj<.05,19,20),col=ifelse(res$padj<.05,"red","black"))
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
res <- subset(res, padj < 0.05)
DEGs_T89_ERF139SRDX_0.5fold <- subset(res, res[, "log2FoldChange"] < -0.5 | res[, "log2FoldChange"] > 0.5)
write.table(DEGs_T89_ERF139SRDX_0.5fold, file="DEGs_T89_ERF139SRDX_0.5fold.txt", row.names=T, col.names=T, sep="\t", quote=F)

### get the vst expression values for your DRGs####
vst <- data.frame(vst)
vst$T89 <- (vst$B.T89.1 + vst$B.T89.2 + vst$B.T89.3)/3
vst$SRDX <- (vst$X681.6 + vst$X681.7 + vst$X681.8 + vst$X681.15)/4

vst_DEGs_SRDX_0.5fold  <- merge(DEGs_T89_ERF139SRDX_0.5fold, vst, by.x = "id", by.y = "row.names")

sum(vst_DEGs_SRDX_0.5fold$log2FoldChange < 0)
[1]173
[1]492

sum(vst_DEGs_SRDX_0.5fold$log2FoldChange > 0)
[1]139
[1]464

heat_DRGs_SRDX <- vst_DEGs_SRDX_0.5fold[,c("id", "T89", "SRDX")]
row.names(heat_DRGs_SRDX) <- heat_DRGs_SRDX[,1]
heat_DRGs_SRDX <- heat_DRGs_SRDX[,-1]
pheatmap(as.matrix(heat_DRGs_SRDX))

min(heat_DRGs_SRDX)
max(heat_DRGs_SRDX)

colors = unique(c(seq(0,1,length=1), seq(1,3,length=2),seq(3,6,length=3), seq(6,9,length=4)))

hmcol <- colorRampPalette(c("white","#fddbc7","#fc9272", "#de2d26")) (length(colors)-1)

pheatmap(as.matrix(heat_DRGs_SRDX), border_color = NA,
         dendrogram="row", clustering_method = "single",
         col=hmcol, display_numbers = F, breaks=colors)


##############################################################################################
############                          ########################################################
############  Analysis for ERF139OE   #######################################
############                          ########################################################
##############################################################################################

count.table <- read.csv("C:/Users/Carolin/Box Sync/1_Projects/ERFs/Bioinformatic part/ERF139/HTSeq_ERF139_wo4.txt", sep="\t", row.names = 1)

count.table <- as.matrix(count.table)
storage.mode(count.table) = "integer"
head (count.table)
meta2<-data.frame(row.names=colnames(count.table),condition=c(rep("T89",2), rep("ERF139",3)))

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=meta2$condition),
  design = ~ condition)

dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
library(pander)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#Check how many genes are never expressed
sel <- rowSums(count.table) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count.table),digits=1),
        sum(sel),
        nrow(count.table))
[1] "26.5% percent (10972) of 41335 genes are not expressed"

#Display the per-gene mean expression i.e. the mean raw count of every gene across samples is calculated and displayed on a log10 scale.
plot(density(log10(count.table[,1])))
lines(density(log10(count.table[,2])),col="#4d4d4d")
lines(density(log10(count.table[,3])),col="#01665e")
lines(density(log10(count.table[,4])),col="#35978f")
lines(density(log10(count.table[,5])),col="#80cdc1")
abline(v=log10(10))

library("DESeq2")
library("pander")
library("vsn")

vsd2 <- varianceStabilizingTransformation(dds, blind=TRUE)
vst2 <- assay(vsd2)
colnames(vst2) <- colnames(count.table)
vst2 <- vst2 - min(vst2)
write.table(vst2, file="20180514_vst_ERF139.txt", row.names=T, col.names=T, sep="\t", quote=F)
#Validate the VST
meanSdPlot(vst2[rowSums(count.table)>0,])
#Visualize the corrected mean - sd relationship. 
#It is fairly linear, meaning we can assume homoscedasticity. The slight initial trend / bump is due to genes having few counts in a few subset of the samples and hence having a higher variability. This is expected.

#First perform a Principal Component Analysis (PCA) of the data to do a quick quality assessment; i.e. replicate should cluster and the first 2-3 dimensions shouldbe explainable by biological means.

pc <- prcomp(t(vst2))
summary(pc)

percent <- round(summary(pc)$importance[2,]*100)

#Plot the PCA first two dimensions
pal <- brewer.pal(8,"Dark2")
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(meta2$condition)],
     main="Principal Component Analysis",sub="variance stabilized counts")


# do boxplot with PC1 and PC2 to see differences between samples
data_for_plot <- cbind(pc$x, meta2)
data_for_plot <- as.data.frame(data_for_plot)

boxplot <- ggplot(data_for_plot, aes(x = condition, y= PC1)) #change y for RDA2
z <- boxplot + geom_boxplot() + scale_x_discrete("")
z

library(corrplot)
cor = cor(count.table)
corrplot.mixed(cor(count.table, method = "pearson"))
col1 <- colorRampPalette(c("navy", "white", "#b2182b"))
corrplot(cor, order="alphabet", method="circle", col=col1(10), title="Correlation TW_NW_ACC_Mock", tl.pos="lt", type="upper", tl.col="black", tl.cex=0.6, tl.srt=45, addCoef.col="black", pch.cex=0.1, insig="blank")

library(ape)
data.dist = as.dist(1 - cor)
out = hclust(data.dist, method="complete")
plot(out)
plot(as.phylo(out), type = "fan")

##Differential expression and Histogram on log2 and p-Values

library("DESeq")
library("gplots")

sum(rowSums(is.na(count.table)) > 0)
cdsFull2<-newCountDataSet(count.table,as.character(meta2$condition))
cdsFull2 <- estimateSizeFactors(cdsFull2)
sizeFactors(cdsFull2)
cdsFullBlind2 <- estimateDispersions(cdsFull2, method="blind")
cdsFull2 <- estimateDispersions(cdsFull2)
vsdFull2 <- varianceStabilizingTransformation(cdsFullBlind2)
res = nbinomTest(cdsFull2, "T89", "ERF139")
write.table(res, file="ress_T89_ERF139_noSelection.txt", row.names=T, col.names=T, sep="\t", quote=F)

plotMA(res)
plot(res$baseMean,res$log2FoldChange,log="x",pch=ifelse(res$padj<.05,19,20),col=ifelse(res$padj<.05,"red","black"))
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
res <- subset(res, padj < 0.05)
DEGs_T89_ERF139_2fold <- subset(res, res[, "log2FoldChange"] < -0.5 | res[, "log2FoldChange"] > 0.5)

write.table(DEGs_T89_ERF139_2fold, file="DEGs_T89_ERF139_0.5fold.txt", row.names=T, col.names=T, sep="\t", quote=F)

### get the vst expression values for your DRGs####
vst2 <- data.frame(vst2)
vst2$T89 <- (vst2$T89.1 + vst2$T89.2)/2
vst2$ERF139 <- (vst2$X250.5 + vst2$X250.8 + vst2$X250.13)/3

vst_DEGs_ERF139_2fold  <- merge(DEGs_T89_ERF139_2fold, vst2, by.x = "id", by.y = "row.names")
write.table(vst_DEGs_ERF139_2fold, file="vst_DEGs_ERF139_2fold.txt", row.names=T, col.names=T, sep="\t", quote=F)

write.table(vst_DEGs_SRDX_2fold, file="vst_DEGs_SRDX_2fold.txt", row.names=T, col.names=T, sep="\t", quote=F)

sum(vst_DEGs_ERF139_2fold$log2FoldChange < 0)
[1]226
sum(vst_DEGs_ERF139_2fold$log2FoldChange > 0)
[1]284

heat_DRGs_ERF139 <- vst_DEGs_ERF139_2fold[,c("id", "T89", "ERF139")]
row.names(heat_DRGs_ERF139) <- heat_DRGs_ERF139[,1]
heat_DRGs_ERF139 <- heat_DRGs_ERF139[,-1]
pheatmap(as.matrix(heat_DRGs_ERF139))

min(heat_DRGs_ERF139)
max(heat_DRGs_ERF139)

colors = unique(c(seq(0,1,length=1), seq(1,3,length=2),seq(3,6,length=3), seq(6,9,length=4)))

hmcol <- colorRampPalette(c("white","#fddbc7","#fc9272", "#de2d26")) (length(colors)-1)

pheatmap(as.matrix(heat_DRGs_ERF139), border_color = NA,
         dendrogram="row", clustering_method = "single",
         col=hmcol, display_numbers = F, breaks=colors)

##############################################################################################
############                                ########################################################
############  Heatmap for primary targets   #######################################
############                                ########################################################
##############################################################################################
PrimaryTargets <- read.csv("C:/Users/Carolin/Box Sync/1_Projects/ERFs/Bioinformatic part/ERF139/DEGs_0.5fold_pAdj05/Venn/Common_negSRDX.txt", sep="\t", row.names = 1)

pheatmap(as.matrix(PrimaryTargets), border_color = NA,
         dendrogram="row", clustering_method = "single",
         display_numbers = F)
