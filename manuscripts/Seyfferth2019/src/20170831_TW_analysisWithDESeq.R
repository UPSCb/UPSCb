library(RColorBrewer)
input.file = file.choose()
count_table = read.delim(input.file, header = T, row.names = 1)

head(count_table)

meta<-data.frame(row.names=colnames(count_table),genotype=c(rep("T89",6),rep("LMX5",6)),condition=c(rep("NW",3),rep("TW",3), rep("NW",3),rep("TW",3)))

conditions<-data.frame(row.names=colnames(count_table),condition=c(rep("T89NW",3),rep("T89TW",3),rep("LMX5NW",3),rep("LMX5TW",3)))

library(DESeq2)
#' Create the dds object, without giving any prior on the design
conditions <- colnames(count_table)
dds <- DESeqDataSetFromMatrix(
  countData = count_table,
  colData = data.frame(condition=conditions,
                       treatment=meta$condition),
  design = ~ condition)

design(dds)=~treatment

dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count_table)
library(pander)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count_table)
vst <- vst - min(vst)
sel <- rowSums(sapply(lapply(
  split.data.frame(t(count_table >= 2),conditions)
  ,colSums), ">=", 2)) >= 1

allMerged <- merge(vst, sel, by.x = "row.names", by.y = "row.names")

write.table(allMerged, file="20190329_vst_TW_wo35S.txt", row.names=T, col.names=T, sep="\t", quote=F)

library("DESeq")
##Quality check on all data incl. normalization
sum(rowSums(is.na(count_table)) > 0)
cdsFull<-newCountDataSet(count_table,as.character(conditions))
cdsFull <- estimateSizeFactors(cdsFull)
sizeFactors(cdsFull)
cdsFullBlind <- estimateDispersions(cdsFull, method="blind")
cdsFull <- estimateDispersions(cdsFull)

vsdFull <- varianceStabilizingTransformation(cdsFullBlind)
library("RColorBrewer")
library("gplots")

select = order(rowMeans(counts(cdsFull)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
png(filename="Heatmap_top30_vst_TWdata.png")
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
dev.off()

pc <- prcomp(t(vst))
summary(pc)

percent <- round(summary(pc)$importance[2,]*100)

##Sample to Sample distance
dists = dist( t( exprs(vsdFull) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(conditions))

postscript(file="sample-to-sample-distances_TWdata.eps", width=25.6,height=14.4)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

##PCA plot
cdsM<-newCountDataSet(count_table,meta)
cdsM <- estimateSizeFactors(cdsM)
sizeFactors(cdsM)
cdsM <- estimateDispersions(cdsM, method="blind")
vsdM <- varianceStabilizingTransformation(cdsM)
postscript(file="PCA_plot_TW-vs-NW.eps", width=25.6,height=16.4)
print(plotPCA(vsdM, intgroup=c("condition", "genotype")))
dev.off()

##go to sample wise analysis 
png(filename="dispersion_estimation_1E-TW-NW.png", width=1920,height=1080,units="px",pointsize=20)
plotDispEsts(cdsFull)
dev.off()

source("C://Users/Carolin/Box Sync/1_Projects_Caro/Scripts/HelperFilesUPSC/plotDispLSD.R")
plotDispLSD(cdsFull)

##Differential expression and Histogram on log2 and p-Values
res_WT = nbinomTest(cdsFull, "T89NW", "T89TW")
write.table(res_WT, file="20171026_Res_T89_NW_TW.txt", row.names=T, col.names=T, sep="\t", quote=F)


res_WT <- subset(res_WT, padj < 0.05)
DEGs_T89_NW_TW_2fold <- subset(res_WT, res_WT[, "log2FoldChange"] < -0.5 | res_WT[, "log2FoldChange"] > 0.5)
write.table(DEGs_T89_NW_TW_2fold, file="20171023_DEGs_T89_NW_TW_05.txt", row.names=T, col.names=T, sep="\t", quote=F)


##Differential expression and Histogram on log2 and p-Values
res_LMX5 = nbinomTest(cdsFull, "LMX5NW", "LMX5TW")
write.table(res, file="20171026_Res_LMX5_NW_TW.txt", row.names=T, col.names=T, sep="\t", quote=F)

#plotMA(res)
#plot(res$baseMean,res$log2FoldChange,log="x",pch=ifelse(res$padj<.01,19,20),col=ifelse(res$padj<.01,"red","black"))
#hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
res_LMX5 <- subset(res_LMX5, padj < 0.05)
DEGs_LMX5_NW_TW_2fold <- subset(res_LMX5, res_LMX5[, "log2FoldChange"] < -0.5 | res_LMX5[, "log2FoldChange"] > 0.5)
write.table(DEGs_LMX5_NW_TW_2fold, file="20171023_DEGs_LMX5_NW_TW_05.txt", row.names=T, col.names=T, sep="\t", quote=F)


#######################################################################################
#############################Heatmaps##################################################
#######################################################################################

input.file = file.choose()
heat = read.delim(input.file, header = T, row.names = 1)

selheat <- merge(vst_l, heat, by.x = "row.names", by.y = "GeneID")

heat <- selheat[,c("T89_TW_1", "T89_TW_2", "T89_TW_3",
                   "LMX5_TW_1", "LMX5_TW_2", "LMX5_TW_3")]

hpal <- colorRampPalette(c("blue","white","red"))(50)

hclustfunc <- function(x) hclust(x, method="ward.D2")
distfunc <- function(x) dist(x, method="euclidean")


heatmap.2(as.matrix(heat),
               trace="none", hclustfun=hclustfunc, distfun = distfunc,
               col=hpal,Colv=F, dendrogram = "row")

heat <- heat[-1]
heatmap.2(as.matrix(heat), col=hpal,
          trace="none",Colv=F, Rowv=F, dendrogram = "none")
