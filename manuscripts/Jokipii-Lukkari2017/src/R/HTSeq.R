
### ==============================
## libs
### ==============================

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(date))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("pander"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(matrixStats))


source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/plotMA.R")
source("~/Git/UPSCb/src/R/volcanoPlot.R")
detach("package:limma")

### ==============================
## set the working dir
### ==============================
setwd("/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/")

### ==============================
## read the sample information
### ==============================
samples <- read.csv("~/Git/UPSCb/projects/spruce-wood-time-series/doc/samplesmod.csv")

### ==============================
## read the HTSeq files in a matrix
### ==============================
res <- mclapply(dir("htseq",pattern="*STAR.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))
names(res) <- gsub("_STAR\\.txt|.*XX_","",dir("htseq",pattern="*STAR.txt"))

### ==============================
## get the count table
### ==============================
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",3))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]
#write.csv(count.table,"analysis/HTSeq/raw-unormalised-data.csv")

## ADD Venn code

rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Phloem" ],colnames(count.table))] >10 ) >0]
rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Xylem" ],colnames(count.table))] >10 ) >0]

pal <- brewer.pal(3,"Dark2")
png("phloem_xylem_Venn.png",width=800,height=800,pointsize = 34)
plot.new()
grid.draw(venn.diagram(list(
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Phloem" ],colnames(count.table))] >10 ) >0],
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Xylem" ],colnames(count.table))] >10 ) >0]
),
filename=NULL,
col=pal[1:2],
category.names=c("Phloem","Xylem")))
dev.off()


Phloem.genes <- setdiff(
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Phloem" ],colnames(count.table))] >10 ) >0],
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Xylem" ],colnames(count.table))] >10 ) >0])

write.csv2(Phloem.genes,"Phloem.genes")

Phloem.genesb <- setdiff(
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Phloem" ],colnames(count.table))] >0 ) >0],
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Xylem" ],colnames(count.table))] >0 ) >0])

write.csv2(Phloem.genesb,"Phloem.genesb")

Xylem.genes <- setdiff(
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Xylem" ],colnames(count.table))] >10 ) >0],
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Phloem" ],colnames(count.table))] >10 ) >0])

write.csv2(Xylem.genes,"Xylem.genes")

Xylem.genesb <- setdiff(
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Xylem" ],colnames(count.table))] >0 ) >0],
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Phloem" ],colnames(count.table))] >0 ) >0])

write.csv2(Xylem.genesb,"Xylem.genesb")

Common.genes <- intersect(
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Phloem" ],colnames(count.table))] >10 ) >0],
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Xylem" ],colnames(count.table))] >10 ) >0])

write.csv2(Common.genes,"Common.genes")

Common.genesb <- intersect(
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Phloem" ],colnames(count.table))] >0 ) >0],
  rownames(count.table)[rowSums(count.table[,match(samples$Filename[ samples$Tissue == "Xylem" ],colnames(count.table))] >0 ) >0])

write.csv2(Common.genesb,"Common.genesb")


### ==============================
## extract the HTSeq stat lines
### ==============================
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- addInfo
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

## as percentages
apply(count.stats,2,function(co){round(co*100/sum(co))})

## read in the unaligned reads
load("/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/STAR-unmapped-reads.rda")
mat
all.stats <- rbind.data.frame(count.stats,
                   unmapped=mat[match(colnames(count.stats),mat[,1]),2])

### ==============================
## plot the stats
### ==============================
pal=brewer.pal(6,"Dark2")[1:nrow(count.stats)]
mar <- par("mar")
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=pal,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,2e+7))
legend("top",fill=pal,legend=rownames(count.stats),bty="n",cex=0.8)
par(mar=mar)

boxplot(unlist(count.stats["aligned",]/colSums(count.stats)),main="aligned reads",ylab="percent aligned",ylim=c(0,1))


###========================================
##Box blots
###========================================
boxplot(t(all.stats))



### ==============================
## 25.7% of the genes are not expressed
## out of a total of 70736 genes
### ==============================
sel <- rowSums(count.table) == 0
sprintf("%s percent",round(sum(sel) * 100/ nrow(count.table),digits=1))
sprintf("of %s genes are not expressed",nrow(count.table))


### ==============================
## display the per-gene mean expression
## i.e. the mean raw count of every
## gene across samples is calculated
## and displayed on a log10 scale
### ==============================
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

### ==============================
## The same is done for the individual
## samples colored by sample type
### ==============================
pal=brewer.pal(8,"Dark2")
plot.multidensity(log10(count.table),col=sample(pal,length(res),TRUE),
                  legend.x="topright",legend.cex=0.5,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")



### =============================================
## For visualization, the data is
## submitted to a variance stabilization
## transformation using DESeq2. The
## dispersion is estimated independently
## of the sample type. Create DESeqDataSet object
### =============================================

conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)


### =======================================
## stabilize variance, comparing rlog ~ vst
### =======================================
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)
vst <- vst - min(vst)

rld <- rlogTransformation(dds, blind=TRUE)
rlt <- assay(rld)
colnames(rlt) <- colnames(count.table)

###===============================
## Check the mRNA pool variance
###===============================

conditions1 <- samples$Tissue[match(colnames(count.table),samples$Filename)]
conditions2 <- samples$Sampling[match(colnames(count.table),samples$Filename)]
count.table2 <- count.table[order(conditions1,conditions2)]

barplot(colSums(count.table2 >0 )) #plot 1#

conditions1 <- samples$Tissue[match(colnames(vst),samples$Filename)]
conditions2 <- samples$Sampling[match(colnames(vst),samples$Filename)]

vst2 <- vst[ , order(conditions1,conditions2)]

library(matrixStats)
barplot(colMeans(vst2)) #plot 2#

barplot(apply(vst2,2,function(co){mean(co[co>0])})) #plot 3#

conditions1 <- samples$Tissue[match(colnames(rlt),samples$Filename)]
conditions2 <- samples$Sampling[match(colnames(rlt),samples$Filename)]


rlt2 <- rlt[ , order(conditions1,conditions2)]


barplot(apply(rlt2,2,function(co){mean(co[co>0])}))# plot 4#



### ==============================
## Visualize the corrected mean - sd
## relationship. It is fairly linear,
## meaning we can assume homoscedasticity.
## the slight initial trend / bump is
## due to genes having few counts in
## a few subset of the samples and hence
## having a higher variability. This is
## expected.
### ==============================


par(mfrow=c(2,1))
meanSdPlot(vst)
meanSdPlot(rlt)
par(mfrow=c(1,1))

write.csv(rlt,"analysis/HTSeq/rlogTransformed_data.csv")

### ==============================
## perform a Principal Component
## Analysis (PCA) of the data
## to do a quick quality assessment
## i.e. replicate should cluster
## and the first dimensions should
## be explainable by biological means.
### ==============================


pc <- prcomp(t(assay(rld)))

#pc <- prcomp(t(vst))
#percent <- round(summary(pc)$importance[2,]*100)
#smpls <- unique(conditions)


### =================================
## plot the PCA 3 first dimensions
### =================================
cols <- brewer.pal(8,"Dark2")
replicates <- factor(samples$Sampling[match(colnames(count.table),samples$Filename)])
replicates <- levels(replicates)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)]
inx <- factor(samples$Sampling[match(colnames(count.table),samples$Filename)],labels=1:23,levels=replicates)


colours <- ifelse(samples$Tissue[match(colnames(count.table),samples$Filename)]=="Phloem",
                  cols[1],
                  cols[2])


mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=colours,
              pch=LETTERS[inx])

scatterplot3d(pc$x[,2],
              pc$x[,3],
              pc$x[,1],
              zlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=colours,
              pch=LETTERS[inx],angle=50)
par(mar=mar)

### =======================
## and the first two dims
### =======================


plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=ifelse(samples$Tissue[match(colnames(count.table),samples$Filename)]=="Phloem",cols[1],cols[2]),
     pch=LETTERS[inx],
     main="Principal Component Analysis",sub="regularized log transformed counts")


### ===============================
## and the 2nd and 3rd dimensions
### ===============================


plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=ifelse(samples$Tissue[match(colnames(count.table),samples$Filename)]=="Phloem",cols[1],cols[2]),
     pch=LETTERS[inx],
     main="Principal Component Analysis",sub="regularized log transformed counts")



###===========================================================
## PCA for June7 and June13 xylem samples, death of early wood
###===========================================================

mysamples <-c("P1789_101","P1789_102","P1789_103","P825_201","P825_202","P825_203")
pc <- prcomp(t(rlt[,mysamples]))
percent <- round(summary(pc)$importance[2,]*100)
smpls <- unique(conditions)


### ==============================
## plot the PCA 3 first dimensions
### ==============================



colours <- ifelse(
  samples$Date[match(mysamples,samples$Filename)]=="June-13-2011",
                  cols[1],
                  cols[3])

colours <- cols[as.integer(factor(as.character(samples$Date[match(mysamples,samples$Filename)])))]

mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=colours,
              pch=19)
par(mar=mar)

### ==============================
## and the first two dims
### ==============================

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=colours,
     pch=19,
     main="Principal Component Analysis",sub="regularized log transformed counts")

### ==============================
## and the 2nd and 3rd dims
### ==============================


plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=colours,
     pch=19,
     main="Principal Component Analysis",sub="regularized log transformed counts")

#==========================
### DE death of early wood
#==========================
mysamples <-c("P1789_101","P1789_102","P1789_103","P825_201","P825_202","P825_203")
conditions <- factor(samples$Sampling[match(mysamples,samples$Filename)])

dds <- DESeqDataSetFromMatrix(
  countData = count.table[,mysamples],
  colData = data.frame(condition=conditions),
  design = ~ condition)

## do the DE
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds)

## how many of every kind is DE at an alpha level of 0.01
alpha <- 0.01


### ==============================
## Check all feature types
### ==============================
## subset to what is expressed
res <- res[!is.na(res$padj),]
resultsNames(dds)
DESeq2::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
hist(res$padj,breaks=seq(0,1,.01))
sum(res$padj<alpha,na.rm=TRUE)

signf <- res[res$padj<alpha,]
signf[order(abs(signf$log2FoldChange),decreasing=TRUE),]
signf[order(signf$padj),]

write.csv2(signf,"signf.csv")

###=========================================================
##PCA for latewood appereance, July 11th vs. July 4th, xylem
###=========================================================

mysamplessec <-c("P825_210","P825_211","P825_212","P825_213","P825_214","P825_215")
pc <- prcomp(t(rlt[,mysamplessec]))
percent <- round(summary(pc)$importance[2,]*100)
smpls <- unique(conditions)

### ==============================
## plot the PCA 3 first dimensions
### ==============================

colours <- ifelse(
  samples$Date[match(mysamplessec,samples$Filename)]=="July-4-2011",
  cols[1],
  cols[3])



mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=colours,
              pch=19)
par(mar=mar)

### ==============================
## and the first two dims
### ==============================


plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=colours,
     pch=19,
     main="Principal Component Analysis",sub="regularized log transformed counts")

### ==============================
## and the 2nd and 3rd dims
### ==============================



plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=colours,
     pch=19,
     main="Principal Component Analysis",sub="regularized log transformed counts")



###==================================================================================================
##PCA for latewood appereance, samples having latewood vs. the ones without, xylem
###==================================================================================================


mysamplesfour <-c("P825_213","P825_214","P825_215","P825_216","P825_217","P825_218","P825_219","P825_220","P825_221")
pc <- prcomp(t(rlt[,mysamplesfour]))
percent <- round(summary(pc)$importance[2,]*100)
smpls <- unique(conditions)


### ==============================
## plot the PCA 3 first dimensions
### ==============================
colours <- ifelse(
  samples$Dev.keyword[match(mysamplesfour,samples$Filename)]=="ALW",
  cols[3],
  cols[1])


mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=colours,
              pch=19)
par(mar=mar)

### ==============================
## and the first two dims
### ==============================

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=colours,
     pch=19,
     main="Principal Component Analysis",sub="regularized log transformed counts")


### ==============================
## and the 2nd and 3rd dims
### ==============================

plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=colours,
     pch=19,
     main="Principal Component Analysis",sub="regularized log transformed counts")



#==============================
### DE latewood vs. not, xylem
#==============================
mysamplesfour <-c("P825_213","P825_214","P825_215","P825_216","P825_217","P825_218","P825_219","P825_220","P825_221")
conditions <- samples$Dev.keyword[match(mysamplesfour,samples$Filename)]

dds <- DESeqDataSetFromMatrix(
  countData = count.table[,mysamplesfour],
  colData = data.frame(condition=conditions),
  design = ~ condition)

## do the DE
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds)

## how many of every kind is DE at an alpha level of 0.01
alpha <- 0.01


### ==============================
## Check all feature types
### ==============================
## subset to what is expressed
res <- res[!is.na(res$padj),]
resultsNames(dds)
DESeq2::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
hist(res$padj,breaks=seq(0,1,.01))
sum(res$padj<alpha,na.rm=TRUE)

signf1 <- as.data.frame(res)[res$padj<alpha,]
signf1[order(abs(signf1$log2FoldChange),decreasing=TRUE),]
signf1[order(signf1$padj),]

write.csv2(signf1,"signf1.csv")



#===============================
### DE July 11th vs. 4th, xylem
#===============================
mysamplessec <-c("P825_210","P825_211","P825_212","P825_213","P825_214","P825_215")
conditions <- samples$Date[match(mysamplessec,samples$Filename)]

dds <- DESeqDataSetFromMatrix(
  countData = count.table[,mysamplessec],
  colData = data.frame(condition=conditions),
  design = ~ condition)

## do the DE
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds,contrast=c("condition","July-11-2011", "July-4-2011"))


## how many of every kind is DE at an alpha level of 0.01
alpha <- 0.01


### ==============================
## Check all feature types
### ==============================
## subset to what is expressed
res <- res[!is.na(res$padj),]
## resultsNames(dds)
DESeq2::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
hist(res$padj,breaks=seq(0,1,.01))
sum(res$padj<alpha,na.rm=TRUE)


signf2 <- as.data.frame(res)[res$padj<alpha,]
signf2[order(abs(signf2$log2FoldChange),decreasing=TRUE),]
signf2[order(signf2$padj),]


write.csv2(signf2,"signf2.csv")



###==========================================================
## 79 overlapping gene models between comparisons
###==========================================================
overlap <- intersect(rownames(signf1),rownames(signf2))

write.csv2(overlap,"overlap.csv")


###=========================================================
## Venn diagram latewood
###=========================================================
library(venneuler)

v <- venneuler(c(A=640, B=467, 'A&B'=79))
plot(v)


###==============================================================
# Heat map, putatively cell death related genes, only xylem
###==============================================================


pcdgenes <- c("MA_10430708g0010","MA_10433349g0010","MA_183122g0010","MA_10430487g0010","MA_9247617g0010","MA_10435064g0010","MA_103463g0010",
              "MA_9371025g0010","MA_34759g0010","MA_96038g0020","MA_247624g0010","MA_210951g0010","MA_10429635g0010",
              "MA_280780g0010","MA_22695g0010","MA_895155g0010","MA_8445917g0010","MA_95383g0020","MA_10433186g0010",
              "MA_536873g0010","MA_75204g0010","MA_9316459g0010","MA_158988g0010")

xylemsamples <- c(
  "P1789_101","P1789_102","P1789_103","P825_201","P825_202","P825_203","P825_204","P825_205","P825_206",
  "P825_207","P825_208","P825_209","P825_210","P825_211","P825_212","P825_213","P825_214","P825_215",
  "P825_216","P825_217","P825_218","P825_219","P825_220","P825_221","P825_222","P825_223","P825_224",
  "P825_225","P825_226","P825_228","P825_229","P825_230","P825_231","P825_232","P825_233","P825_234",
  "P825_236","P825_237","P825_238","P825_239","P825_240","P825_241","P825_242","P825_243","P825_244",
  "P825_245","P825_246","P825_247","P825_248","P825_249","P825_250","P825_251","P825_252","P825_253",
  "P825_254","P825_255","P825_256","P825_257","P825_258","P825_259","P825_260","P825_261","P825_262",
  "P825_263","P825_264")

pcdmatrix <- scale(rlt[pcdgenes,xylemsamples])


colnames(pcdmatrix) <- samples[match(colnames(pcdmatrix),samples$Filename),"Date"]

# creates a own color palette from blue to red
my_palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 99)

pdf("~/Git/UPSCb/projects/spruce-wood-time-series/doc/heat_PCD.pdf", width = 12, height = 12, onefile = FALSE)

# create heatmap
heatmap.2(pcdmatrix,Colv=NA,col=my_palette,keysize=0.7,density.info="none",
          trace="none",margins=c(7,9))

dev.off()


###==============================================================================================================================
## Heat map, all full length lignin genes. Genes without expression in xylem discarded
###==============================================================================================================================
ligninge2 <- c("MA_15852g0010","MA_10429279g0010","MA_44561g0010","MA_123220g0010","MA_118702g0010",
               "MA_130482g0010","MA_14460g0010","MA_10426159g0010","MA_10430515g0010",
               "MA_8946476g0010","MA_70509g0010","MA_633407g0010","MA_56692g0010","MA_109119g0010","MA_106573g0010",
               "MA_109033g0010","MA_10435234g0010","MA_175723g0010","MA_137442g0010",
               "MA_10055697g0010","MA_334200g0010","MA_9103099g0010","MA_10433066g0010","MA_7247276g0010","MA_10432446g0030",
               "MA_9992421g0010","MA_59638g0010","MA_10435850g0010","MA_13550g0010","MA_10427193g0010","MA_10435550g0030",
               "MA_28222g0010","MA_10427075g0010","MA_10434709g0010","MA_53288g0010","MA_10429973g0020",
               "MA_10434848g0010","MA_24539g0010","MA_10432870g0010","MA_158072g0010","MA_10430254g0010",
               "MA_52987g0020","MA_61321g0010","MA_17680g0010","MA_54958g0010","MA_109548g0010",
               "MA_99750g0010","MA_10425888g0010","MA_10432610g0020","MA_949441g0010","MA_16731g0010","MA_10427515g0010",
               "MA_362678g0010","MA_6931g0010","MA_9065834g0010","MA_222430g0010","MA_832501g0010","MA_137109g0010",
               "MA_166604g0010","MA_68461g0010","MA_10435631g0010","MA_9446650g0010","MA_10435810g0010","MA_10436663g0010",
               "MA_29397g0010","MA_46269g0010","MA_10426788g0020","MA_76956g0010","MA_165197g0010","MA_363801g0010",
               "MA_10432099g0010","MA_10436313g0010","MA_104401g0010","MA_10427985g0010",
               "MA_423264g0010","MA_10436231g0020","MA_18380g0010","MA_135446g0010","MA_135525g0010","MA_10432110g0010",
               "MA_93032g0010","MA_9470692g0010","MA_62368g0010","MA_87599g0010","MA_47348g0010",
               "MA_811078g0010","MA_170004g0010","MA_67291g0010","MA_132857g0010","MA_14717g0010","MA_10434090g0010",
              "MA_185257g0010","MA_192464g0010","MA_57140g0010","MA_192698g0010",
              "MA_57140g0020","MA_76578g0010","MA_10434864g0010","MA_10432715g0010","MA_10427033g0010","MA_118833g0010",
              "MA_75861g0010","MA_33057g0010","MA_66348g0010","MA_122769g0010","MA_28768g0010","MA_115430g0010",
              "MA_14240g0010","MA_205016g0010","MA_6005856g0010","MA_10430775g0010","MA_10435976g0010","MA_186345g0010",
              "MA_10432379g0020","MA_256822g0010","MA_66808g0010","MA_185360g0010","MA_74620g0010","MA_656798g0010",
              "MA_124869g0010",
              "MA_140596g0010","MA_108425g0010","MA_109058g0010","MA_10227622g0010","MA_125040g0020","MA_45356g0010",
              "MA_10431506g0010","MA_235485g0010","MA_111431g0010","MA_10434028g0020","MA_3486g0010",
              "MA_96757g0010","MA_91956g0010","MA_10425995g0010","MA_10437025g0010",
              "MA_177017g0010","MA_10427936g0010","MA_10432809g0020","MA_195775g0010","MA_6833274g0010",
              "MA_412995g0010","MA_183305g0010","MA_41416g0010","MA_10g0010",
              "MA_166386g0010","MA_9873048g0010","MA_941794g0010","MA_119796g0010","MA_132958g0020",
              "MA_10428075g0010","MA_10431507g0020",
              "MA_10435488g0020","MA_444g0010","MA_5316g0010",
              "MA_10433564g0020","MA_83406g0010","MA_91294g0020","MA_93867g0010",
              "MA_10426168g0010","MA_10434825g0010","MA_10429595g0010")

xylemsamples <- c("P1789_101","P1789_102","P1789_103","P825_201","P825_202","P825_203","P825_204","P825_205","P825_206",
                  "P825_207","P825_208","P825_209","P825_210","P825_211","P825_212","P825_213","P825_214","P825_215",
                  "P825_216","P825_217","P825_218","P825_219","P825_220","P825_221","P825_222","P825_223","P825_224",
                  "P825_225","P825_226","P825_228","P825_229","P825_230","P825_231","P825_232","P825_233","P825_234",
                  "P825_236","P825_237","P825_238","P825_239","P825_240","P825_241","P825_242","P825_243","P825_244",
                  "P825_245","P825_246","P825_247","P825_248","P825_249","P825_250","P825_251","P825_252","P825_253",
                  "P825_254","P825_255","P825_256","P825_257","P825_258","P825_259","P825_260","P825_261","P825_262",
                  "P825_263","P825_264")

ligningmatrix2 <- scale(rlt[ligninge2,xylemsamples])

colnames(ligningmatrix2) <- samples[match(colnames(ligningmatrix2),samples$Filename),"Date"]

# creates a own color palette from blue to red
my_palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 99)

pdf("~/Git/UPSCb/projects/spruce-wood-time-series/doc/heat_lignin.pdf", width = 24, height = 13, onefile = FALSE)

my.breaks <- c(seq(-3,3,length.out=100))

# create heatmap
heatmap.2(ligningmatrix2,Colv=NA,col=my_palette,keysize=0.7,density.info="none",  
          trace="none",cexRow=0.7,cexCol=0.7,margins=c(7,7),breaks=my.breaks)

dev.off()




###====================================
#Expression of all lignin clusters
###====================================

lgenes <-data.frame("lgenes" = c("MA_6931g0010","MA_109548g0010","MA_130482g0010","MA_118702g0010",
                                 "MA_10432099g0010","MA_9446650g0010","MA_87599g0010","MA_106573g0010",
                                 "MA_362678g0010","MA_123220g0010","MA_56692g0010","MA_10436231g0020","MA_10435810g0010","MA_10427515g0010","MA_104401g0010",
                                 "MA_10433066g0010","MA_61321g0010","MA_62368g0010","MA_91956g0010","MA_6833274g0010","MA_412995g0010","MA_76578g0010","MA_57140g0010",
                                 "MA_192464g0010","MA_33057g0010","MA_6005856g0010","MA_175723g0010","MA_118833g0010",
                                 "MA_205016g0010","MA_10434090g0010","MA_10432715g0010"),
                    "family" = c("CCoAOMT","C3H","C4H","C4H","COMT","CCR","CSE","HCT","CCoAOMT","PAL","4CL","CAD","CCR2","C3H2",
                                 "CAD","C3H2","F5H","CSE2","PRX","PRX","PRX","LAC","LAC","LAC","LAC","LAC","CYP","LAC","LAC","LAC","LAC"),
                    "cluster" = paste("Cluster", c(rep(1, 11), rep(2, 7), rep(3, 13))),
                    stringsAsFactors=FALSE)

cluter_family_map <- tapply(lgenes$family, lgenes$cluster, function(x){unique(x)})
lgenes$lty <- apply(lgenes,1,function(x){
  which(cluter_family_map[[x[3]]] %in% x[2])
})
samples <- samples[match(colnames(count.table),samples$Filename),]


cbPalette <- c("gray47","saddlebrown","blue3", "brown3", "chartreuse3", "cadetblue3", "cyan3", "paleturquoise4", "aquamarine3",
               "deeppink3","darkseagreen3","darkorchid3","red3","gray22","yellow3","dodgerblue3",
               "burlywood4","plum3","salmon","violetred3","springgreen3","slategray4","peru","darkslategrey","darkorange3",
               "midnightblue","yellowgreen","forestgreen","palevioletred","goldenrod4","royalblue")


plotCandidatesDate <- function(lgenes){
  subdf <- as.data.frame(rlt[lgenes[,1],])
  colnames(subdf) <- samples$ID
  subdf <- as.data.frame(cbind(subdf,lgenes))
  sub.rlt <- melt(subdf, id.vars = c("lgenes", "family", "cluster", "lty")) #Convert to row major data frame format
  sub.rlt[,5] <- as.character(sub.rlt[,5])
  d <- mdy(as.character(samples$Date)[match(sub.rlt[,5],samples$ID)])
  o <- order(d)
  d.dmy <- paste(day(d),month(d),year(d),sep="-")
  sub.rlt$Date <- factor(d.dmy, levels=unique(d.dmy[o]))
  sub.rlt$Tissue <- factor(samples[match(sub.rlt[,5],samples$ID),"Tissue"])
  sub.rlt$Filename <- samples[match(sub.rlt[,5],samples$ID),"Filename"]
  colnames(sub.rlt) <- c("Gene","Family","Cluster","Line","Sample","Expression",
                         "Date","Tissue","Filename") #Set column names
  line_labels <- lapply(
    lapply(
      tapply(
        sub.rlt$Family, sub.rlt$Line, paste),
      unique),
    function(f){
      l<-length(f)
      paste("C", seq(l),":",f,sep="")
    }
  )
  line_labels <- unlist(lapply(line_labels,paste,collapse=", "))
  sub.rlt <- sub.rlt[sub.rlt$Tissue=="Xylem",]
  sub.rlt$Family_Gene <- paste(sub.rlt$Family,sub.rlt$Gene,sep="_")
  pl <- ggplot(sub.rlt, aes(x=Date,group=Gene,y=Expression)) + #Base layer for plot
    stat_summary(fun.y = mean,geom = "line", lwd=1,
                 aes(color = Gene, lty = factor(Line))) +
    scale_colour_manual(values=cbPalette) +
    scale_linetype(name = "Cluster:Family", labels = line_labels) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    facet_grid(Cluster~., scales="free_y") +
    guides(col = guide_legend(nrow = 3 ),
           lty = guide_legend(nrow = 2, keywidth = unit(2, "cm"))) +
    theme(legend.position = "bottom", panel.background = element_blank(),
          panel.grid = element_blank(), panel.border = element_rect(fill=NA)) +
    scale_x_discrete(expand = c(0, 0))
  return(pl)
}

plotCandidatesDate(lgenes)


###=======================================================
#Expression of all lignin clusters without gene info
###=======================================================

lgenes <-data.frame("lgenes" = c("MA_6931g0010","MA_109548g0010","MA_130482g0010","MA_118702g0010",
                                 "MA_10432099g0010","MA_9446650g0010","MA_87599g0010","MA_106573g0010",
                                 "MA_362678g0010","MA_123220g0010","MA_56692g0010","MA_10436231g0020","MA_10435810g0010","MA_10427515g0010","MA_104401g0010",
                                 "MA_10433066g0010","MA_61321g0010","MA_62368g0010","MA_91956g0010","MA_6833274g0010","MA_412995g0010","MA_76578g0010","MA_57140g0010",
                                 "MA_192464g0010","MA_33057g0010","MA_6005856g0010","MA_175723g0010","MA_118833g0010",
                                 "MA_205016g0010","MA_10434090g0010","MA_10432715g0010"),
                    "family" = c("CCoAOMT","C3H","C4H","C4H","COMT","CCR","CSE","HCT","CCoAOMT","PAL","4CL","CAD","CCR2","C3H2",
                                 "CAD","C3H2","F5H","CSE2","PRX","PRX","PRX","LAC","LAC","LAC","LAC","LAC","CYP","LAC","LAC","LAC","LAC"),
                    "cluster" = paste("Cluster", c(rep(1, 11), rep(2, 7), rep(3, 13))),
                    stringsAsFactors=FALSE)

samples <- samples[match(colnames(count.table),samples$Filename),]


cbPalette <- c("gray47","saddlebrown","midnightblue", "chartreuse3", "dodgerblue3",
               "deeppink3","darkorchid3","coral3","gray22","goldenrod3",
               "plum4","aquamarine3","khaki3","azure3","darkorange3",
               "cyan3","royalblue")


plotCandidatesDate <- function(lgenes){
  subdf <- as.data.frame(rlt[lgenes[,1],])
  colnames(subdf) <- samples$ID
  subdf <- as.data.frame(cbind(subdf,lgenes))
  sub.rlt <- melt(subdf, id.vars = c("lgenes", "family", "cluster")) #Convert to row major data frame format
  sub.rlt[,4] <- as.character(sub.rlt[,4])
  d <- mdy(as.character(samples$Date)[match(sub.rlt[,4],samples$ID)])
  o <- order(d)
  d.dmy <- paste(day(d),month(d),year(d),sep="-")
  sub.rlt$Date <- factor(d.dmy, levels=unique(d.dmy[o]))
  sub.rlt$Tissue <- factor(samples[match(sub.rlt[,4],samples$ID),"Tissue"])
  sub.rlt$Filename <- samples[match(sub.rlt[,4],samples$ID),"Filename"]
  colnames(sub.rlt) <- c("Gene","Family","Cluster","Sample","Expression",
                         "Date","Tissue","Filename") #Set column names
  sub.rlt <- sub.rlt[sub.rlt$Tissue=="Xylem",]
  sub.rlt$Family_Gene <- paste(sub.rlt$Family,sub.rlt$Gene,sep="_")
  pl <- ggplot(sub.rlt, aes(x=Date,group=Gene,y=Expression)) + #Base layer for plot
    stat_summary(fun.y = mean,geom = "line", lwd=1,
                 aes(color = Family)) +
    scale_colour_manual(values=cbPalette) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1), panel.background = element_blank(),
          panel.grid = element_blank(), panel.border = element_rect(fill=NA)) +
    scale_x_discrete(expand = c(0, 0)) +
    facet_grid(Cluster~., scales="free_y")
  return(pl)
}


pdf("~/Git/UPSCb/projects/spruce-wood-time-series/doc/plot_lignin.pdf", width = 8, height = 13, onefile = FALSE)  
plotCandidatesDate(lgenes)
dev.off()


###===============================================
## Total lignin in living wood, linear regression
###===============================================

totlignin = data.frame(lignin = c(20.62453,
                                  25.25811,
                                  25.48501,
                                  25.0127,
                                  25.67378,
                                  22.68023,
                                  25.49668,
                                  24.72025,
                                  26.28046,
                                  19.71202,
                                  23.84836,
                                  21.26355,
                                  26.84147,
                                  23.36291,
                                  23.99029,
                                  26.25053,
                                  24.63209,
                                  24.90654,
                                  24.88728,
                                  25.60002,
                                  24.32007,
                                  29.56139,
                                  25.52557,
                                  27.0483,
                                  25.4805,
                                  23.26183,
                                  30.42753,
                                  25.148),
                       daynro = c(1,
                                  1,
                                  1,
                                  16,
                                  16,
                                  16,
                                  30,
                                  30,
                                  30,
                                  57,
                                  57,
                                  87,
                                  87,
                                  87,
                                  115,
                                  115,
                                  115,
                                  155,
                                  155,
                                  155,
                                  184,
                                  184,
                                  184,
                                  219,
                                  219,
                                  247,
                                  247,
                                  247
                       )
)


lig.mod1 = lm(lignin ~ daynro, data = totlignin)


###===============================================
## Total lignin dead wood, linear regression
###===============================================

totlignin = data.frame(lignin = c(28.98292,
                                  29.93209,
                                  29.56214,
                                  29.92202,
                                  31.63774,
                                  28.73988,
                                  31.18209,
                                  31.2291,
                                  31.0205,
                                  28.40422,
                                  31.1398,
                                  30.20112,
                                  27.27941,
                                  28.30992,
                                  28.76324,
                                  27.94533,
                                  27.34749,
                                  28.74167,
                                  27.46825,
                                  26.82445,
                                  26.51981,
                                  25.6461,
                                  27.2398,
                                  27.35525,
                                  32.74898,
                                  28.80224,
                                  27.01786,
                                  28.14662,
                                  30.69162,
                                  31.16474,
                                  29.54658,
                                  27.76631,
                                  30.06107,
                                  29.52142,
                                  25.75068
),
                       daynro = c(1,
                                  1,
                                  1,
                                  15,
                                  15,
                                  15,
                                  30,
                                  30,
                                  30,
                                  45,
                                  45,
                                  72,
                                  72,
                                  72,
                                  102,
                                  102,
                                  102,
                                  130,
                                  130,
                                  130,
                                  170,
                                  170,
                                  170,
                                  199,
                                  199,
                                  199,
                                  234,
                                  234,
                                  234,
                                  262,
                                  262,
                                  262,
                                  274,
                                  274,
                                  274
                                  )
)


lig.mod2 = lm(lignin ~ daynro, data = totlignin)


###====================================================
## Expression profiles of marker genes, phloem/cambium
###====================================================
markers <-data.frame("markers" = c("MA_9897324g0010",
                                   "MA_10020g0010",
                                   "MA_915221g0010",
                                   "MA_106416g0010",
                                   "MA_920994g0010",
                                   "MA_325398g0010"
                                   
),

"family" = c("APL","EXPA","EXPA","CDKB","CDKB","CDKB"
),

stringsAsFactors=FALSE)

samples <- samples[match(colnames(count.table),samples$Filename),]

cbPalette <- c("azure3","saddlebrown","blue3", "salmon","deeppink3","darkorchid3"
)

conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)

## Two functions to define the error bar y positions for the standard error of the mean
sem_ymin <- function(x){
  m <- mean(x)
  s <- sd(x)
  se <- s/sqrt(length(x))
  return(m - se)
}
sem_ymax <- function(x){
  m <- mean(x)
  s <- sd(x)
  se <- s/sqrt(length(x))
  return(m + se)
}


plotCandidatesDate <- function(markers){
  subdf <- as.data.frame(rlt[markers[,1],])
  colnames(subdf) <- samples$ID
  subdf <- as.data.frame(cbind(subdf,markers))
  sub.rlt <- melt(subdf) #Convert to row major data frame format
  sub.rlt[,3] <- as.character(sub.rlt[,3])
  d <- mdy(as.character(samples$Date)[match(sub.rlt[,3],samples$ID)])
  o <- order(d)
  d.dmy <- paste(day(d),month(d),year(d),sep="-")
  sub.rlt$Date <- factor(d.dmy, levels=unique(d.dmy[o]))
  sub.rlt$Tissue <- factor(samples[match(sub.rlt[,3],samples$ID),"Tissue"])
  sub.rlt$Filename <- samples[match(sub.rlt[,3],samples$ID),"Filename"]
  colnames(sub.rlt) <- c("Gene","Family","Sample","Expression","Date","Tissue","Filename") #Set column names
  sub.rlt$Family_Gene <- paste(sub.rlt$Family,sub.rlt$Gene,sep="_")
  pl <- ggplot(sub.rlt, aes(x=Date,group=Gene,y=Expression)) + #Base layer for plot
    stat_summary(fun.y = mean,geom = "line", lwd=1, aes(color = Gene)) + #Add line for the mean
    stat_summary(fun.y = mean, fun.ymin = sem_ymin, fun.ymax = sem_ymax)+ #Add error bars and mean dots#Put all genes in separate facets
    scale_colour_manual(values=cbPalette) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    facet_grid(Family~Tissue,scales="free")#Rotate axis text
  return(pl)
}

plotCandidatesDate(markers)


###=============================================
## Expression profiles of marker genes, xylem
###=============================================
markers <-data.frame("markers" = c("MA_10427616g0010",
                                   "MA_10020g0010",
                                   "MA_915221g0010",
                                   "MA_183130g0010",
                                   "MA_140410g0010",
                                   "MA_10429177g0010"
),

"family" = c("AHP6","EXPA","EXPA","CESA","CESA","CESA"
),
stringsAsFactors=FALSE)

samples <- samples[match(colnames(count.table),samples$Filename),]

cbPalette <- c("yellow3","gray47","plum4","cyan3","forestgreen","darkorange3"
)

conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)


plotCandidatesDate <- function(markers){
  subdf <- as.data.frame(rlt[markers[,1],])
  colnames(subdf) <- samples$ID
  subdf <- as.data.frame(cbind(subdf,markers))
  sub.rlt <- melt(subdf) #Convert to row major data frame format
  sub.rlt[,3] <- as.character(sub.rlt[,3])
  d <- mdy(as.character(samples$Date)[match(sub.rlt[,3],samples$ID)])
  o <- order(d)
  d.dmy <- paste(day(d),month(d),year(d),sep="-")
  sub.rlt$Date <- factor(d.dmy, levels=unique(d.dmy[o]))
  sub.rlt$Tissue <- factor(samples[match(sub.rlt[,3],samples$ID),"Tissue"])
  sub.rlt$Filename <- samples[match(sub.rlt[,3],samples$ID),"Filename"]
  colnames(sub.rlt) <- c("Gene","Family","Sample","Expression","Date","Tissue","Filename") #Set column names
  sub.rlt$Family_Gene <- paste(sub.rlt$Family,sub.rlt$Gene,sep="_")
  pl <- ggplot(sub.rlt, aes(x=Date,group=Gene,y=Expression)) + #Base layer for plot
    stat_summary(fun.y = mean,geom = "line", lwd=1, aes(color = Gene)) + #Add line for the mean
    stat_summary(fun.y = mean, fun.ymin = sem_ymin, fun.ymax = sem_ymax)+ #Add error bars and mean dots#Put all genes in separate facets
    scale_colour_manual(values=cbPalette) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    facet_grid(Family~Tissue,scales="free")#Rotate axis text
  return(pl)
}

plotCandidatesDate(markers)



###=======================================
## Expression profiles of lignin TF genes
###=======================================

gen <-data.frame("gen" = c("MA_10434782g0020",
                             "MA_9483804g0010",
                             "MA_101790g0010",
                             "MA_17311g0010",
                             "MA_7091g0010"
                             
),

"Family" = c("AS2","MYB16","MYB4","MYB4","MYB4"
),
stringsAsFactors=FALSE)

samples <- samples[match(colnames(count.table),samples$Filename),]

cbPalette <- c("gray47","deeppink3","dodgerblue3","peru","yellowgreen"
)

conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)

## Two functions to define the error bar y positions for the standard error of the mean
sem_ymin <- function(x){
  m <- mean(x)
  s <- sd(x)
  se <- s/sqrt(length(x))
  return(m - se)
}
sem_ymax <- function(x){
  m <- mean(x)
  s <- sd(x)
  se <- s/sqrt(length(x))
  return(m + se)
}


plotCandidatesDate <- function(gen){
  subdf <- as.data.frame(rlt[gen[,1],])
  colnames(subdf) <- samples$ID
  subdf <- as.data.frame(cbind(subdf,gen))
  sub.rlt <- melt(subdf) #Convert to row major data frame format
  sub.rlt[,3] <- as.character(sub.rlt[,3])
  d <- mdy(as.character(samples$Date)[match(sub.rlt[,3],samples$ID)])
  o <- order(d)
  d.dmy <- paste(day(d),month(d),year(d),sep="-")
  sub.rlt$Date <- factor(d.dmy, levels=unique(d.dmy[o]))
  sub.rlt$Tissue <- factor(samples[match(sub.rlt[,3],samples$ID),"Tissue"])
  sub.rlt$Filename <- samples[match(sub.rlt[,3],samples$ID),"Filename"]
  colnames(sub.rlt) <- c("Gene","Family","Sample","Expression","Date","Tissue","Filename") #Set column names
  sub.rlt <- sub.rlt[sub.rlt$Tissue=="Xylem",]
  sub.rlt$Family_Gene <- paste(sub.rlt$Family,sub.rlt$Gene,sep="_")
  pl <- ggplot(sub.rlt, aes(x=Date,group=Gene,y=Expression)) + #Base layer for plot
    stat_summary(fun.y = mean,geom = "line", lwd=1, aes(color = Gene)) + #Add line for the mean
    stat_summary(fun.y = mean, fun.ymin = sem_ymin, fun.ymax = sem_ymax)+ #Add error bars and mean dots#Put all genes in separate facets
    scale_colour_manual(values=cbPalette) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1), panel.background = element_blank(),
          panel.grid = element_blank(), panel.border = element_rect(fill=NA),aspect.ratio=0.6) +
    scale_x_discrete(expand = c(0, 0.6))+
  facet_grid(.~ Family, scales="free_y")
  return(pl)
} 


plotCandidatesDate(gen)

###=======================================
## Expression profiles of AIL1 and FTLs
###=======================================

gen <-data.frame("gen" = c("MA_121578g0010",
                           "MA_400747g0010",
                           "MA_5386467g0010"
                           
),

"Family" = c("AIL1","FTL1","FTL2"
),
stringsAsFactors=FALSE)

samples <- samples[match(colnames(count.table),samples$Filename),]

cbPalette <- c("gray47","deeppink3","dodgerblue3","peru"
)

conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)

## Two functions to define the error bar y positions for the standard error of the mean
sem_ymin <- function(x){
  m <- mean(x)
  s <- sd(x)
  se <- s/sqrt(length(x))
  return(m - se)
}
sem_ymax <- function(x){
  m <- mean(x)
  s <- sd(x)
  se <- s/sqrt(length(x))
  return(m + se)
}


plotCandidatesDate <- function(gen){
  subdf <- as.data.frame(rlt[gen[,1],])
  colnames(subdf) <- samples$ID
  subdf <- as.data.frame(cbind(subdf,gen))
  sub.rlt <- melt(subdf) #Convert to row major data frame format
  sub.rlt[,3] <- as.character(sub.rlt[,3])
  d <- mdy(as.character(samples$Date)[match(sub.rlt[,3],samples$ID)])
  o <- order(d)
  d.dmy <- paste(day(d),month(d),year(d),sep="-")
  sub.rlt$Date <- factor(d.dmy, levels=unique(d.dmy[o]))
  sub.rlt$Tissue <- factor(samples[match(sub.rlt[,3],samples$ID),"Tissue"])
  sub.rlt$Filename <- samples[match(sub.rlt[,3],samples$ID),"Filename"]
  colnames(sub.rlt) <- c("Gene","Family","Sample","Expression","Date","Tissue","Filename") #Set column names
  sub.rlt <- sub.rlt[sub.rlt$Tissue=="Phloem",]
  sub.rlt$Family_Gene <- paste(sub.rlt$Family,sub.rlt$Gene,sep="_")
  pl <- ggplot(sub.rlt, aes(x=Date,group=Gene,y=Expression)) + #Base layer for plot
    stat_summary(fun.y = mean,geom = "line", lwd=1, aes(color = Gene)) + #Add line for the mean
    stat_summary(fun.y = mean, fun.ymin = sem_ymin, fun.ymax = sem_ymax)+ #Add error bars and mean dots#Put all genes in separate facets
    scale_colour_manual(values=cbPalette) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1), panel.background = element_blank(),
          panel.grid = element_blank(), panel.border = element_rect(fill=NA),aspect.ratio=0.6) +
    scale_x_discrete(expand = c(0, 0.6))+
    facet_grid(.~ Family, scales="free_y")
  return(pl)
} 


plotCandidatesDate(gen)

###========================================================================
#Expression of the putative latewood partition 2:2:1:X of the core network 
###========================================================================
   
gen <-data.frame("gen" = c("MA_131708g0010","MA_470613g0010","MA_10430216g0010",
                                                                 "MA_10973g0010",
                                                                 "MA_10432073g0010",
                                                                 "MA_112831g0010",
                                                                 "MA_90487g0010",
                                                                 "MA_379710g0010",
                                                                 "MA_371148g0010",
                                                                 "MA_10358874g0010",
                                                                 "MA_10435604g0010",
                                                                 "MA_401686g0010",
                                                                 "MA_9761124g0010",
                                                                 "MA_10350802g0010",
                                                                 "MA_10305558g0010",
                                                                 "MA_10432173g0020",
                                                                 "MA_2720g0020",
                                                                 "MA_10274939g0010",
                                                                 "MA_191789g0010",
                                                                 "MA_23596g0020",
                                                                 "MA_2720g0010",
                                                                 "MA_4044619g0010",
                                                                 "MA_4509352g0010",
                                                                 "MA_92862g0010",
                                                                 "MA_10053838g0010",
                                                                 "MA_23596g0010",
                                                                 "MA_9147880g0010",
                                                                 "MA_9198822g0010",
                                                                 "MA_10435961g0020",
                                                                 "MA_104992g0010",
                                                                 "MA_36564g0010",
                                                                 "MA_8059983g0010",
                                                                 "MA_10436313g0020",
                                                                 "MA_9337419g0010",
                                                                 "MA_10430511g0010",
                                                                 "MA_10429479g0010",
                                                                 "MA_8553148g0010",
                                                                 "MA_4719g0010",
                                                                 "MA_10435508g0010",
                                                                 "MA_68461g0010",
                                                                 "MA_6025923g0010",
                                                                 "MA_827389g0010",
                                                                 "MA_9822402g0010",
                                                                 "MA_38411g0010",
                                                                 "MA_9922406g0010",
                                                                 "MA_10435508g0020",
                                                                 "MA_517867g0010",
                                                                 "MA_10431130g0010",
                                                                 "MA_18929g0020",
                                                                 "MA_14041g0010",
                                                                 "MA_10437278g0030",
                                                                 "MA_10428184g0010",
                                                                 "MA_10432594g0020",
                                                                 "MA_1115518g0010",
                                                                 "MA_176882g0010",
                                                                 "MA_10432701g0010",
                                                                 "MA_147583g0010",
                                                                 "MA_413510g0010",
                                                                 "MA_11806g0010",
                                                                 "MA_57712g0010",
                                                                 "MA_673539g0010",
                                                                 "MA_10435107g0010",
                                                                 "MA_904g0010",
                                                                 "MA_9171580g0010",
                                                                 "MA_81147g0010",
                                                                 "MA_436575g0010",
                                                                 "MA_10435847g0010",
                                                                 "MA_10427283g0030",
                                                                 "MA_10436587g0020",
                                                                 "MA_66608g0010",
                                                                 "MA_10426545g0010",
                                                                 "MA_10426545g0020",
                                                                 "MA_10144850g0010",
                                                                 "MA_1854g0020",
                                                                 "MA_216719g0010",
                                                                 "MA_19115g0010",
                                                                 "MA_90677g0010",
                                                                 "MA_62006g0010",
                                                                 "MA_1356g0010",
                                                                 "MA_444g0010",
                                                                 "MA_224485g0010",
                                                                 "MA_10435055g0010",
                                                                 "MA_53796g0010",
                                                                 "MA_10435495g0010",
                                                                 "MA_216312g0010",
                                                                 "MA_158751g0010",
                                                                 "MA_10426903g0010",
                                                                 "MA_131887g0010",
                                                                 "MA_133989g0010",
                                                                 "MA_10437025g0010",
                                                                 "MA_6150g0020",
                                                                 "MA_782885g0010",
                                                                 "MA_496233g0010",
                                                                 "MA_180698g0010",
                                                                 "MA_10430663g0010",
                                                                 "MA_10429174g0010",
                                                                 "MA_10431010g0010",
                                                                 "MA_5534358g0010",
                                                                 "MA_9540939g0010",
                                                                 "MA_49713g0010",
                                                                 "MA_2012364g0010",
                                                                 "MA_46805g0010",
                                                                 "MA_629925g0010",
                                                                 "MA_896685g0010",
                                                                 "MA_408380g0010",
                                                                 "MA_96669g0010",
                                                                 "MA_946576g0010",
                                                                 "MA_55891g0010",
                                                                 "MA_309477g0010",
                                                                 "MA_10429296g0020",
                                                                 "MA_10433716g0010",
                                                                 "MA_10436513g0010",
                                                                 "MA_379082g0010",
                                                                 "MA_563604g0010",
                                                                 "MA_10430629g0010",
                                                                 "MA_890224g0010",
                                                                 "MA_109285g0010",
                                                                 "MA_5590688g0010",
                                                                 "MA_183426g0010",
                                                                 "MA_81001g0010",
                                                                 "MA_14242g0010",
                                                                 "MA_10428802g0010"
                                                          
                                  ),
                                          "family" = c("MA_131708g0010",rep("other",121)),
                                          stringsAsFactors=FALSE)
 
 samples <- samples[match(colnames(count.table),samples$Filename),]
   
 cbPalette <- c("deeppink3",rep("azure3",121))
 
   
 plotCandidatesDate <- function(gen){
   subdf <- as.data.frame(rlt[gen[,1],])
   colnames(subdf) <- samples$ID
   subdf <- as.data.frame(cbind(subdf,gen))
   sub.rlt <- melt(subdf) #Convert to row major data frame format
   sub.rlt[,3] <- as.character(sub.rlt[,3])
   d <- mdy(as.character(samples$Date)[match(sub.rlt[,3],samples$ID)])
   o <- order(d)
   d.dmy <- paste(day(d),month(d),year(d),sep="-")
      sub.rlt$Date <- factor(d.dmy, levels=unique(d.dmy[o]))
      sub.rlt$Tissue <- factor(samples[match(sub.rlt[,3],samples$ID),"Tissue"])
      sub.rlt$Filename <- samples[match(sub.rlt[,3],samples$ID),"Filename"]
      colnames(sub.rlt) <- c("Gene","Family","Sample","Expression","Date","Tissue","Filename") #Set column names
      sub.rlt <- sub.rlt[sub.rlt$Tissue=="Xylem",]
      sub.rlt$Family_Gene <- paste(sub.rlt$Family,sub.rlt$Gene,sep="_")
       pl <- ggplot(sub.rlt, aes(x=Date,group=Gene,y=Expression)) + #Base layer for plot
           stat_summary(fun.y = mean,geom = "line", lwd=0.7, aes(color = Family)) + #Add line for the mean
         scale_colour_manual(values=cbPalette)+ 
           theme(axis.text.x = element_text(angle = 70, hjust = 1), panel.background = element_blank(),
          panel.grid = element_blank(), panel.border = element_rect(fill=NA),aspect.ratio=0.6) 
         
         return(pl)
 }
 
 plotCandidatesDate(gen)
 
 
#' ---
#' # Aim
#' The goal is to process the graph obtained from merging 8 gene correlation
#' networking tools.
#'
#'
#' Load the libraries
library(data.table)
library(igraph)
library(LSD)
library(pander)

#' Source helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/WgcnaClusterPlot.R")


###===========================
## Cutoffs
###===========================
lenient.cutoff <- 0.89
stringent.cutoff <- 0.95
core.cutoff <- 0.99


#' # Process
#' ## Edges
#' First, we read in the edges
lenient.edges <- fread("/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/Subset/ACGGNPST-lenient-el.tsv")
stringent.edges <- fread("/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/Subset/ACGGNPST-stringent-el.tsv")
core.edges <- fread("/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/Subset/ACGGNPST-core-el.tsv")

#' And take a look at the overall edges properties
#'
#' The third column contains the ranked weights
plot(density(lenient.edges[[3]]))
plot(density(stringent.edges[[3]]))
plot(density(core.edges[[3]]))

#' The 4th column contains the number of times an edge as been reported by a
#' different gene correlation networking algorithms. Since we have a mix
#' of complete (Narromi, Anova) and sparse methods (GENIE3, TIGs, CLR), the plot
#' is as expected
barplot(table(lenient.edges[[4]]))
barplot(table(stringent.edges[[4]]))
barplot(table(core.edges[[4]]))


#' ## Graphs
#' We create an undirected graph
#'
#' for the lenient

graf <- graph.edgelist(cbind(lenient.edges[[1]],
                             lenient.edges[[2]]),
                       directed = FALSE)

# for stringent
stringent.graf <- graph.edgelist(cbind(stringent.edges[[1]],
                             stringent.edges[[2]]),
                       directed = FALSE)


#' and core thresholds
core.graf <- graph.edgelist(cbind(core.edges[[1]],
                                  core.edges[[2]]),
                            directed = FALSE)


#' And check theirs characteristics, such as the number of clusters and their sizes
#'
#' lenient graph
pander(clusters(graf)$no)
pander(clusters(graf)$csize)

#' stringent graph
pander(clusters(stringent.graf)$no)
pander(clusters(stringent.graf)$csize)



#' core graph
#'
#' The second table list the size of the clusters (first row) and the number
#' of clusters of that size (second row)
pander(clusters(core.graf)$no)
pander(table(clusters(core.graf)$csize))

#' Next we simplify the grafs (just in case we have duplicated edges)
graf <- simplify(graf)
stringent.graf <- simplify(stringent.graf)
core.graf <- simplify(core.graf)

#' ### Lenient graph characteristics
#' We check the centrality of the nodes in terms of degree and betweeness
#'
#' Degree:
g.degree <- centr_degree(graf)
plot(density(g.degree$res))
abline(v=200)
pander(sum(g.degree$res > 200))

#plot(density(g.degree$res),log="x")
plot(density(log10(g.degree$res)))
mean(g.degree$res)
median(g.degree$res)

#' Betweeness:
betweenness <- centr_betw(graf)
plot(density(betweenness$res))
plot(density(log10(betweenness$res)))

#' Let's compare their relationship. This was an attempt to see if the graph
#' could be splitted into relevant subgraph using a cutoff on the degree or
#' betweeness of the nodes. This is not the case, so we moved to looking more
#' specifically at genes of interests. This is done in the grafSubset.R script.
comparisonplot(log10(g.degree$res),log10(betweenness$res),xlim=c(0,4))
heatscatter(log10(g.degree$res),log10(betweenness$res),xlim=c(0,4))

#' Save the graph
save(graf,file="/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/Data/ACGGNPST_lenient_graph.rda")

#' And export the graph
write_graph(graf,format="graphml",
            file="/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/Data/ACGGNPST_lenient_graph.graphml")


#' ### Stringent graph characteristics
#' We check the centrality of the nodes in terms of degree and betweeness
#'
#' Degree:
g.degree <- centr_degree(stringent.graf)
plot(density(g.degree$res))
mean(g.degree$res)
median(g.degree$res)

#' Betweeness:
betweenness <- centr_betw(stringent.graf)
plot(density(betweenness$res))
plot(density(log10(betweenness$res)))

#' Let's compare their relationship.
comparisonplot(log10(g.degree$res),log10(betweenness$res),xlim=c(0,2.2))
heatscatter(log10(g.degree$res),log10(betweenness$res),xlim=c(0,2.2))

#' Save the graph
save(stringent.graf,file="/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/Data/ACGGNPST_stringent_graph.rda")

#' And export the graph
write_graph(stringent.graf,format="graphml",
            file="/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/Data/ACGGNPST_stringent_graph.graphml")

#' ### Core graph characteristics
#' We check the centrality of the nodes in terms of degree and betweeness
#'
#' Degree:
g.degree <- centr_degree(core.graf)
plot(density(g.degree$res))
mean(g.degree$res)
median(g.degree$res)

#' Betweeness:
betweenness <- centr_betw(core.graf)
plot(density(betweenness$res))
plot(density(log10(betweenness$res)))

#' Let's compare their relationship.
comparisonplot(log10(g.degree$res),log10(betweenness$res),xlim=c(0,2.2))
heatscatter(log10(g.degree$res),log10(betweenness$res),xlim=c(0,2.2))

#' Save the graph
save(core.graf,file="/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/Data/ACGGNPST_core_graph.rda")

#' And export the graph
write_graph(core.graf,format="graphml",
            file="/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/Data/ACGGNPST_core_graph.graphml")



#' Average expression of monolignol genes

gene.list <- c("MA_43667g0010",
               "MA_5898g0020",
               "MA_5898g0010",
               "MA_5904g0010",
               "MA_130482g0010",
               "MA_10436313g0010",
               "MA_9483804g0010",
               "MA_514764g0010",
               "MA_202753g0010",
               "MA_101790g0010",
               "MA_17311g0010",
               "MA_1570097g0010",
               "MA_117567g0010",
               "MA_10432784g0010",
               "MA_10433667g0010",
               "MA_480885g0020",
               "MA_480885g0010",
               "MA_91138g0020",
               "MA_825097g0010",
               "MA_10428168g0010",
               "MA_601781g0010",
               "MA_56692g0010",
               "MA_161608g0010",
               "MA_177306g0010",
               "MA_106573g0010",
               "MA_10436001g0020",
               "MA_10436001g0010",
               "MA_10436001g0030",
               "MA_8923441g0010",
               "MA_137442g0010",
               "MA_10429279g0010",
               "MA_10106144g0010",
               "MA_4911g0010",
               "MA_19223g0010",
               "MA_10429529g0010",
               "MA_10432052g0010",
               "MA_10162068g0010",
               "MA_10389075g0010",
               "MA_9392760g0010",
               "MA_7091g0010",
               "MA_125690g0020",
               "MA_10061g0010",
               "MA_6189912g0010",
               "MA_10434782g0020",
               "MA_10159692g0010",
               "MA_10430727g0010",
               "MA_10432514g0010",
               "MA_10434782g0010",
               "MA_197296g0010",
               "MA_118702g0010",
               "MA_10427334g0010",
               "MA_123220g0010",
               "MA_10427064g0010",
               "MA_138523g0010",
               "MA_8772058g0010",
               "MA_903046g0010",
               "MA_445633g0010"
)

rlto <- read.csv("analysis/HTSeq/rlogTransformed_data.csv",row.names = 1)
conditions1 <- samples$Tissue[match(colnames(rlto),samples$Filename)]
conditions2 <- samples$Sampling[match(colnames(rlto),samples$Filename)]
rltonew <- rlto[order(conditions1,conditions2)]

#Average expression values for core.graf

avg.exp <- rowMeans(rltonew)

core.graf <- set.vertex.attribute(core.graf,"average.expression",
                                  get.vertex.attribute(core.graf,"name"),
                                  avg.exp[get.vertex.attribute(core.graf,"name")])


# FOR this to work, we need a graph with the "average.expression" attribute set

WgcnaClusterPlot(matrix(get.vertex.attribute(core.graf,"average.expression")[get.vertex.attribute(core.graf,"name") %in% gene.list,],nrow=length(gene.list),byrow=TRUE))

# This is the rlt expression
WgcnaClusterPlot(as.matrix(rltonew[get.vertex.attribute(core.graf,"name"),]
                           [get.vertex.attribute(core.graf,"name") %in% gene.list,]),
                 las=2)

#' ---
#' Stringent graph neighbours selection for MA_6931g0010"
#' Set the working dir
setwd("/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/27_SpruceTimeSeries/CoExp_Networks")
#' ```
#' Load libraries
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(pander))

#' Source a helper function
source("~/Git/UPSCb/src/R/WgcnaClusterPlot.R")

#' # Process
#' ## Read the different data we need
#' Read the gene of interest
goi.filename="~/Git/UPSCb/projects/spruce-wood-time-series/doc/monolignol-genes.txt"
goi <- scan(goi.filename,
            what="character")

#' Read the graph
load(file="/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/CoExp_Networks/Data/ACGGNPST_stringent_graph.rda")

#' Read the sample info
samples <- read.csv("~/Git/UPSCb/projects/spruce-wood-time-series/doc/samplesmod.csv")


#' ## 2nd degree neighbours
#' Extract 2nd degree neighbours as subgraphs
subgraf <- make_ego_graph(stringent.graf,2,
                          get.vertex.attribute(stringent.graf,"name") %in% goi)

#' Merge the graph
second.degree.graf <- Reduce("%u%",subgraf)

#' Get some stats
#' Number of clusters
pander(clusters(second.degree.graf)$no)
#' Cluster sizes
pander(table(clusters(second.degree.graf)$csize))

#' Create an id <-> name mapping for the genes
id.mapping <- data.frame(
  id=1:length(get.vertex.attribute(second.degree.graf,"name")),
  name=get.vertex.attribute(second.degree.graf,"name"))
pander(id.mapping)

#' So that the plot remains read-able
plot(second.degree.graf,
     vertex.label=id.mapping$id,
     vertex.label.cex=1,
     vertex.size=12)


#' ## 3rd degree neighbours
#' Extract 3rd degree neighbours as subgraphs
subgraf <- make_ego_graph(stringent.graf,3,
                          get.vertex.attribute(stringent.graf,"name") %in% goi)

#' Merge the graph
third.degree.graf <- Reduce("%u%",subgraf)

#' Get some stats
#' Number of clusters
pander(clusters(third.degree.graf)$no)
#' Cluster sizes
pander(table(clusters(third.degree.graf)$csize))

#' Create an id <-> name mapping for the genes
id.mapping <- data.frame(
  id=1:length(get.vertex.attribute(third.degree.graf,"name")),
  name=get.vertex.attribute(third.degree.graf,"name"))
pander(id.mapping)

#' So that the plot remains read-able
plot(third.degree.graf,
     vertex.label=id.mapping$id,
     vertex.label.cex=1,
     vertex.size=12)


#' # SessionInfo
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```



#' ---
#' title: "Exonerate Potri missing genes"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Working directory
setwd("/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/CoExp_Networks/Infomap")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/CoExp_Networks/Infomap")
#' ```

#' Libraries
library(S4Vectors)

#' # Process
#' ## CORE

#' Read the infomap tree
core.infomap <- read.table("ACGGNPST_core_graph.tree",stringsAsFactors = FALSE)[,1:3]

#' Read the nodes
core.nodes <- scan("ACGGNPST_core_graph_nodes.txt",what="character")

#' Update the infomap
core.infomap[,3] <- core.nodes[core.infomap[,3]]

#' Split by levels
#' first level first - and ordered by the partition size: largest -> smallest
core.first.level.partitions <- split(core.infomap[,3],sub(":.*","",core.infomap[,1]))
core.first.level.partitions <- core.first.level.partitions[order(elementLengths(core.first.level.partitions),decreasing = TRUE)]
barplot(elementLengths(core.first.level.partitions))

#' ## STRINGENT
stringent.infomap <- read.table("ACGGNPST_stringent_graph.tree",stringsAsFactors = FALSE)[,1:3]
stringent.nodes <- scan("ACGGNPST_stringent_graph_nodes.txt",what="character")
stringent.infomap[,3] <- stringent.nodes[stringent.infomap[,3]]
stringent.first.level.partitions <- split(stringent.infomap[,3],sub(":.*","",stringent.infomap[,1]))
stringent.first.level.partitions <- stringent.first.level.partitions[order(elementLengths(stringent.first.level.partitions),decreasing = TRUE)]
barplot(elementLengths(stringent.first.level.partitions))

#' ## LENIENT
lenient.infomap <- read.table("ACGGNPST_lenient_graph.tree",stringsAsFactors = FALSE)[,1:3]
lenient.nodes <- scan("ACGGNPST_lenient_graph_nodes.txt",what="character")
lenient.infomap[,3] <- lenient.nodes[lenient.infomap[,3]]
lenient.first.level.partitions <- split(lenient.infomap[,3],sub(":.*","",lenient.infomap[,1]))
lenient.first.level.partitions <- lenient.first.level.partitions[order(elementLengths(lenient.first.level.partitions),decreasing = TRUE)]
barplot(elementLengths(lenient.first.level.partitions))

#' # SessionInfo
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
grep("MA_131708g0010",core.first.level.partitions)

lapply(core.first.level.partitions,grep,pattern="MA_131708g0010")

which(elementLengths(lapply(core.first.level.partitions,grep,pattern="MA_131708g0010"))>0)

grep("MA_131708g0010",core.first.level.partitions[[1]])

write.csv(core.first.level.partitions[[1]], "/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/CoExp_Networks/core-first-level-part-first-elem.csv")

core.second.level.partitions <- split(core.infomap[,"V3"],substr(core.infomap[,"V1"],1,3))

grep("MA_131708g0010",core.second.level.partitions)

which(elementLengths(lapply(core.second.level.partitions,grep,pattern="MA_131708g0010"))>0)

which(elementLengths(lapply(core.second.level.partitions,grep,pattern="MA_927631g0010"))>0)

write.csv(core.second.level.partitions[[116]], "/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/CoExp_Networks/core-second-level-part-116-elem.csv")

grep("MA_10196643g0010",core.second.level.partitions)

grep("MA_10430035g0010",core.second.level.partitions)

write.csv(core.second.level.partitions[[115]], "/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/CoExp_Networks/core-second-level-part-115-elem.csv")

write.csv(core.second.level.partitions[[415]], "/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/CoExp_Networks/core-second-level-part-415-elem.csv")

write.csv(core.infomap, "/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries/CoExp_Networks/core.infomap.soile.csv")

plot(density(elementLengths(core.second.level.partitions)))


