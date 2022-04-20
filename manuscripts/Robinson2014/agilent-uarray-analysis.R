### ============
## setup
### ============
## libs
library(Biobase)
library(limma)
library(RColorBrewer)
library(LSD)
library(scatterplot3d)
library(xtable)
library(arrayQualityMetrics)
library(sva)

## wdir
setwd("/mnt/picea/storage/projects/03_Aspen_Project")

## source some helper
source("~/Git/UPSCb/src/R/plot.multidensity.R")

## read the annot
annot <- read.csv("~/Git/UPSCb/projects/sex/doc/agilent-uarray.csv",stringsAsFactors=FALSE)
annot$Sample <- sub("\\.0$","",as.character(annot$Sample))

## Make a vector of file names
fnames <- dir(path="uarray/gpr/",pattern=".gpr$",full.names=TRUE)

## get the sample names
smpl.sel <- match(paste(annot$Barcode,sub("1_","",annot$Block),sep="_"),sub("-4\\.gpr","",basename(fnames)))
sample <- annot[smpl.sel,"Sample"]

##Define some pretty colours
cols <-brewer.pal(8, "Dark2")

### ============
## read
### ============
##Define which foreground and background values to use
Cy3 <- "F532 Median"
Cy3b <- "B532 Median"

##Read in the gpr files and only the green channel
rawData <- read.maimages(fnames,
                    source="genepix",
                    columns=list(R=Cy3,Rb=Cy3b),green.only=TRUE)

## check the boxplot of every sample
png(filename="uarray/analysis/raw-boxplot.png",
    width=800, height=600, pointsize=16)
boxplot(log2(rawData$E+1),col=rep(cols[1:4],5),
        main="Raw intensities per array block",
        names=sample,las=2,ylab="log2(I+1)")
dev.off()

rawList <- lapply(1:length(sample),function(i){log2(rawData$E[,i]+1)})
names(rawList) <- sample
png(filename="uarray/analysis/raw-density.png",
    width=800, height=600, pointsize=16)
plot.multidensity(rawList,xlab="raw intensity",legend.x="topright",legend.cex=0.8)
dev.off()

##Spatial plots of raw foreground and background signals
sapply(1:length(sample),function(k){
  png(filename=paste("uarray/analysis/spatial_plot", 
                     sample[k], "png", sep="."), 
      width=800, height=600, pointsize=16)
  par(mfcol=c(1,2))
  imageplot(rawData$E[,k], rawData$printer, 
            zerocenter=TRUE, low="white",
            high="green", main="Green FG", zlim=c(500,40000), mar=c(2,2,3,2))
  imageplot(rawData$Eb[,k],rawData$printer, 
            zerocenter=TRUE, low="white", 
            high="green", main="Green BG", zlim=c(50,5000), mar=c(2,2,3,2))
  dev.off()
})

### ============
## preprocess
### ============
## background correct
bgData <- backgroundCorrect(rawData,method="normexp")
png(filename="uarray/analysis/bg-corrected-boxplot.png",
    width=800, height=600, pointsize=16)
boxplot(log2(bgData$E),col=rep(cols[1:4],5),
        main="BG corrected intensities per array block",
        names=sample,las=2,ylab="log2(I+1)")
dev.off()

bgList <- lapply(1:length(sample),function(i){log2(bgData$E[,i])})
names(bgList) <- sample
png(filename="uarray/analysis/bg-corrected-density.png",
    width=800, height=600, pointsize=16)
plot.multidensity(bgList,xlab="bg corrected intensity",legend.x="topright",legend.cex=0.8)
dev.off()

## normalize between array
nData <- normalizeBetweenArrays(bgData,method="quantile")
png(filename="uarray/analysis/normalized-boxplot.png",
    width=800, height=600, pointsize=16)
boxplot(nData$E,col=rep(cols[1:4],5),
        main="Btw-array normalized intensities per array block",
        names=sample,las=2,ylab="log2(I+1)")
dev.off()

nList <- lapply(1:length(sample),function(i){nData$E[,i]})
names(nList) <- sample
png(filename="uarray/analysis/normalized-density.png",
    width=800, height=600, pointsize=16)
plot.multidensity(nList,xlab="normalized intensity",legend.x="topright",legend.cex=0.8)
dev.off()

### cyclicloess show no significant difference
# nData <- normalizeBetweenArrays(bgData,method="cyclicloess")
# png(filename="uarray/analysis/cyclicloess-normalized-boxplot.png",
#     width=800, height=600, pointsize=16)
# boxplot(nData$E,col=rep(cols[1:4],5),
#         main="Btw-array normalized intensities per array block",
#         names=sample,las=2,ylab="log2(I+1)")
# dev.off()

## arrayQualityMetrics
arrayQualityMetrics(ExpressionSet(assayData=nData$E),
                    outdir="uarray/analysis/qa",
                    reporttitle="Agilent by Sex")

## as we had seen from the previous analysis, sample 17 is inexplicably different
## discard 17 and redo the btw array norm
bgData$targets <- bgData$targets[-17]
bgData$E <- bgData$E[,-17]
sample <- sample[-17]
smpl.sel <- smpl.sel[-17]

nData <- normalizeBetweenArrays(bgData,method="quantile")
png(filename="uarray/analysis/normalized-boxplot-19-sample.png",
    width=800, height=600, pointsize=16)
boxplot(nData$E,col=rep(cols[1:4],5)[-17],
        main="Btw-array normalized intensities per array block",
        names=sample,las=2,ylab="log2(I+1)")
dev.off()

### ============
## QA
### ============
## arrayQualityMetrics
arrayQualityMetrics(ExpressionSet(assayData=nData$E),
                    outdir="uarray/analysis/qa-19-samples",
                    reporttitle="Agilent 19 sex samples")

##MA plot of average male v's female signal
male<-which(annot[smpl.sel,"Sex"]=="M")
female<-which(annot[smpl.sel,"Sex"]=="F")
Male<-rowMeans(nData$E[,male])
Female<-rowMeans(nData$E[,female])
png(filename="uarray/analysis/MA-19-samples-plot.png", width=800, height=600, pointsize=12, bg="white")
plot(.5*(Male+Female),Male-Female, xlab="A", ylab="M",pch=19,cex=.8)
abline(0,0, col="grey", lwd=3)
abline(mean(Male - Female),0, col="blue", lwd=3)
abline(1,0, col="red", lwd=3)
abline(-1,0, col="green", lwd=3)
abline(v=8, col="grey", lwd=3)
dev.off()

png(filename="uarray/analysis/MA-heatscatter-19-samples-plot.png", width=800, height=600, pointsize=12, bg="white")
heatscatter(.5*(Male+Female),Male-Female, xlab="A", ylab="M")
abline(h=0, col=1, lwd=2,lty=3)
abline(mean(Male-Female),0, col="blue", lwd=2, lty=3)
abline(1,0, col="red", lwd=2, lty=3)
abline(-1,0, col="green", lwd=2, lty=3)
abline(v=7, col="grey", lwd=2, lty=3)
legend("topright",col=c(1,4,2,3),lty=3,legend=c("baseline","mean difference","1 log2FC","-1 log2FC"))
dev.off()

### ============
## filter
### ============
## filter low expression
## Now filter out control probes and low expressed probes. 
## To get an idea of how bright expression probes should be, 
## we compute the 95% percentile of the negative control probes 
## on each array. We keep probes that are at least 10% brighter
## than the negative controls on at least four arrays 
## (because there are four replicates):
## the negative probes are most probably called "DarkCorner"

## 1. check there are some NA, but only 2 for 606 probes
colSums(is.na(nData$E[nData$genes$ID=="DarkCorner",]))

## 2. then filter
neg95 <- apply(nData$E[nData$genes$ID=="DarkCorner",],
               2,function(x){quantile(x,p=0.95,na.rm=TRUE)})
cutoff <- matrix(1.1*neg95,nrow(nData),ncol(nData),byrow=TRUE)
isexpr <- rowSums(nData$E > cutoff) >= 4
table(isexpr)

## Regular probes have names starting with "pt" in the ID column:
fData <- nData[grepl("pt_",nData$genes$ID)  & isexpr,]

png(filename="uarray/analysis/filtered-boxplot.png",
    width=800, height=600, pointsize=16)
boxplot(fData$E,col=rep(cols[1:4],5)[-17],
        main="filtered normalized intensities per array block",
        names=sample,las=2,ylab="expression value")
dev.off()

fList <- lapply(1:length(sample),function(i){fData$E[,i]})
names(fList) <- sample
png(filename="uarray/analysis/filtered-density.png",
    width=800, height=600, pointsize=16)
plot.multidensity(fList,xlab="expression",legend.x="topright",legend.cex=0.8,
                  col=rep(cols[1:8],length.out=19)[-17],lty=rep(1:2,length.out=19))
dev.off()

### ============
## Annot
### ============
## read in the array annotation information
anno<-read.table("uarray/annotation/AgilentAnnotation121510.txt", sep="\t", 
                 header=TRUE, fill=TRUE, 
                 quote="",stringsAsFactors=FALSE)
## str(anno)

## match probe IDs and add gene model to RG.final$genes$Name
## using the genomic location instead
fData$genes$Name<-anno$Genomic[match(fData$genes$ID, anno$Probe.Id)]

## some are unknown
sel <- fData$genes$Name == "" 
sum(sel)
fData$genes$Name[sel]<-sprintf("Unknown%04d",1:sum(sel))

## 589 probes match the same gene
## and that's not caring about the multi-gene lines, those comma-separated
sum(duplicated(fData$genes$Name))
## with a few occuring several times
table(table(fData$genes$Name[duplicated(fData$genes$Name)]))

## save the data
write.csv(cbind(ID=fData$genes$Name,fData$E),
          file="uarray/analysis/agilent-microarray_normalized-data.csv",
          row.names=FALSE)

### ============
## PCA
### ============
pca <- prcomp(t(fData$E))
percent <- round(summary(pca)$importance[2,]*100)
sex.cols <- c("pink","lightblue")

pdf(file="uarray/analysis/PCA-by-Sex-19-samples.pdf", width=10, height=10)
scatterplot3d(pca$x[,1],pca$x[,2],pca$x[,3],
              xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3(",percent[3],"%)",sep=""),
              color=sex.cols[as.integer(pData$sex)],
              pch=19,main="Principal Components Analysis")
legend("topleft",pch=19,col=sex.cols,
       legend=levels(factor(pData$sex)),bty="n")
dev.off()

png(filename="uarray/analysis/PCA-by-Sex-19-samples.png", width=800, height=600, pointsize=12, bg="white")
scatterplot3d(pca$x[,1],pca$x[,2],pca$x[,3],
              xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3(",percent[3],"%)",sep=""),
              color=sex.cols[as.integer(pData$sex)],
              pch=19,main="Principal Components Analysis")
legend("topleft",pch=19,col=sex.cols,
       legend=levels(factor(pData$sex)),bty="n")
dev.off()

## Color by barcode
pdf(file="uarray/analysis/PCA-by-Barcode-19-samples.pdf", width=10, height=10)
scatterplot3d(pca$x[,1],pca$x[,2],pca$x[,3],
              xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3(",percent[3],"%)",sep=""),
              color=cols[as.integer(factor(annot$Barcode[-17]))],
              pch=19,main="Principal Components Analysis")
legend("topleft",pch=c(NA,rep(19,5)),col=c(NA,cols[1:5]),
       legend=c("Barcode",levels(factor(annot$Barcode))),bty="n")
dev.off()

png(filename="uarray/analysis/PCA-by-Barcode-19-samples.png", width=800, height=600, pointsize=12, bg="white")
scatterplot3d(pca$x[,1],pca$x[,2],pca$x[,3],
              xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3(",percent[3],"%)",sep=""),
              color=cols[as.integer(factor(annot$Barcode))][-17],
              pch=19,main="Principal Components Analysis")
legend("topleft",pch=c(NA,rep(19,5)),col=c(NA,cols[1:5]),
       legend=c("Barcode",levels(factor(annot$Barcode))),bty="n")
dev.off()

## Color by block
pdf(file="uarray/analysis/PCA-by-Block-19-samples.pdf", width=10, height=10)
scatterplot3d(pca$x[,1],pca$x[,2],pca$x[,3],
              xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3(",percent[3],"%)",sep=""),
              color=cols[as.integer(factor(annot$Block))][-17],
              pch=19,main="Principal Components Analysis")
legend("topleft",pch=c(NA,rep(19,4)),col=c(NA,cols[1:4]),
       legend=c("Block",levels(factor(annot$Block))),bty="n")
dev.off()

png(filename="uarray/analysis/PCA-by-Block-19-samples.png", width=800, height=600, pointsize=12, bg="white")
scatterplot3d(pca$x[,1],pca$x[,2],pca$x[,3],
              xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3(",percent[3],"%)",sep=""),
              color=cols[as.integer(factor(annot$Block))][-17],
              pch=19,main="Principal Components Analysis")
legend("topleft",pch=c(NA,rep(19,4)),col=c(NA,cols[1:4]),
       legend=c("Block",levels(factor(annot$Block))),bty="n")
dev.off()

## 2D
png(filename="uarray/analysis/PCA-Sex_1D-2D-19-samples.png", width=800, height=600, pointsize=12, bg="white")
plot(pca$x[,1],pca$x[,2],
     xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
     col=sex.cols[as.integer(pData$sex)],pch=19,
     main="Principal Components Analysis")
legend("bottomright",pch=19,
       col=sex.cols,
       legend=levels(factor(pData$sex)))
dev.off()

## Added Sept 23rd.
## data has been mixed up with limma-for-knitr.R script
## i.e. 
## load("Robinson-2013-microarray.rda")
## pData <- as(sample.info[-17,],"AnnotatedDataFrame")
## setwd("/mnt/picea/storage/projects/03_Aspen_Project/uarray/analysis")
png(file="Fig3A.png",height=1200,width=1200,res=300,pointsize=6)
plot(pca$x[,1],pca$x[,2],
     col=sex.cols[as.integer(factor(pData$Sex))],pch=19,
     xlab="",ylab="",main="",axes=FALSE,frame.plot=TRUE)
axis(side=1,at=seq(-40,40,20),labels=rep("",5))
axis(side=2,at=seq(-80,20,20),labels=rep("",6))
dev.off()

pdf(file="uarray/analysis/PCA-Sex_1D-2D-19-samples.pdf", width=10, height=10)
plot(pca$x[,1],pca$x[,2],
     xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
     col=sex.cols[as.integer(pData$sex)],pch=19,
     main="Principal Components Analysis")
legend("bottomright",pch=19,
       col=sex.cols,
       legend=levels(factor(pData$sex)))
dev.off()

## 1D-2D by Barcode
png(filename="uarray/analysis/PCA-Barcode_1D-2D-19-samples.png", width=800, height=600, pointsize=12, bg="white")
plot(pca$x[,1],pca$x[,2],
     xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
     col=cols[as.integer(factor(annot$Barcode))][-17],pch=19,
     main="Principal Components Analysis")
legend("bottomright",pch=c(NA,rep(19,5)),col=c(NA,cols[1:5]),
       legend=c("Barcode",levels(factor(annot$Barcode))),bty="n")
dev.off()

pdf(file="uarray/analysis/PCA-Barcode_1D-2D-19-samples.pdf", width=10, height=10)
plot(pca$x[,1],pca$x[,2],
     xlab=paste("Comp. 1(",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2(",percent[2],"%)",sep=""),
     col=cols[as.integer(factor(annot$Barcode))][-17],pch=19,
     main="Principal Components Analysis")
legend("bottomright",pch=c(NA,rep(19,5)),col=c(NA,cols[1:5]),
       legend=c("Barcode",levels(factor(annot$Barcode))),bty="n")
dev.off()

png(filename="uarray/analysis/PCA-Block_2D-3D-19-samples..png", width=800, height=600, pointsize=12, bg="white")
plot(pca$x[,2],pca$x[,3],
     xlab=paste("Comp. 2(",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3(",percent[3],"%)",sep=""),
     col=cols[as.integer(factor(annot$Block))][-17],pch=19,
     main="Principal Components Analysis")
legend("bottomright",pch=c(NA,rep(19,4)),col=c(NA,cols[1:4]),
       legend=c("Block",levels(factor(annot$Block))),bty="n")
dev.off()

pdf(file="uarray/analysis/PCA-Block_2D-3D-19-samples..pdf", width=10, height=10)
plot(pca$x[,2],pca$x[,3],
     xlab=paste("Comp. 2(",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3(",percent[3],"%)",sep=""),
     col=cols[as.integer(factor(annot$Block))][-17],pch=19,
     main="Principal Components Analysis")
legend("bottomright",pch=c(NA,rep(19,4)),col=c(NA,cols[1:4]),
       legend=c("Block",levels(factor(annot$Block))),bty="n")
dev.off()

### ============
## correct batch
### ============
colnames(fData$E) <- sample
pData <- as(annot[smpl.sel,],"AnnotatedDataFrame")
sampleNames(pData) <- sample
eSet <- ExpressionSet(assayData=fData$E,
                      phenoData=pData,
                      featureData=as(fData$genes,"AnnotatedDataFrame"))

## ComBat / sva won't work as the system is exactly singular
mod = model.matrix(~as.factor(Sex), data=pData(pData))
mod = model.matrix(~Block+Sex, data=pData(pData))
mod0 = model.matrix(~1,data=pData(pData))
n.sv = num.sv(exprs(eSet),mod,method="leek")
n.sv
eSetB <- ComBat(dat=exprs(eSet),batch=as.character(pData$Barcode),mod=design)


### ============
## limma
### ============
## pheno data
f<-factor(pData$Sex, levels=unique(pData$Sex))

mod = model.matrix(~as.character(Barcode)+Block+Sex, data=pData(pData))
fitBBS <- lmFit(fData, mod)
fitBBS <- eBayes(fitBBS)
topTable(fitBBS,coef="SexM",p.value=0.01)

datBBS <- topTable(fitBBS,coef="SexM",n=nrow(fitBBS))
write.csv(datBBS,file="uarray/analysis/agilent-microarray_limma-array-block-blocked.csv",row.names=FALSE,quote=FALSE)


mod = model.matrix(~as.character(Barcode)+Sex, data=pData(pData))
fitBS <- lmFit(fData, mod)
fitBS <- eBayes(fitBS)
topTable(fitBS,coef="SexM",p.value=0.01)

mod = model.matrix(~Sex, data=pData(pData))
fit <- lmFit(fData, mod)
fit <- eBayes(fit)
topTable(fit,coef="SexM",p.value=0.01)

## is equivalent to that
# design = model.matrix(~0+Sex, data=pData(pData))
# colnames(design) <- levels(f) 
# contrast.matrix <- makeContrasts("M-F", levels=design)
# fit <- lmFit(fData, design)
# fit <- contrasts.fit(fit,contrast.matrix)
# fit <- eBayes(fit)
# topTable(fit,p.value=0.01)

### ============
## volcano
### ============
png(filename="uarray/analysis/male-vs-female_volcanoplot.png", 
    width=800, height=600, pointsize=16)
heatscatter(fit$coefficients[,"SexM"],fit$lods[,"SexM"],
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="Log Odds")
dev.off()

pdf(file="uarray/analysis/male-vs-female_volcanoplot.pdf", 
    width=10, height=10)
heatscatter(fit$coefficients[,"SexM"],fit$lods[,"SexM"],
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="Log Odds")
dev.off()

png(filename="uarray/analysis/male-vs-female_batch-blocked_volcanoplot.png", 
    width=800, height=600, pointsize=16)
heatscatter(fitBS$coefficients[,"SexM"],fitBS$lods[,"SexM"],
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="Log Odds",
            sub="Batch effect blocked",cex.sub=0.8)
dev.off()

pdf(file="uarray/analysis/male-vs-female_batch-blocked_volcanoplot.pdf", 
    width=10, height=10)
heatscatter(fitBS$coefficients[,"SexM"],fitBS$lods[,"SexM"],
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="Log Odds",
            sub="Batch effect blocked",cex.sub=0.8)
dev.off()

png(filename="uarray/analysis/male-vs-female_batch-block-blocked_volcanoplot.png", 
    width=800, height=600, pointsize=16)
heatscatter(fitBBS$coefficients[,"SexM"],
            fitBBS$lods[,"SexM"],
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="Log Odds",
            sub="Batch and Block effect blocked",cex.sub=0.8)
dev.off()

pdf(file="uarray/analysis/male-vs-female_batch-block-blocked_volcanoplot.pdf", 
    width=10, height=10)
heatscatter(fitBBS$coefficients[,"SexM"],fitBBS$lods[,"SexM"],
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="Log Odds",
            sub="Batch and Block effect blocked",cex.sub=0.8)
dev.off()

png(filename="uarray/analysis/male-vs-female_batch-block-blocked_p-value-volcanoplot.png", 
    width=800, height=600, pointsize=16)
heatscatter(fitBBS$coefficients[,"SexM"],
            -log10(fitBBS$p.value[,"SexM"]),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) p-value",
            sub="Batch and Block effect blocked",cex.sub=0.8)
abline(h=2,lty=2,col="gray")
legend("topleft",bty="n","1% p-value cutoff",lty=2,col="gray")
dev.off()

pdf(file="uarray/analysis/male-vs-female_batch-block-blocked_p-value-volcanoplot.pdf", 
    width=10, height=10)
heatscatter(fitBBS$coefficients[,"SexM"],
            -log10(fitBBS$p.value[,"SexM"]),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) p-value",
            sub="Batch and Block effect blocked",cex.sub=0.8)
abline(h=2,lty=2,col="gray")
legend("topleft",bty="n","1% adj p-value cutoff",lty=2,col="gray")
dev.off()

png(filename="uarray/analysis/male-vs-female_batch-block-blocked_-log10-p-value-volcanoplot.png", 
    width=800, height=600, pointsize=16)
heatscatter(fitBBS$coefficients[,"SexM"],
            -log10(1- (exp(fitBBS$lods[,"SexM"]) / (1+exp(fitBBS$lods[,"SexM"])))),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) adj. p-value",
            sub="Batch and Block effect blocked",cex.sub=0.8)
dev.off()

pdf(file="uarray/analysis/male-vs-female_batch-block-blocked_-log10-p-value-volcanoplot.pdf", 
    width=10, height=10)
heatscatter(fitBBS$coefficients[,"SexM"],
            -log10(1- (exp(fitBBS$lods[,"SexM"]) / (1+exp(fitBBS$lods[,"SexM"])))),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) adj. p-value",
            sub="Batch and Block effect blocked",cex.sub=0.8)
dev.off()

## Added  Sept. 23
## see previous comment added on the same date
png(file="Fig3B.png",height=1200,width=1200,res=300,pointsize=6)
heatscatter(fitBBS$coefficients[,"SexM"],
            -log10(1- (exp(fitBBS$lods[,"SexM"]) / (1+exp(fitBBS$lods[,"SexM"])))),
            main="",cor=FALSE,
            xlab="", ylab="",axes=FALSE,frame.plot=TRUE)
axis(side=1,at=seq(-3,2,1),label=rep("",6))
axis(side=2,at=seq(.004,.01,.002),label=rep("",4))
dev.off()

## final save
save.image("uarray/analysis/analysis_20130918.rda")
