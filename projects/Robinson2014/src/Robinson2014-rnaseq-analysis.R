### ==============================
## load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(arrayQualityMetrics))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(gplots))

source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/densityPlot.R")
source("~/Git/UPSCb/src/R/plotMA.R")
source("~/Git/UPSCb/src/R/volcanoPlot.R")

### ==============================
## set the working dir
### ==============================
# setwd("/mnt/picea/projects/aspseq/sex")
setwd("/Users/delhomme/Desktop/tmp/sex")

### ==============================
## read the HTSeq files in a matrix
### ==============================
# res <- mclapply(dir("HTSeq/Ptrichocarpa",pattern="^[2,3].*_STAR\\.txt",full.names=TRUE),function(fil){
res <- mclapply(dir(".",pattern="^[2,3].*_STAR\\.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=2)
# names(res) <- gsub("_.*_STAR\\.txt","",dir("HTSeq/Ptrichocarpa",pattern="^[2,3].*_STAR\\.txt"))
names(res) <- gsub("_.*_STAR\\.txt","",dir(".",pattern="^[2,3].*_STAR\\.txt"))

### ==============================
## get the count table
### ==============================
addInfo <- c("__no_feature","__ambiguous",
             "__too_low_aQual","__not_aligned",
             "__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]
write.csv(count.table,"analysis/sex-raw-unormalised-data.csv")

### ==============================
## extract the HTSeq stat lines
### ==============================
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

## as percentages
apply(count.stats,2,function(co){round(co*100/sum(co))})

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

### ==============================
## 11.4% of the genes are not expressed
## out of a total of 41335 genes
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
plot.multidensity(log10(count.table),col=rep(pal,each=3),
                  legend.x="topright",legend.cex=0.5,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

### ==============================
## For visualization, the data is
## submitted to a variance stabilization
## transformation using DESeq2. The 
## dispersion is estimated independently
## of the sample type
### ==============================
samples <- read.csv("~/Git/UPSCb/projects/sex/doc/sex-samples.csv")
conditions <- names(res)
sex <- samples$sex[match(names(res),samples$sample)]

## correct the sex based on sex determination gene
sex[names(res) == "226.1"] <- "M"

## get the date
date <- factor(samples$date[match(names(res),samples$sample)])

## create the dds object
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       sex=sex,
                       date=date),
  design = ~ condition)

### ==============================
## check the size factors
## no big variation, can go for vsd
## over rld
### ==============================
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- names(res)
sizes

## do the vsd
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)
vst <- vst - min(vst)
write.csv(vst,"analysis/sex-library-size-normalized_variance-stabilized_data.csv")

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
meanSdPlot(assay(vsd)[rowSums(count.table)>0,], ylim = c(0,2.5))

### ==============================
## perform a Principal Component
## Analysis (PCA) of the data
## to do a quick quality assessment
## i.e. replicate should cluster
## and the first dimensions should 
## be explainable by biological means.
### ==============================
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
smpls <- conditions

## color by date and sex
sex.cols<-c("pink","lightblue")
sex.names<-c(F="Female",M="Male")

## is the separation sex?
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=sex.cols[as.integer(factor(sex))],
              pch=19)
legend("bottomright",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
       legend=c("Color:",sex.names[levels(factor(sex))]))
par(mar=mar)

## or date?
dat <- date
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(factor(dat))],
              pch=19)
legend("topleft",pch=rep(c(19,23),each=10),col=rep(pal,2),legend=levels(factor(dat)),bty="n")
par(mar=mar)

### final plot
### color = sex and sample = date
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=sex.cols[as.integer(factor(sex))],
              pch=c(19,17)[as.integer(factor(dat))])
legend("bottomright",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
       legend=c("Color:",sex.names[levels(factor(sex))]))
legend("topleft",pch=c(NA,21,24),col=c(NA,1,1),
       legend=c("Symbol:",levels(factor(dat))),cex=0.85)
par(mar=mar)

### ==============================
## and the first two dims
### ==============================
plot(pc$x[,1],  
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=sex.cols[as.integer(factor(sex))],
     pch=c(19,17)[as.integer(factor(dat))],
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottomleft",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
       legend=c("Color:",sex.names[levels(factor(sex))]))
legend("topright",pch=c(NA,21,24),col=c(NA,1,1),
       legend=c("Symbol:",levels(factor(dat))),cex=0.85)
text(pc$x[,1],  
     pc$x[,2],
     labels=conditions,cex=.5,adj=-1)


### ==============================
## and the 2nd and 3rd dims
### ==============================
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=sex.cols[as.integer(factor(sex))],
     pch=c(19,17)[as.integer(factor(dat))],     
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottomleft",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
       legend=c("Color:",sex.names[levels(factor(sex))]))
legend("topright",pch=c(NA,21,24),col=c(NA,1,1),
       legend=c("Symbol:",levels(factor(dat))),cex=0.85)

### ==============================
## plot 
### ==============================
## do the rld
rld <- rlogTransformation(dds)
rlt <- assay(rld)
colnames(rlt) <- colnames(count.table)
write.csv(rlt,"analysis/sex-library-size-normalized_regularized-log-transformed_data.csv")

## get the PCA from RLT data
rpc <- prcomp(t(rlt))
percent <- round(summary(rpc)$importance[2,]*100)

## and plot the first 2 dim
plot(rpc$x[,1],
     rpc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=sex.cols[as.integer(factor(sex))],
     pch=c(19,17)[as.integer(factor(dat))],
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottomleft",pch=c(NA,15,15),col=c(NA,sex.cols[1:2]),
       legend=c("Color:",sex.names[levels(factor(sex))]))
legend("topright",pch=c(NA,21,24),col=c(NA,1,1),
       legend=c("Symbol:",levels(factor(dat))),cex=0.85)
text(pc$x[,1],  
     pc$x[,2],
     labels=conditions,cex=.5,adj=1)
par(mar=mar)

### ==============================
## heatmap and hc (samples and samples x genes)
### ==============================
#heatmap.2(vst, col=redgreen(75), scale="row", ColSideColors=patientcolors,
#          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)

### ==============================
## arrayQualityMetrics
### ==============================
#suppressMessages(arrayQualityMetrics(ExpressionSet(assayData=vst),
#                                     outdir="analysis/aQM"))

### ==============================
## do the Diff. Exp. -- with corrected sex
### ==============================
design(dds) <- ~date+sex
dds <- DESeq(dds)

swap <- counts(dds,normalized=TRUE)["Potri.019G047300.0",]
names(swap) <- sex

plotDispEsts(dds)
res <- results(dds,contrast=c("sex","F","M"))

res["Potri.019G047300.0",]

alpha=0.01
plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
hist(res$padj,breaks=seq(0,1,.01))
sum(res$padj<alpha,na.rm=TRUE)
UmAspG <- rownames(res[order(res$padj,na.last=TRUE),][1:sum(res$padj<alpha,na.rm=TRUE),])

### ==============================
## do the Diff. Exp. -- with the original sample
### ==============================
sex[conditions=="226.1"] <- "F"
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       sex=sex,
                       date=date),
  design = ~ date+sex)

dds <- DESeq(dds)
plotDispEsts(dds)
res.orig <- results(dds,contrast=c("sex","F","M"))
plotMA(res.orig,alpha=alpha)
volcanoPlot(res.orig,alpha=alpha)
hist(res.orig$padj,breaks=seq(0,1,.01))
sum(res.orig$padj<alpha,na.rm=TRUE)
UmAspG.orig <- rownames(res.orig[order(res.orig$padj,na.last=TRUE),][1:sum(res.orig$padj<alpha,na.rm=TRUE),])

## the same using the LRT, more similar to DESeq
## even more sensitive to an outlier
# dds <- DESeq(dds,test="LRT",reduced= ~date)
# res.lrt <- results(dds,name="sex_M_vs_F")
# plotMA(res.lrt,alpha=alpha)
# volcanoPlot(res.lrt,alpha=alpha)
# hist(res.lrt$padj,breaks=seq(0,1,.01))
# sum(res.lrt$padj<alpha,na.rm=TRUE)
# UmAspG.lrt <- rownames(res.orig[order(res.lrt$padj,na.last=TRUE),][1:sum(res.lrt$padj<alpha,na.rm=TRUE),])

### ==============================
## do the Diff. Exp. -- without the sample
### ==============================
sel <- conditions != "226.1"
dds <- DESeqDataSetFromMatrix(
  countData = count.table[,sel],
  colData = data.frame(condition=conditions[sel],
                       sex=sex[sel],
                       date=date[sel]),
  design = ~ date+sex)

dds <- DESeq(dds)
plotDispEsts(dds)
res.sel <- results(dds,contrast=c("sex","F","M"))
plotMA(res.sel,alpha=alpha)
volcanoPlot(res.sel,alpha=alpha)
hist(res.sel$padj,breaks=seq(0,1,.01))
sum(res.sel$padj<alpha,na.rm=TRUE)
UmAspG.sel <- rownames(res.sel[order(res.sel$padj,na.last=TRUE),][1:sum(res.sel$padj<alpha,na.rm=TRUE),])

### ==============================
## get the SwAsp results 
## ==============================
#SwAsp <- read.csv("analysis/DESeq2/mock_sex_deseq2_result.csv",row.names=1)
SwAsp <- read.csv("mock_sex_deseq2_result.csv",row.names=1)
plotMA(DataFrame(SwAsp),alpha=alpha)
volcanoPlot(DataFrame(SwAsp),alpha=alpha)
SwAspG <- rownames(SwAsp[order(SwAsp$padj,na.last=TRUE),][1:sum(SwAsp$padj<alpha,na.rm=TRUE),])

### ==============================
## run DESeq
### ==============================
suppressPackageStartupMessages(library(DESeq))
yearSexDesign <- samples[grep("[2,3].*",samples$sample),c("sex","date")]

## we want to look at the sex by blocking the year effect
cdsFull <- newCountDataSet(count.table,yearSexDesign)
cdsFull = estimateSizeFactors( cdsFull )
cdsFull = estimateDispersions( cdsFull )

## create both models
fit1 = suppressWarnings(fitNbinomGLMs( cdsFull, count ~ date + sex ))
fit0 = suppressWarnings(fitNbinomGLMs( cdsFull, count ~ date))

## keep only the converged ones
sel <- !is.na(fit1$converged) & !is.na(fit0$converged) & fit1$converged & fit0$converged
fit1 <- fit1[sel,]
fit0 <- fit0[sel,]

## calculate the GLM p-value and adjust them
pvalsGLM = nbinomGLMTest( fit1, fit0 )
padjGLM = p.adjust( pvalsGLM, method="BH" )

## filter those with a way too high FC
## they are due to genes with very few read in very few samples
boxplot(rowSums(count.table[rownames(fit1[fit1$sexM > 20,]),]>0),
        ylab="number of sample with expression > 0",
        main="distribution of the # of samples with\nexpression > 0 for genes with extreme high log2FC"
)
boxplot(rowSums(count.table[rownames(fit1[fit1$sexM < -20,]),]>0),
        ylab="number of sample with expression > 0",
        main="distribution of the # of samples with\nexpression > 0 for genes with extreme low log2FC"
)
plot(density(padjGLM[fit1$sexM > 20]),
     main="Adj. p-value for genes with extreme log2FC",
     xlab="Adj. p-value")
plot(density(padjGLM[fit1$sexM > -20]),
     main="Adj. p-value for genes with extreme log2FC",
     xlab="Adj. p-value")

sel <- fit1$sexM > -20 & fit1$sexM < 20
fit1 <- fit1[sel,]
padjGLM <- padjGLM[sel]
pvalsGLM <- pvalsGLM[sel]

## one p-value is equal to 0
## so change it to be 1 log10 fold change smaller than
## the smallest one for plotting
padjGLM[padjGLM==0] <- min(padjGLM[padjGLM!=0])/10
pvalsGLM[pvalsGLM==0] <- min(pvalsGLM[pvalsGLM!=0])/10

## volcano plot with the non-adjusted p-values
heatscatter(fit1$sexM,-log10(pvalsGLM),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) p-value")
legend("topleft",bty="n","1% p-value cutoff",lty=2,col="gray")
points(fit1[pvalsGLM<.01,"sexM"],-log10(pvalsGLM[pvalsGLM<.01]),col="lightblue",pch=19)
## circle the points for the dot plot
pos <- c(order(pvalsGLM)[1:4],order(fit1$sexM,decreasing=TRUE)[1:4],order(fit1$sexM)[1:4])
points(fit1[pos,"sexM"],-log10(pvalsGLM[pos]),col="red")
abline(h=2,lty=2,col="gray")

## volcano plot with adjusted p-values
heatscatter(fit1$sexM,-log10(padjGLM),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) adj. p-value")
legend("topleft",bty="n","10% FDR cutoff",lty=2,col="gray")
points(fit1[padjGLM<.1,"sexM"],-log10(padjGLM[padjGLM<.1]),col="lightblue",pch=19)
## circle the points for the dot plot
pos <- c(order(padjGLM)[1:4],order(fit1$sexM,decreasing=TRUE)[1:4],order(fit1$sexM)[1:4])
points(fit1[pos,"sexM"],-log10(padjGLM[pos]),col="red")
abline(h=1,lty=2,col="gray")

## the genes
UmAspD <- rownames(fit1[padjGLM<.01,])

### ==============================
## compare all
### ==============================
plot.new()
grid.draw(venn.diagram(list(
  UmAsp=UmAspG,
  UmAsp.orig=UmAspG.orig,
  UmAsp.sel=UmAspG.sel,
  UmAspD = UmAspD,
  SwAsp=SwAspG),
                       filename=NULL,
                       col=pal[1:5],
                       category.names=c("UmAsp - corrected",
                                        "UmAsp - original",
                                        "UmAsp - removed",
                                        "UmAsp - DESeq",
                                        "SwAsp")))

sort(intersect(UmAspG,SwAspG))
sort(intersect(UmAspG.orig,SwAspG))
sort(intersect(UmAspG.orig,UmAspG))
sort(intersect(UmAspG.sel,UmAspG))
sort(intersect(UmAspG.orig,UmAspG.sel))

### ==============================
## For Mike to explain why DESeq2
## is so sensitive to misclassification
### ==============================

## mis-classified
sex <- samples$sex[match(colnames(count.table),samples$sample)]
date <- factor(samples$date[match(colnames(count.table),samples$sample)])
ddsM <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       sex=sex,
                       date=date),
  design = ~ date+sex)
ddsM <- DESeq(ddsM)
resM <- results(ddsM,contrast=c("sex","F","M"))

counts(ddsM,normalized=TRUE)["Potri.019G047300.0",]
resM["Potri.019G047300.0",]
mcols(ddsM)[names(rowData(ddsM)) == "Potri.019G047300.0",]

## reclassified
sexR <- sex
sexR[5] <- "M"
ddsR <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       sex=sexR,
                       date=date),
  design = ~ date+sex)
ddsR <- DESeq(ddsR)
resR <- results(ddsR,contrast=c("sex","F","M"))

counts(ddsR,normalized=TRUE)["Potri.019G047300.0",]
resR["Potri.019G047300.0",]
mcols(ddsR)[names(rowData(ddsR)) == "Potri.019G047300.0",]

## samples details
sex
date

## the model should not be affected by the reclassification
lapply(split(sex,date),table)
lapply(split(sexR,date),table)

## sessionInfo
sessionInfo()

