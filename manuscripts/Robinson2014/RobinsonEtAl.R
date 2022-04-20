#' ---
#' title: "RNA-Seq tutorial"
#' author: "Nicolas Delhomme, Bastian Schiffthaler"
#' date: "`r Sys.Date()`"
#' output:
#'  BiocStyle::html_document: 
#'     toc: true
#'     number_sections: true
#' ---

#' # Introduction
#' As a running example, we use a dataset derived from a study performed
#' in Populus tremula (Robinson, Delhomme et al., BMC Plant Biology, 2014}. 
#' The authors test the evolutionary theory suggesting
#' that males and females may evolve sexually dimorphic phenotypic and biochemical 
#' traits concordant with each sex having different optimal strategies of resource 
#' investment to maximise reproductive success and fitness. Such sexual dimorphism
#' should result in sex biased gene expression patterns in non-floral organs for
#' autosomal genes associated with the control and development of such phenotypic 
#' traits. Hence, the authors, among other approaches have looked for gene expression 
#' differences between sex. This was achieved by an RNA-Seq differential expression 
#' analysis between samples grouped by sex. The samples have been sequenced on an
#' Illumina HiSeq 2000, using TruSeq adapters, through a 101 cycles paired-end protocol,
#' which yielded between 11 million and 34 million reads per sample. 
#' 
#' In the following sections, we look at a subset of the data,
#' corresponding to reads obtained from lanes of the RNA-seq
#' experiment, and aligned to the Populus trichocarpa reference genome. The reason why the reads
#' were aligned to the Populus trichocarpa genome is that there is no reference genome for Populus tremula
#' and Populus trichocarpa is a closely related species (the latter grows mostly in northern america 
#' while the former grows in northern and central europe). 
#' 

#' # Reproducing the DE analysis

#' ## Setup - libPaths
.libPaths(c("/home/training/R-packages",.libPaths()))

#' ## Setup - loading the libraries
suppressPackageStartupMessages(library(DESeq))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(arrayQualityMetrics))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(RnaSeqTutorial))

#' ## Setup - defining the working directory
#' Change the following value to
#' /media/sf_shared/EMBO/preprocessing/
setwd("/Users/delhomme/Desktop/tmp/Day2/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/Users/delhomme/Desktop/tmp/Day2/")
#' ```

#' ## Processing the data
#' ### Reading in the data
#' First we read the HTSeq files in a matrix. The DESeq2 package now
#' actually has a function to ease that process: DESeqDataSetFromHTSeqCount. 
#' Here we just process the samples in parallel using mclapply instead.
res <- mclapply(dir("HTSeq",pattern="^[2,3].*_STAR\\.txt",
full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=2)
names(res) <- gsub("_.*_STAR\\.txt","",dir("HTSeq",pattern="^[2,3].*_STAR\\.txt"))

#' Then we extract the additional info that HTSeq writes at the end of every
#' file detailing how many reads
addInfo <- c("__no_feature","__ambiguous",
             "__too_low_aQual","__not_aligned",
             "__alignment_not_unique")

#' and we then extract the reads
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]

#' ### the HTSeq stat lines
#' Here we aggreagte the information about how many
#' reads aligned, how many were ambiguous, etc
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' Then we convert them as percentages and check them all
apply(count.stats,2,function(co){round(co*100/sum(co))})

#' As can be seen an average 82% of the reads aligned uniquely
#' unambiguously to features. About 13% were mapping multiple 
#' locations, 2% were on ambigous features (as the data is non
#' strand specific, we have used the "Union" counting scheme) and
#' finally 3% align to no features. This is not at all unexpected
#' given the fact that the alignments were done against the 
#' related species P. trichocarpa.

#' ### Visualizing the HTSeq stats
pal=brewer.pal(6,"Dark2")[1:nrow(count.stats)]
mar <- par("mar")
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=pal,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,2e+7))
legend("top",fill=pal,legend=rownames(count.stats),bty="n",cex=0.8,horiz=TRUE)
par(mar=mar)

#' Here, one can clearly see the library size effect; the samples 229.1 and 349.2
#' have been sequenced deeper as these are the parents of a large scale F1 population
#' established for another study. Appart from these, the amount of reads per library
#' seems fairly equally distributed.

#' ### Deriving more information
#' We can estimate how much of the readsare not expressed
sel <- rowSums(count.table) == 0
sprintf("%s percent",round(sum(sel) * 100/ nrow(count.table),digits=1))
sprintf("of %s genes are not expressed",nrow(count.table))
#' So 14.2% of the genes are not expressed out of a total of 41335 genes

#' ## Biological QA
#' To assess the validity of the replicates, we can look at
#' different metrics. 

#' ###  Raw count distribution
#' We start by looking at the
#' per-gene mean expression, i.e. the mean raw count of every 
#' gene across samples is calculated and displayed on a log10 scale
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' Then the same is done for the individual samples colored by sample type
pal=brewer.pal(8,"Dark2")
plot.multidensity(log10(count.table),col=rep(pal,each=3),
                  legend.x="topright",legend.cex=0.5,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' As one can see, most samples show the similar trends, a few samples are
#' shifted to the right - those that were sequenced deeper, but this does not
#' affect the global shape of the curve.

#' ## Variance stabilisation for better visualisation
#' For visualization, the data is submitted to a variance 
#' stabilization transformation using DESeq2. The 
#'  dispersion is estimated independently of the sample type
#' We nonetheless define the metadata fo importance, i.e. every
#' sample sex and year of sampling
samples <- read.csv("~/Git/UPSCb/projects/sex/doc/sex-samples.csv")
conditions <- names(res)
sex <- samples$sex[match(names(res),samples$sample)]
date <- factor(samples$date[match(names(res),samples$sample)])

#' And then we create the DESeq2 object, using the sample name
#' as condition (which is hence irrelevant to any comparison)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       sex=sex,
                       date=date),
  design = ~ condition)

#' Once this is done, we first check the size factors and
#' as there is no big variation, we decide to go for a variance
#' stabilisation transformation over a regularised log transformation.
#' Check the DESeq2 vignette for more details about this.
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- names(res)
sizes

#' Let us do the vst
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)

#' The vst introduces an offset, i.e. non expressed gene will have a value set
#' as the minimal value of the vst. To avoid lengthy unnecessary explanations, we
#' simply shift all values so that non expressed genes have a vst value of 0
vst <- vst - min(vst)

#' It is then essential to validate the vst effect and ensure its validity. To
#' that extend, we visualize the corrected mean - sd relationship. As it is 
#' fairly linear, meaning we can assume homoscedasticity.
#' the slight initial trend / bump is due to genes having few 
#' counts in a few subset of the samples and hence having a 
#' higher variability. This is expected.
meanSdPlot(assay(vsd)[rowSums(count.table)>0,], ylim = c(0,2.5))

#' ### PCA analysis
#' We perform a Principal Component Analysis (PCA) of the data
#' to do a quick quality assessment, i.e. replicate should 
#' cluster and the first dimensions should be explainable by 
#' biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
smpls <- conditions

#' We color the samples by date and sex
sex.cols<-c("pink","lightblue")
sex.names<-c(F="Female",M="Male")

#' And then check wheter the observed separation is due to the sample sex
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

#' This does not seem to be the case at all. There are 2 clusters of points, but
#' both equally contain pink and blue dots. So, we next color by sampling date
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

#' And here we are... The sampling data is a STRONG CONFOUNDING FACTOR in our
#' analysis. So let's do a final plot with the color = sex and shape = date
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

#' Sometimes the 3D PCA are not so easy to read and we may be just
#' as well off looking at 2 dimensions at a time. First the 2 first dims
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

#' Clearly the sampling date separate most samples. It is interesting to see that
#' the 2 samples sequenced deeper are in between the 2 clusters. We could hypothesize
#' that the environmental change between sampling years are affecting less important 
#' genes, hence genes which expression levels are lower; i.e. somewhat looking at 
#' transcriptional "noise". Since the 2 parent samples have been sequenced much deeper, 
#' the fact that they group with neither cluster might be due to that fact that they too have 
#' an increased proportion of transcriptional noise.

#' Looking at the 2nd and 3rd dims does not reveal any obvious effects
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

#' ## Doing the differential expression analysis.
#' The QA revealed that we have a confounding factor. Luckily, enough
#' the grouping of sex is balanced between sampling years and hence 
#' we will be able to block the year effect to investigate the sex effect.

yearSexDesign <- samples[grep("[2,3].*",samples$sample),c("sex","date")]

#' ### DESeq
#' The primary reason to use DESeq over DESeq2 or edgeR or 
#' any other is that DESeq is very conservative. The second
#' reason is that DESeq give less weight to outliers than
#' DESeq2 does and based on the study by Pakull et al., which
#' identify a candidate gene for sexual determination - the
#' same we simultaneously, independantly originally identified
#' in our analyses - it seems that our sample 226.1 may be
#' mis-classified. However, the sex phenotyping has been done
#' thoroughly, so it may also be that the sex determination
#' is a more complex trait and that Potri.019G047300 is only
#' one factor influencing it

#' As introduced we want to look at the sex by blocking the year effect.
#' We start by estimating the size factors and the dispersion
cdsFull <- newCountDataSet(count.table,yearSexDesign)
cdsFull = estimateSizeFactors( cdsFull )
cdsFull = estimateDispersions( cdsFull )

#' We can now check the dispersion estimated by DESeq, which is, obviously, conservative.
#' I.e. most dispersion values lie below the actual dispersion fit by about 1 fold.
plotDispLSD(cdsFull)

#' Next we create both models (one considering the date only and on the date and sex)
fit1 = suppressWarnings(fitNbinomGLMs( cdsFull, count ~ date + sex ))
fit0 = suppressWarnings(fitNbinomGLMs( cdsFull, count ~ date))

#' For the rest of the analysis, we ignore the genes that did not converge in the 
#' previous step
sel <- !is.na(fit1$converged) & !is.na(fit0$converged) & fit1$converged & fit0$converged
fit1 <- fit1[sel,]
fit0 <- fit0[sel,]

#' We next calculate the GLM p-value and adjust them for multiple testing
pvalsGLM = nbinomGLMTest( fit1, fit0 )
padjGLM = p.adjust( pvalsGLM, method="BH" )

#' We finally visualize the obtained results. A first insight shows that a number
#' of genes have very high fold changes, which is due the fact that these genes
#' have very few read in very few samples. For clarity, in the following plots,
#' we filter those genes with a too high FC.
boxplot(rowSums(count.table[rownames(fit1[fit1$sexM > 20,]),]>0),
        ylab="number of sample with expression > 0",
        main="distribution of the # of samples with\nexpression > 0 for genes with extreme high log2FC"
)
boxplot(rowSums(count.table[rownames(fit1[fit1$sexM < -20,]),]>0),
        ylab="number of sample with expression > 0",
        main="distribution of the # of samples with\nexpression > 0 for genes with extreme low log2FC"
)

#' The vast majority of these genes are anyway not significant.
plot(density(padjGLM[fit1$sexM > 20]),
     main="Adj. p-value for genes with extreme log2FC",
     xlab="Adj. p-value")
plot(density(padjGLM[fit1$sexM > -20]),
     main="Adj. p-value for genes with extreme log2FC",
     xlab="Adj. p-value")

#' So we filter them out
sel <- fit1$sexM > -20 & fit1$sexM < 20
fit1 <- fit1[sel,]
padjGLM <- padjGLM[sel]
pvalsGLM <- pvalsGLM[sel]

#' A further look into the data reveals that two p-values are equal to 0. As
#' this is inadequate for plotting log odds (-log of the p-value), we change
#' these 0s to be 1 log10 fold change smaller than the otherwise smallest p-value
#' (for plotting purposes only!)
padjGLM[padjGLM==0] <- min(padjGLM[padjGLM!=0])/10
pvalsGLM[pvalsGLM==0] <- min(pvalsGLM[pvalsGLM!=0])/10

#' ### Visualize DESeq results
#' volcano plot with the non-adjusted p-values
heatscatter(fit1$sexM,-log10(pvalsGLM),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) p-value")
legend("topleft",bty="n","1% p-value cutoff",lty=2,col="gray")
points(fit1[pvalsGLM<.01,"sexM"],-log10(pvalsGLM[pvalsGLM<.01]),col="lightblue",pch=19)
pos <- c(order(pvalsGLM)[1:4],order(fit1$sexM,decreasing=TRUE)[1:4],order(fit1$sexM)[1:4])
points(fit1[pos,"sexM"],-log10(pvalsGLM[pos]),col="red")
abline(h=2,lty=2,col="gray")

#' volcano plot with adjusted p-values
heatscatter(fit1$sexM,-log10(padjGLM),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) adj. p-value")
legend("topleft",bty="n","1% FDR cutoff",lty=2,col="gray")
points(fit1[padjGLM<.01,"sexM"],-log10(padjGLM[padjGLM<.01]),col="lightblue",pch=19,cex=1.5)
pos <- c(order(padjGLM)[1:4],order(fit1$sexM,decreasing=TRUE)[1:4],order(fit1$sexM)[1:4])
points(fit1[pos,"sexM"],-log10(padjGLM[pos]),col="red",cex=c(2,2,rep(1,length(pos)-2)))
abline(h=2,lty=2,col="gray")

#' As DESeq is an "older" implementation of this type of analysis, we will 
#' redo the analysis using DESeq2. To ultimately be able to compare the
#' obtained set of differential expression candidate genes, we save the list
#' of genes significant at a 1% adjusted p-value cutoff.
UmAspD <- rownames(fit1[padjGLM<.01,])

#' which are not so many
UmAspD

#' ### DESeq2 differential expression analyses
#' We redo the analyses performed above. As you will see, it has been made
#' more intuitive in DESeq2 and the same can be achieved in fewer commands.
design(dds) <- ~date+sex
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds,contrast=c("sex","M","F"))

#' In 4 commands we have reproduced the analysis and we have observed that the
#' dispersion estimation looked different. DESeq2 introduces a so called 
#' "shrinkage", which you can learn more about in the DESeq2 vignettes. 
#' We, next, extract the results using the same 1% cutoff and plot similar
#' validation plots
alpha=0.01
plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
#' The volcano plot look similar to what we observed previously (except that we
#' have by mistake inverted the ratio calculation, we now did F-M), although
#' slightly less dead that in DESeq. We also do not get any genes with a p-value
#' of 0
hist(res$padj,breaks=seq(0,1,.01))
sum(res$padj<alpha,na.rm=TRUE)
#' 4 genes are candidates for differential expression at a 1% adjusted p-value
#' cutoff, which we also save for later comparison
UmAspD2 <- rownames(res[order(res$padj,na.last=TRUE),][1:sum(res$padj<alpha,na.rm=TRUE),])
UmAspD2

#' We can though already observe that one gene is common with the previous list,
#' but that the chromosome 19 candidate has disappeared. This is surprising, given
#' the Pakull et al. publication, that Potri.019G047300.0 does not come out anymore as
#' a differentially expressed candidate gene at a 1% FDR cutoff with
#' DESeq2. But it's adjusted p-value is 2.1%. Moreover it's
#' cook distance - a measure of the homogeneity of the gene
#' expression within a condition - is relatively high (2x
#' the average) and might indicate the presence of an outlier. 
#' Looking into more details at the count values and sex 
#' assignment shows that the sample 226.1 behaves like a male 
#' sample with regards to that gene, so either the sample has been mis-sexed
#' or the Pakull et al. results in P. tremuloides and P.tremula are only a partial
#' view of the true biology. 
res["Potri.019G047300.0",]
data.frame(
  counts=counts(dds,normalized=TRUE)["Potri.019G047300.0",],
  sex,
  date,
  conditions)

#' ### DESeq2 with the sample re-sexed
#' Nevertheless, while waiting for the next time the tree
#' will flower and we can confirm its sex with certainty, we can assume a that the 
#' sample was mis-sexed and see what effect a sex correction would have.
#' Given the data from Pakull et al. about that gene; 
#' assuming that the sample 226.1 was mis-sexed, we redo the 
#' DESeq2 analysis after swaping the sex of that sample
sex[conditions=="226.1"] <- "M"
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       sex=sex,
                       date=date),
  design = ~ date+sex)
dds <- DESeq(dds)
res.swap <- results(dds,contrast=c("sex","M","F"))
plotMA(res.swap,alpha=alpha)
volcanoPlot(res.swap,alpha=alpha)
hist(res.swap$padj,breaks=seq(0,1,.01))
sum(res.swap$padj<alpha,na.rm=TRUE)
UmAspD2.swap <- rownames(res.swap[order(res.swap$padj,na.last=TRUE),][1:sum(res.swap$padj<alpha,na.rm=TRUE),])
UmAspD2.swap
#' fair enough, the gene made it back in the list: Potri.019G047300.0
#' But, still, to assess the importance of the 226.1 sample,
# we redo the DESeq2 DE analysis without it (since we have enough replicate for
#' that removal not to affect the power of our analysis)
sel <- conditions != "226.1"
dds <- DESeqDataSetFromMatrix(
  countData = count.table[,sel],
  colData = data.frame(condition=conditions[sel],
                       sex=sex[sel],
                       date=date[sel]),
  design = ~ date+sex)
dds <- DESeq(dds)
plotDispEsts(dds)
res.sel <- results(dds,contrast=c("sex","M","F"))
plotMA(res.sel,alpha=alpha)
volcanoPlot(res.sel,alpha=alpha)
hist(res.sel$padj,breaks=seq(0,1,.01))
sum(res.sel$padj<alpha,na.rm=TRUE)
UmAspD2.sel <- rownames(res.sel[order(res.sel$padj,na.last=TRUE),][1:sum(res.sel$padj<alpha,na.rm=TRUE),])
UmAspD2.sel

#' ## Results comparison
#' We compare the results from the 4 approaches, the easiest way being to plot
#' a Vann Diagram
plot.new()
grid.draw(venn.diagram(list(
  UmAsp=UmAspD2,
  UmAsp.swap=UmAspD2.swap,
  UmAsp.sel=UmAspD2.sel,
  UmAspD = UmAspD),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("UmAsp - DESeq2",
                                        "UmAsp - corrected",
                                        "UmAsp - removed",
                                        "UmAsp - DESeq")))

#' Potri.014G155300.0 is constantly found by all 4 approaches
#' Potri.019G047300.0 is found by DESeq and DESeq2, when
#' the putative mis-sexed sample is removed or re-classified
#' the sample sex correction actually leads to the identification
#' of 16 genes unique to that comparison! It is still unclear if
#' these are the results of the technical error or if they would be
#' worth investigating further.
sort(intersect(UmAspD2,UmAspD))
sort(intersect(UmAspD2.swap,UmAspD))
sort(intersect(UmAspD2.sel,UmAspD))

#' This conclude this differential expression analysis. The following details
#' the R version that was used to perform the analyses.

#' # sessionInfo
sessionInfo()

#' # Appendix
#' Analyses performed when asking Mike Love (DESeq2 developer)
#' why DESeq2 seem so sensitive to misclassification
#' The short answer was to use the Cook distance to further
#' investigate that and that the 1% FDR cutoff was 
#' conservative. These are commented out on purpose.

## mis-classified
# sex <- samples$sex[match(colnames(count.table),samples$sample)]
# date <- factor(samples$date[match(colnames(count.table),samples$sample)])
# ddsM <- DESeqDataSetFromMatrix(
#   countData = count.table,
#   colData = data.frame(condition=conditions,
#                        sex=sex,
#                        date=date),
#   design = ~ date+sex)
# ddsM <- DESeq(ddsM)
# resM <- results(ddsM,contrast=c("sex","F","M"))
# 
# counts(ddsM,normalized=TRUE)["Potri.019G047300.0",]
# resM["Potri.019G047300.0",]
# mcols(ddsM)[names(rowData(ddsM)) == "Potri.019G047300.0",]
# 
# ## reclassified
# sexR <- sex
# sexR[5] <- "M"
# ddsR <- DESeqDataSetFromMatrix(
#   countData = count.table,
#   colData = data.frame(condition=conditions,
#                        sex=sexR,
#                        date=date),
#   design = ~ date+sex)
# ddsR <- DESeq(ddsR)
# resR <- results(ddsR,contrast=c("sex","F","M"))
# 
# counts(ddsR,normalized=TRUE)["Potri.019G047300.0",]
# resR["Potri.019G047300.0",]
# mcols(ddsR)[names(rowData(ddsR)) == "Potri.019G047300.0",]
# 
# ## samples details
# sex
# date
# 
# ## the model is not be affected by the reclassification
# lapply(split(sex,date),table)
# lapply(split(sexR,date),table)
# 