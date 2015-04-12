#' ---
#' title: "RNA-Seq tutorial"
#' author: "Nicolas Delhomme, Bastian Schiffthaler"
#' date: "`r Sys.Date()`"
#' output:
#'  BiocStyle::html_document:
#'     toc: true
#'     number_sections: true
#' bibliography: ~/Git/UPSCb-public/tutorial/Robinson-Delhomme-et-al_BMC-Plant-Biology_2014.bib
#' ---

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="~/Git/UPSCb-public")
#' ```

#' # Introduction
#' This tutorial introduces and RNA-Seq differential expression analysis performed
#' using R and Bioconductor[@ref:R, @Gentleman:2004p2013].
#' As a running example, we use a dataset derived from a study performed
#' in _Populus tremula_ [@Robinson:2014p6362].
#' The study looks at the evolutionary theory suggesting
#' that males and females may evolve sexually dimorphic phenotypic and biochemical
#' traits concordant with each sex having different optimal strategies of resource
#' investment to maximise reproductive success and fitness. Such sexual dimorphism
#' should result in sex biased gene expression patterns in non-floral organs for
#' autosomal genes associated with the control and development of such phenotypic
#' traits. Hence, the study among other approaches used gene expression to look for
#' differences between sex. This was achieved by an RNA-Seq differential expression
#' analysis between samples grouped by sex. The samples have been sequenced on an
#' Illumina HiSeq 2000, using TruSeq adapters, through a 101 cycles paired-end protocol,
#' which yielded between 11 million and 34 million reads per sample.
#'
#' In the following sections, we use the RNA-Seq count data that has been obtained
#' by aligning the pre-processed (adaptor and quality trimmed, rRNA filtered) reads
#' to the _Populus trichocarpa_ reference genome. The reason why the reads
#' were aligned to the _Populus trichocarpa_ genome is that there is no reference
#' genome for _Populus tremula_ and _Populus trichocarpa_ is a closely related
#' species (the latter grows in north-western America while the former grows
#' in northern and central Europe).
#'

#' # Reproducing the DE analysis

#' ## Setup - loading the libraries
suppressPackageStartupMessages(library(DESeq))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(easyRNASeq))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(LSD))
source("src/R/plot.multidensity.R")
source("src/R/volcanoPlot.R")

#' ## Processing the data
#' ### Reading in the data
#' First we read the count files produced by HTSeq[@Anders:2014p6365] in a matrix.
#' The DESeq2[@Love:2014p6358] package now actually has a function to ease that process: **DESeqDataSetFromHTSeqCount**.
#' Here we just process the samples in parallel using mclapply instead.
res <- mclapply(dir("data/htseq-count",pattern="^[2,3].*_STAR\\.txt",
full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=2)
names(res) <- gsub("_.*_STAR\\.txt","",dir("data/htseq-count",pattern="^[2,3].*_STAR\\.txt"))

#' Then we extract the additional information that HTSeq writes at the end of every
#' file detailing the number of reads that were not taken into account while
#' counting summarized by the reason they were rejected.
addInfo <- c("__no_feature","__ambiguous",
             "__too_low_aQual","__not_aligned",
             "__alignment_not_unique")

#' Finally, we extract the reads
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]

#' ### The HTSeq stat lines
#' Here we aggregate the information about how many
#' reads aligned together with the information gathered above.
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' And convert them as percentages to check them all
apply(count.stats,2,function(co){round(co*100/sum(co))})

#' As can be seen, an average 82% of the reads aligned uniquely
#' unambiguously to features. About 13% were mapping multiple
#' locations, 2% were on ambigous features (as the data is non
#' strand specific, we have used the "Union" counting scheme) and
#' finally 3% align to no features. This is not at all unexpected
#' given the fact that the alignments were done against the
#' related species P. trichocarpa.
#'```{r empty, eval=FALSE,echo=FALSE}
#'```
#' ### Visualizing the HTSeq stats
pal=brewer.pal(6,"Dark2")[1:nrow(count.stats)]
mar <- par("mar")
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=pal,beside=TRUE,
        las=2,main="read proportion",
        ylim=range(count.stats) + c(0,2e+7))
legend("top",fill=pal,legend=rownames(count.stats),
       bty="n",cex=0.8,horiz=TRUE)
par(mar=mar)

#' Here, one can clearly see the library size effect; the samples 229.1 and 349.2
#' have been sequenced deeper as these are the parents of a large scale F1 population
#' established for another study. Appart from these, the amount of reads per library
#' seems fairly equally distributed.
#'```{r empty, eval=FALSE,echo=FALSE}
#'```
#' ### Deriving more information
#' It is also important to estimate how much of the reads are not expressed
sel <- rowSums(count.table) == 0
sprintf("%s percent of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count.table),digits=1)
        ,nrow(count.table))
#' So 14.2% of the genes are not expressed out of a total of 41,335 genes
#'```{r empty, eval=FALSE,echo=FALSE}
#'```
#' ## Biological QA
#' To assess the validity of the replicates, we can look at
#' different metrics.
#'```{r empty, eval=FALSE,echo=FALSE}
#'```
#' ###  Raw count distribution
#' We start by looking at the
#' per-gene mean expression, _i.e._ the mean raw count of every
#' gene across samples is calculated and displayed on a log10 scale
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' This shows an average of ~300X coverage. Next, the same is done for the
#' individual samples colored by sample type
pal=brewer.pal(8,"Dark2")
plot.multidensity(log10(count.table),col=rep(pal,each=3),
                  legend.x="topright",legend.cex=0.5,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' As one can see, the samples show a similar trend, although few samples are
#' shifted to the right - those that were sequenced deeper, but this does not
#' affect the global shape of the curve.
#'```{r empty, eval=FALSE,echo=FALSE}
#'```
#' ## Variance stabilisation for better visualisation
#' For visualization, the data is submitted to a variance
#' stabilization transformation using DESeq2[@Love:2014p6358]. The
#'  dispersion is estimated independently of the sample type
#' We nonetheless define the metadata of importance, _i.e._ every
#' sample sex and year of sampling
samples <- read.csv("~/Git/UPSCb/projects/sex/doc/sex-samples.csv")
conditions <- names(res)
sex <- samples$sex[match(names(res),samples$sample)]
date <- factor(samples$date[match(names(res),samples$sample)])

#' Finally, we create the DESeq2 object, using the sample name
#' as condition (which is hence irrelevant to any comparison)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       sex=sex,
                       date=date),
  design = ~ condition)

#' Once this is done, we first check the sequencing library sizes (a.k.a. the
#' size factors) variation.
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- names(res)
sizes

#' As there is no big variation, we decide to go for a variance
#' stabilisation transformation over a regularised log transformation.
#' Check the DESeq2 [vignette](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf) for more details about this.
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)

#' The vst introduces an offset, _i.e._ non expressed gene will have a value set
#' as the minimal value of the vst. To avoid confusion, we
#' simply shift all values so that non expressed genes have a vst value of 0.
#' This is **acceptable** because we can now assume a linear relationship in the
#' expression values, but really because we only **do this for visualisation** and
#' **not for statistical analyses**.
vst <- vst - min(vst)

#' It is then essential to validate the vst effect and ensure its validity. To
#' that extend, we visualize the corrected mean - sd relationship.
meanSdPlot(assay(vsd)[rowSums(count.table)>0,], ylim = c(0,2.5))

#' As it is
#' fairly linear, we can assume homoscedasticity.
#' the slight initial trend / bump is due to genes having few
#' counts in a few subset of the samples and hence having a
#' higher variability. This is expected.
#'```{r empty, eval=FALSE,echo=FALSE}
#'```
#' ### PCA analysis
#' We perform a Principal Component Analysis (PCA) of the data
#' to do a quick quality assessment, _i.e._ replicate samples should
#' cluster together and the first 2-3 dimensions should be explainable by
#' biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
smpls <- conditions

#' At first, we color the samples by date and sex
sex.cols<-c("pink","lightblue")
sex.names<-c(F="Female",M="Male")

#' And check whether the observed separation is due to the sample sex
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

#' This is obviously not the case at all. There are 2 clusters of points, but
#' both contain pink and blue dots. So, we next color by sampling date
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

#' And here we **realise that the sampling date is a STRONG CONFOUNDING FACTOR** in our
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

#' Sometimes 3D PCA are not easy to interpret, so we may just
#' as well look at 2 dimensions at a time.
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

#' Clearly the sampling date separate most samples, _i.e._ first dimension equals "time".
#' It is interesting to see that
#' the 2 samples sequenced deeper are in between the 2 clusters. We could hypothesize
#' that the environmental change between sampling years are affecting a subset of
#' genes, which expression levels are presumably lower; _i.e._ somewhat looking at
#' transcriptional "noise" environmentaly-induced (difference in temperature, moisture,_etc._).
#' Since the 2 parent samples have been sequenced much deeper, the fact that they do not group
#' with either cluster might be due to the fact that they, too, have
#' an increased proportion of transcriptional "noise" possibly technically-induced; _i.e._
#' sequencing so deep that we start to measure non-mature mRNA, mis-splice mRNA, _etc._ - after
#' all transcription has to be a stochastic process.

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
#' The primary reason to use DESeq[Anders:2010p1659] over DESeq2 or edgeR[Robinson:2010p775] or
#' any other is that DESeq is very conservative (see Soneson & Delorenzi [-@Soneson:2013p5778]). The second
#' reason is that DESeq give less weight to outliers than
#' DESeq2 does. Based on the study by Pakull _et al._ [-@Pakull], which
#' identified a candidate gene: **Potri.019G047300** for sexual determination - the
#' same we simultaneously, independantly originally identified
#' in our analyses - it seems that our sample 226.1 may be
#' mis-classified, so our first approach uses DESeq.
#'
#' However, the sex phenotyping had been done
#' thoroughly, so it may also be that the sex determination
#' is a more complex trait and that Potri.019G047300 is only
#' one factor influencing it
#'
#' As introduced we want to look at the sex effect blocking for the year effect.
#' We start by estimating the size factors and the dispersion
cdsFull <- newCountDataSet(count.table,yearSexDesign)
cdsFull = estimateSizeFactors( cdsFull )
cdsFull = estimateDispersions( cdsFull )

#' We can now check the dispersion estimated by DESeq, which is, obviously, conservative.
#' I.e. most dispersion values lie below the actual dispersion fit by about 1 fold.
plotDispLSD(cdsFull)

#' Next, we create both models (one considering the date only and one the date and sex combination)
fit1 = suppressWarnings(fitNbinomGLMs( cdsFull, count ~ date + sex ))
fit0 = suppressWarnings(fitNbinomGLMs( cdsFull, count ~ date))

#' For the rest of the analysis, we ignore the genes that did not converge in the
#' previous step
sel <- !is.na(fit1$converged) & !is.na(fit0$converged) & fit1$converged & fit0$converged
fit1 <- fit1[sel,]
fit0 <- fit0[sel,]

#' We calculate the GLM p-value and adjust them for multiple testing
pvalsGLM = nbinomGLMTest( fit1, fit0 )
padjGLM = p.adjust( pvalsGLM, method="BH" )

#' Finally, we visualize the obtained results. A first insight shows that a number
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
range(padjGLM[fit1$sexM > 20])
plot(density(padjGLM[fit1$sexM > 20]),
     main="Adj. p-value for genes with extreme log2FC",
     xlab="Adj. p-value")
range(padjGLM[fit1$sexM < -20])
plot(density(padjGLM[fit1$sexM < -20]),
     main="Adj. p-value for genes with extreme log2FC",
     xlab="Adj. p-value")

#' So we filter them out
sel <- fit1$sexM > -20 & fit1$sexM < 20
plot(density(padjGLM[sel]),
     main="Adj. p-value for genes without extreme log2FC",
     xlab="Adj. p-value")
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
#' An intuitively easy to intrepret visualisation is the (volcanoplot)[http://en.wikipedia.org/wiki/Volcano_plot_(statistics)].
#' We plot both the non-adjusted and adjusted p-values to also
#' visualise the multiple correction effect.
#' First, the volcano plot with the non-adjusted p-values
heatscatter(fit1$sexM,-log10(pvalsGLM),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) p-value")
legend("topleft",bty="n","1% p-value cutoff",lty=2,col="gray")
points(fit1[pvalsGLM<.01,"sexM"],-log10(pvalsGLM[pvalsGLM<.01]),col="lightblue",pch=19)
pos <- c(order(pvalsGLM)[1:4],order(fit1$sexM,decreasing=TRUE)[1:4],order(fit1$sexM)[1:4])
points(fit1[pos,"sexM"],-log10(pvalsGLM[pos]),col="red")
abline(h=2,lty=2,col="gray")

#' and then the volcano plot with the adjusted p-values
heatscatter(fit1$sexM,-log10(padjGLM),
            main="Male vs. Female Differential Expression",cor=FALSE,
            xlab="Log2 Fold Change", ylab="- log(10) adj. p-value")
legend("topleft",bty="n","1% FDR cutoff",lty=2,col="gray")
points(fit1[padjGLM<.01,"sexM"],-log10(padjGLM[padjGLM<.01]),col="lightblue",pch=19,cex=1.5)
pos <- c(order(padjGLM)[1:4],order(fit1$sexM,decreasing=TRUE)[1:4],order(fit1$sexM)[1:4])
points(fit1[pos,"sexM"],-log10(padjGLM[pos]),col="red",cex=c(2,2,rep(1,length(pos)-2)))
abline(h=2,lty=2,col="gray")

#' As DESeq is an "older" implementation for this type of analysis, we
#' redo the analysis using DESeq2. To ultimately be able to compare the
#' obtained set of differential expression candidate genes, we save the list
#' of genes significant at a 1% adjusted p-value cutoff.
UmAspD <- rownames(fit1[padjGLM<.01,])

#' These are not so many
UmAspD

#' ### DESeq2 differential expression analyses
#' As you see, doing that type of analysis has been made
#' more intuitive in DESeq2 (_i.e._ the same can be achieved in fewer commands).
design(dds) <- ~date+sex
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds,contrast=c("sex","F","M"))

#' In solely 4 commands we have reproduced the analysis and we have observed that the
#' dispersion estimation looked different. DESeq2 introduces a so called
#' "shrinkage", which you can read more about in the DESeq2 [vignette](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf).
#' Next, we extract the results using the same 1% cutoff and plot similar
#' validation plots
alpha=0.01
plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
#' The volcano plot look similar to what we observed previously, although
#' slightly less dead that in DESeq. We also do not get any genes with a p-value
#' of 0
hist(res$padj,breaks=seq(0,1,.01))
sum(res$padj<alpha,na.rm=TRUE)
#' 4 genes are candidates for differential expression at a 1% adjusted p-value
#' cutoff, which we also save for later comparison
UmAspD2 <- rownames(res[order(res$padj,na.last=TRUE),][1:sum(res$padj<alpha,na.rm=TRUE),])
UmAspD2

#' We can observe that one gene is common with the previous list,
#' but that the chromosome 19 candidate has disappeared. This is surprising, given
#' the Pakull _et al._ publication, that Potri.019G047300.0 does not come out anymore as
#' a differentially expressed candidate gene at a 1% FDR cutoff with
#' DESeq2. Its adjusted p-value is actually 2.1%. Moreover it's
#' cook distance - a measure of the homogeneity of the gene
#' expression within a condition - is relatively high (2x
#' the average) and might indicate the presence of an outlier.
#' Looking closely at the count values and sex
#' assignment shows that the sample 226.1 behaves like a male
#' sample with regards to that gene, so either the sample has been mis-sexed
#' or the Pakull _et al._ results in P. tremuloides and P.tremula are only a partial
#' view of the true biology.
res["Potri.019G047300.0",]
data.frame(
  counts=counts(dds,normalized=TRUE)["Potri.019G047300.0",],
  sex,
  date,
  conditions)
#' Nevertheless, while waiting for the next time the tree
#' will flower and we can confirm its sex with certainty, we can assume that the
#' sample was mis-sexed and see what effect a sex correction would have.
#'```{r empty, eval=FALSE,echo=FALSE}
#'```
#' ### DESeq2 with the sample re-sexed
#' Given the data from Pakull _et al._ about that gene;
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
#' Fair enough, the gene Potri.019G047300.0 made it back in the list.
#' But, still, to assess the importance of the 226.1 sample,
#' we redo the DESeq2 DE analysis without it (since we have enough replicate for
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
#' a [Venn Diagram](http://en.wikipedia.org/wiki/Venn_diagram)
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
#' the putative mis-sexed sample is removed or re-classified.
#'
#' The sample sex correction actually leads to the identification
#' of 16 genes unique to that comparison! It is unclear if
#' these are the results of the technical error or if they would be
#' worth investigating further.
sort(intersect(UmAspD2,UmAspD))
sort(intersect(UmAspD2.swap,UmAspD))
sort(intersect(UmAspD2.sel,UmAspD))

#' This conclude this differential expression analysis. The following details
#' the R version that was used to perform the analyses.
#'```{r empty, eval=FALSE,echo=FALSE}
#'```
#' # Session Info
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
#' # References
