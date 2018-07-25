rm(list=ls(all=TRUE))
#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme & Karen Kloth"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' # Environment
#' Set the working dir
setwd("/mnt/picea/projects/arabidopsis/balbrectsen/arabidopsis-aphids")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/arabidopsis/balbrectsen/arabidopsis-aphids")
#' ```

#' Libs
#'```{r drop tcltk, echo=FALSE}
#' options(gsubfn.engine = "R")
#'```
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(plotrix))
suppressPackageStartupMessages(library(sciplot))

#' Helper files
suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))

#' Load DESeq object
load("analysis/kallisto/DESeq-object_combtechrep_batch.rda")

#' Quickly VST it
vst <- varianceStabilizingTransformation(dds,blind=FALSE)

#' Extract the normalised counts
vsd <- assay(vst)

#' Which we remove, so that 0 means no expression
vsd <- vsd - min(vsd)
colnames(vsd)

#' Setup graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

# ----
#' # Process
#' Differential Expression
dds <- DESeq(dds)

#' Join the two factors
dds$Group <- factor(paste0(dds$Genotype, dds$Treatment))

#' # CONTRAST: compare genotypes within control treatment, Col-0 c as reference 
dds$Group <- relevel(dds$Group,"Col-0c")
design(dds) <- ~Batch + Group
dds <- DESeq(dds)
resGenoc <- results(dds,name="Group_pae9.2c_vs_Col.0c")

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selGenoc <- resGenoc$padj <= padj &
  abs(resGenoc$log2FoldChange) >= lfc &
  ! is.na(resGenoc$padj)

#' Get Gene names
gnames <- read.csv("analysis/kallisto/DEgenes_Genoc_genenames.csv")
colnames(gnames)
gnames <- gnames[order(gnames$Genes), c(1,3,4)]

#' Heatmap with gene clustering based on pae9-2 
#' mean Z score per genotype*treatment
colnames(vsd)
selGenoc_colc <- vsd[selGenoc,grep("Col-0_c", colnames(vsd))]
selGenoc_col4 <- vsd[selGenoc,grep("Col-0_4hpi", colnames(vsd))]
selGenoc_col8 <- vsd[selGenoc,grep("Col-0_8hpi", colnames(vsd))]
selGenoc_paec <- vsd[selGenoc,grep("pae9-2_c", colnames(vsd))]
selGenoc_pae4 <- vsd[selGenoc,grep("pae9-2_4hpi", colnames(vsd))]
selGenoc_pae8 <- vsd[selGenoc,grep("pae9-2_8hpi", colnames(vsd))]
selGenoc_mean <- cbind(apply(selGenoc_paec, 1, mean), apply(selGenoc_pae4, 1, mean), apply(selGenoc_pae8, 1, mean),
                       apply(selGenoc_colc, 1, mean), apply(selGenoc_col4, 1, mean), apply(selGenoc_col8, 1, mean))
colnames(selGenoc_mean) <- c("pae9-2_c", "pae9-2_4hpi", "pae9-2_8hpi", "Col-0_c", "Col-0_4hpi", "Col-0_8hpi")
selGenoc_meanscale <- t(scale(t(selGenoc_mean)))

#' row dendrogram pertubations based on pae9-2
clusterr     <-  hclust(pearson.dist(selGenoc_meanscale[,1:3]),method = "ward.D") # col 1:3 is pae9-2
dendrogramr  <-  as.dendrogram(clusterr)

#' Plot heatmap
pdf("analysis/kallisto/Heatmap_DEG_control.pdf", width = 7, height = 8)
par(oma = c(4, 0, 0, 8))
heatmap.2(selGenoc_mean, 
          trace="none", cexRow = 0.6,
          col=hpal, scale="row",
          hclustfun = function(X){hclust(X,method="ward.D")},
          labCol = colnames(selGenoc_mean),
          dendrogram = "row", Colv = FALSE, Rowv = dendrogramr,
          labRow = paste(rownames(selGenoc_mean), gnames$Nameshort),
          rowsep=c(20, 40), sepwidth = c(0.3, 0.3))
dev.off()

#' # Gene clusters and timeplots
#' z scale for gene counts (difference between mean expression and sample expression)
selGenoc_meanscale <- t(scale(t(selGenoc_mean)))
hc <- hclust(pearson.dist(selGenoc_meanscale[,1:3]),method = "ward.D")

#' from dendrogram heatmap, roughly 3 groups
tc <- cutree(hc,k=3)
nams <- split(names(tc),tc) 
barplot(elementNROWS(nams))
sapply(1:length(nams),function(i){
  clusterplot(selGenoc_meanscale[nams[[i]],4:6],
              xlabels = colnames(selGenoc_meanscale)[4:6],
              colpal = "blues",las=2)
  write(nams[[i]],file=paste0("heat-cluster",i,".txt"))
})

#' Gene names
clus <- data.frame(rbind(cbind(nams[[1]],rep("Cluster III",length(nams[[1]]))),
                         cbind(nams[[2]],rep("Cluster II",length(nams[[2]]))),
                         cbind(nams[[3]],rep("Cluster I",length(nams[[3]])))))
colnames(clus) <- c("Genes", "Cluster")
clus <- merge(clus, gnames, by = "Genes", all = TRUE)
write.csv(clus, "analysis/kallisto/Clusters_Genoc.csv")

#' Order of columns (pae9-2 c-4-8, then col-0 c-4-8)
colnames(selGenoc_meanscale)

#' Prepare timeplots of some interesting genes:
AT1G21120<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT1G21120", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT1G21120", grep("Col-0",colnames(vsd))])))
colnames(AT1G21120) <- c("Genotype", "Treatment", "normexpr")
AT1G21120$normexpr <- as.numeric(as.character(AT1G21120$normexpr))
AT1G21120$Treatment <- relevel(AT1G21120$Treatment, "control")

AT3G14440 <- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT3G14440", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT3G14440", grep("Col-0",colnames(vsd))])))
colnames(AT3G14440) <- c("Genotype", "Treatment", "normexpr")
AT3G14440$normexpr <- as.numeric(as.character(AT3G14440$normexpr))
AT3G14440$Treatment <- relevel(AT3G14440$Treatment, "control")

AT3G26830<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT3G26830", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT3G26830", grep("Col-0",colnames(vsd))])))
colnames(AT3G26830) <- c("Genotype", "Treatment", "normexpr")
AT3G26830$normexpr <- as.numeric(as.character(AT3G26830$normexpr))
AT3G26830$Treatment <- relevel(AT3G26830$Treatment, "control")

AT1G51780<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("c",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT1G51780", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT1G51780", grep("Col-0",colnames(vsd))])))
colnames(AT1G51780) <- c("Genotype", "Treatment", "normexpr")
AT1G51780$normexpr <- as.numeric(as.character(AT1G51780$normexpr))
AT1G51780$Treatment <- relevel(AT1G51780$Treatment, "control")

AT5G04950<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT5G04950", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT5G04950", grep("Col-0",colnames(vsd))])))
colnames(AT5G04950) <- c("Genotype", "Treatment", "normexpr")
AT5G04950$normexpr <- as.numeric(as.character(AT5G04950$normexpr))
AT5G04950$Treatment <- relevel(AT5G04950$Treatment, "control")

AT1G75750<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("c",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT1G75750", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT1G75750", grep("Col-0",colnames(vsd))])))
colnames(AT1G75750) <- c("Genotype", "Treatment", "normexpr")
AT1G75750$normexpr <- as.numeric(as.character(AT1G75750$normexpr))
AT1G75750$Treatment <- relevel(AT1G75750$Treatment, "control")

AT1G53480 <- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4 hpi",3), rep("8 hpi",3)),
                   c(vsd[rownames = "AT1G53480", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT1G53480", grep("Col-0",colnames(vsd))])))
colnames(AT1G53480) <- c("Genotype", "Treatment", "normexpr")
AT1G53480$normexpr <- as.numeric(as.character(AT1G53480$normexpr))
AT1G53480$Treatment <- relevel(AT1G53480$Treatment, "control")

AT5G22610<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("c",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT5G22610", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT5G22610", grep("Col-0",colnames(vsd))])))
colnames(AT5G22610) <- c("Genotype", "Treatment", "normexpr")
AT5G22610$normexpr <- as.numeric(as.character(AT5G22610$normexpr))
AT5G22610$Treatment <- relevel(AT5G22610$Treatment, "control")

AT1G75945<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT1G75945", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT1G75945", grep("Col-0",colnames(vsd))])))
colnames(AT1G75945) <- c("Genotype", "Treatment", "normexpr")
AT1G75945$normexpr <- as.numeric(as.character(AT1G75945$normexpr))
AT1G75945$Treatment <- relevel(AT1G75945$Treatment, "control")

AT3G48520<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT3G48520", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT3G48520", grep("Col-0",colnames(vsd))])))
colnames(AT3G48520) <- c("Genotype", "Treatment", "normexpr")
AT3G48520$normexpr <- as.numeric(as.character(AT3G48520$normexpr))
AT3G48520$Treatment <- relevel(AT3G48520$Treatment, "control")

AT1G75750<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4 hpi",3), rep("8 hpi",3)),
                   c(vsd[rownames = "AT1G75750", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT1G75750", grep("Col-0",colnames(vsd))])))
colnames(AT1G75750) <- c("Genotype", "Treatment", "normexpr")
AT1G75750$normexpr <- as.numeric(as.character(AT1G75750$normexpr))
AT1G75750$Treatment <- relevel(AT1G75750$Treatment, "control")


AT5G01900<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT5G01900", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT5G01900", grep("Col-0",colnames(vsd))])))
colnames(AT5G01900) <- c("Genotype", "Treatment", "normexpr")
AT5G01900$normexpr <- as.numeric(as.character(AT5G01900$normexpr))
AT5G01900$Treatment <- relevel(AT5G01900$Treatment, "control")

AT5G64810<- 
  data.frame(cbind(c(rep("pae9-2", 9),rep("Col-0", 9)),
                   c(rep("control",3), rep("4hpi",3), rep("8hpi",3)),
                   c(vsd[rownames = "AT5G64810", grep("pae9-2",colnames(vsd))],
                     vsd[rownames = "AT5G64810", grep("Col-0",colnames(vsd))])))
colnames(AT5G64810) <- c("Genotype", "Treatment", "normexpr")
AT5G64810$normexpr <- as.numeric(as.character(AT5G64810$normexpr))
AT5G64810$Treatment <- relevel(AT5G64810$Treatment, "control")

#' Plot clusters and timeplots
pdf("analysis/kallisto/Timeplots_clusters_DEG.pdf", width = 8, height = 6)
par(mfrow=c(3,3), mar = c(3,4,2,2))
#1
clusterplot(selGenoc_meanscale[nams[[3]],4:6],
            xlabels = NA, xaxt= 'n',
            colpal = "greens",las=1,ylim=c(-1.8,1.6),
            main = "Cluster I")
lineplot.CI(AT3G14440$Treatment,legend = FALSE,
            AT3G14440$normexpr, cex.main = 1.2, cex.axis = 1.2,
            group=AT3G14440$Genotype, xaxt = 'n',
            type="b", ylab="", lwd=2,
            err.width=0.002, ylim=c(1,4), 
            xlab = NA, col= c("gray40","chartreuse2"), pch=c(19,19),
            main="NCED3", bty = 'n') #NINE-CIS-EPOXYCAROTENOID DIOXYGENASE
legend(2,4,c("Col-0", expression(paste(italic("pae9-2")))),
       col=c("gray40","chartreuse2"),lty=c(2,1),pch=c(19,19),bty='n',cex=1.2)
lineplot.CI(AT3G48520$Treatment,legend = FALSE,x.leg = 2, y.leg = 2,
            AT3G48520$normexpr, cex.main = 1.2, xaxt = 'n',
            group=AT3G48520$Genotype,ylim=c(0.3,2), cex.axis = 1.2,
            type="b", ylab="", 
            err.width=0.002, lwd=2,
            xlab = "", col= c("gray40","chartreuse2"), pch=c(19,19),
            main="CYP94B3", bty = 'n')

#2
clusterplot(selGenoc_meanscale[nams[[2]],4:6],
            xlabels = NA, xaxt = 'n', cex = 1.2,
            colpal = "blues",las=2,ylim=c(-1.8,1.6),
            main = "Cluster II")
mtext("Z score", side=2, line=2, cex.axis = 1.2)
lineplot.CI(AT1G21120$Treatment, legend = F, x.leg = 2, y.leg = 1.8,
            AT1G21120$normexpr, cex.main = 1.2, cex.axis = 1.2,
            group=AT1G21120$Genotype,lwd=2, xaxt = 'n',
            type="b", ylab="Normalized expression", cex.lab = 1.2,
            err.width=0.002,#ylim=c(2,5), 
            xlab = NA,col= c("gray40","dodgerblue"), pch=c(19,19),
            main="IGMT2", bty = 'n') #INDOLE GLUCOSINOLATE O-METHYLTRANSFERASE 2
legend(2,1.8,c("Col-0", expression(paste(italic("pae9-2")))),
       col=c("gray40","dodgerblue"),lty=c(2,1),pch=c(19,19),bty='n',cex=1.2)
lineplot.CI(AT3G26830$Treatment,legend = FALSE,
            AT3G26830$normexpr, cex.main = 1.3, cex.axis = 1.2,
            group=AT3G26830$Genotype,xaxt = 'n',
            type="b", ylab="", lwd=2,
            err.width=0.002,ylim=c(1,3), 
            xlab = "",col= c("gray40","dodgerblue"), pch=c(19,19),
            main="PAD3", bty = 'n') #PHYTOALEXIN DEFICIENT 3
#3
clusterplot(selGenoc_meanscale[nams[[1]],4:6],
            xlabels = c("control", "4 hpi", "8 hpi"),
            ylab = "Z score", cex.axis = 1.2, 
            colpal = "reds",las=1,ylim=c(-1.8,1.6),
            main = "Cluster III")
lineplot.CI(AT1G53480$Treatment, legend = FALSE,
            AT1G53480$normexpr, cex.main = 1.2,cex.axis = 1.2,
            group=AT1G53480$Genotype,lwd=1.5,
            type="b", ylab="", 
            err.width=0.002,ylim=c(0,6.5), xlab = NA,col= c("gray40","red"), pch=c(19,19),
            main="MRD1", bty = 'n') #"MTO 1 RESPONDING DOWN 1
legend(2,6.5,c("Col-0", expression(paste(italic("pae9-2")))),
       col=c("gray40","red"),lty=c(2,1),pch=c(19,19),bty='n',cex=1.2)
lineplot.CI(AT1G75750$Treatment,legend = FALSE,x.leg = 2, y.leg = 2,
            AT1G75750$normexpr, cex.main = 1.2, cex.axis = 1.2,
            group=AT1G75750$Genotype,#ylim=c(0.4,2),
            type="b", ylab="", 
            err.width=0.002, lwd=2,
            xlab = "", col= c("gray40","red"), pch=c(19,19),
            main="GASA1", bty = 'n') #GA-RESPONSIVE GAST1
dev.off()

#---------------------------------------------------------
# Number of DEG
x <- read.csv("analysis/kallisto/DEgenes_allconstrast_raw.csv")
# control
genoc <- x[, c(2,3)]
genoc <- genoc[complete.cases(genoc),]
genoc <- subset(genoc, abs(genoc$LFC_Genoc)>=0.5)
genoc_down <- nrow(subset(genoc, genoc$LFC_Genoc < 0))
genoc_up <- nrow(subset(genoc, genoc$LFC_Genoc > 0))

# 4hpi
geno4 <- x[, c(2,5)]
geno4 <- geno4[complete.cases(geno4),]
geno4 <- subset(geno4, abs(geno4$LFC_Geno4)>=0.5)
geno4_down <- nrow(subset(geno4, geno4$LFC_Geno4 < 0))
geno4_up <- nrow(subset(geno4, geno4$LFC_Geno4 > 0))

# 8 hpi
geno8 <- x[, c(2,7)]
geno8 <- geno8[complete.cases(geno8),]
geno8 <- subset(geno8, abs(geno8$LFC_Geno8)>=0.5)
geno8_down <- nrow(subset(geno8, geno8$LFC_Geno8 < 0))
geno8_up <- nrow(subset(geno8, geno8$LFC_Geno8 > 0))

numDE <- data.frame(cbind(rbind(genoc_up, geno4_up, geno8_up),
                          rbind(genoc_down, geno4_down, geno8_down),
                          c("control", "4 hpi", "8 hpi")))
colnames(numDE) <- c("up", "down", "treatm")

numDE$up <- as.numeric(as.character(numDE$up))
numDE$down <- as.numeric(as.character(numDE$down))
numDE

pdf("analysis/kallisto/NumDEgenes_Geno.pdf", width=6, height=5)
par(mfrow = c(1,1), oma=c(0,0,0,0))
# total number of DE genes
x <- c(1:3, 3:1)
yup <- c(numDE[,1], rep(0,3))
ydown <- c(-numDE[,2], rep(0,3))
plot   (x, yup, ylim = c(-45, 40), xlim = c(1,3.2), type = "n", axes=F, xlab="", ylab="")
polygon(x, yup, col = rgb(1,0,0,0.6), border = "black")
polygon(x, ydown, col = "white", border = "black")
## Proportion of GO term based on number of significant enriched terms 
## in downregulated DEG inferred from BinGO (Maere et al. 2005)
# Up: no enriched GO terms. 
# Down:
# Biotic stress
y_biotic <- c((-34*9/12), (-36*10/17), -3, 0, 0, 0)
polygon(x, y_biotic, col = rgb(0,0,1,0.4), border = "black")
# Secondary metabolism
y_meta <- c(((-34*9/12)+(-34*1/12)), ((-36*10/17)+(-36*7/17)), -3, -3, (-36*10/17), (-34*9/12))
polygon(x, y_meta, col = rgb(0,0,1,0.6), border = "black")
# Other metabolism
y_oth <- c(((-34*9/12)+(-34*1/12)+(-34*2/12)), ((-36*10/17)+(-36*7/17)), -3, -3, ((-36*10/17)+(-36*7/17)), ((-34*9/12)+(-34*1/12)))
polygon(x, y_oth, col = rgb(0,0,1,0.8), border = "black")
# Add axes
axis(2, ylim=c(-45,40), at=c(-40,-20,0,20,40), labels=c(40,20,0,20,40),col="black", lwd=2, cex.axis=1)
mtext(2,text="# DEG",line=2.5,font=2, cex=1)
text(c(1,2,3), rep(-45,3),c("c", "4 hpi", "8 hpi"), cex=1)
abline(h=0, col="black", lwd=2)
# legend
legend(x=2,y=40, c("response to biotic stimulus", "secondary metabolism", "other"),
       fill=c(rgb(0,0,1,0.4),rgb(0,0,1,0.6),rgb(0,0,1,0.8)),bty='n')
dev.off()


