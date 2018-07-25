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
#suppressPackageStartupMessages(library(Glimma))
#suppressPackageStartupMessages(library(LSD))
#suppressPackageStartupMessages(library(Mfuzz))
#suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(VennDiagram))

#' Helper files
#suppressMessages(source("~/Git/UPSCb/src/R/densityPlot.R"))
suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))
#suppressMessages(source("~/Git/UPSCb/src/R/plot.multidensity.R"))

#' Load saved data
load("analysis/kallisto/DESeq-object_combtechrep_batch.rda")

#' Quickly VST it
vst <- varianceStabilizingTransformation(dds,blind=FALSE)

#' Extract the normalised counts
vsd <- assay(vst)

#' Which we remove, so that 0 means no expression
vsd <- vsd - min(vsd)


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


#' # CONTRAST group A: Col-0 c as reference ======================================
dds$Group <- relevel(dds$Group,"Col-0c")
#' Old design
design(dds) #~Batch + Treatment * Genotype
#' New design:
design(dds) <- ~Batch + Group
dds <- DESeq(dds)

#' Dispersion Estimation
#'
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' What contrasts are available
resultsNames(dds)
#' "Intercept"                  
#' "Batch_2_vs_1"               "Group_Col.04hpi_vs_Col.0c"  "Group_Col.08hpi_vs_Col.0c" 
#' "Group_pae9.24hpi_vs_Col.0c" "Group_pae9.28hpi_vs_Col.0c" "Group_pae9.2c_vs_Col.0c"  

#--------------------------------------
#' ## A1. Col-0c vs pae9-2c 
resGenoc <- results(dds,name="Group_pae9.2c_vs_Col.0c")

#' volcano plot
volcanoPlot(resGenoc)

#' a look at the results
# baseMean (linear expression, library size corrected)
# lfcSE (standard error on the log2 FC)
# padj == FDR (multiple testing correction, "P value adjusted")
head(resGenoc)

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selGenoc <- resGenoc$padj <= padj &
  abs(resGenoc$log2FoldChange) >= lfc &
  ! is.na(resGenoc$padj)

write.csv(resGenoc,file="analysis/kallisto/Group_pae9.2c_vs_Col.0c_batch.csv")
write(rownames(resGenoc[selGenoc,]),file="analysis/kallisto/Group_pae9.2c_vs_Col.0c_batch.txt")

#--------------------------------------
#' ## A2. Col-0c vs Col-04hpi
resColc4 <- results(dds,name="Group_Col.04hpi_vs_Col.0c")

#' volcano plot
volcanoPlot(resColc4)

#' a look at the results
# baseMean (linear expression, library size corrected)
# lfcSE (standard error on the log2 FC)
# padj == FDR (multiple testing correction, "P value adjusted")
head(resColc4)

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selColc4 <- resColc4$padj <= padj &
  abs(resColc4$log2FoldChange) >= lfc &
  ! is.na(resColc4$padj)

write.csv(resColc4,file="analysis/kallisto/Group_Col.04hpi_vs_Col.0c_batch.csv")
write(rownames(resColc4[selColc4,]),file="analysis/kallisto/Group_Col.04hpi_vs_Col.0c_batch.txt")

#--------------------------------------
#' ## A3. Col-0c vs Col-08hpi
resColc8 <- results(dds,name="Group_Col.08hpi_vs_Col.0c")

#' volcano plot
volcanoPlot(resColc8)

#' a look at the results
# baseMean (linear expression, library size corrected)
# lfcSE (standard error on the log2 FC)
# padj == FDR (multiple testing correction, "P value adjusted")
head(resColc8)

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selColc8 <- resColc8$padj <= padj &
  abs(resColc8$log2FoldChange) >= lfc &
  ! is.na(resColc8$padj)

write.csv(resColc8,file="analysis/kallisto/Group_Col.08hpi_vs_Col.0c_batch.csv")
write(rownames(resColc4[selColc8,]),file="analysis/kallisto/Group_Col.08hpi_vs_Col.0c_batch.txt")

#' # CONTRAST group B: Col-0 4hpi as reference ======================================
dds$Group <- relevel(dds$Group,"Col-04hpi")
#' New design:
design(dds) <- ~Batch + Group
dds <- DESeq(dds)

#' What contrasts are available
resultsNames(dds)
#[1] "Intercept"                     "Batch_2_vs_1"                  "Group_Col.0c_vs_Col.04hpi"    
#[4] "Group_Col.08hpi_vs_Col.04hpi"  "Group_pae9.24hpi_vs_Col.04hpi" "Group_pae9.28hpi_vs_Col.04hpi"
#[7] "Group_pae9.2c_vs_Col.04hpi"

#--------------------------------------
#' ## B1. Col-0 4hpi vs pae9-2 4hpi 
resGeno4 <- results(dds,name="Group_pae9.24hpi_vs_Col.04hpi")

#' volcano plot
volcanoPlot(resGeno4)

#' a look at the results
# baseMean (linear expression, library size corrected)
# lfcSE (standard error on the log2 FC)
# padj == FDR (multiple testing correction, "P value adjusted")
head(resGeno4)

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selGeno4 <- resGeno4$padj <= padj &
  abs(resGeno4$log2FoldChange) >= lfc &
  ! is.na(resGeno4$padj)

write.csv(resGeno4,file="analysis/kallisto/Group_pae9.24hpi_vs_Col.04hpi_batch.csv")
write(rownames(resGeno4[selGeno4,]),file="analysis/kallisto/Group_pae9.24hpi_vs_Col.04hpi_batch.txt")

#--------------------------------------
#' ## B2 Col-0 4hpi vs Col-0 8hpi 
resCol48 <- results(dds,name="Group_Col.08hpi_vs_Col.04hpi")

#' volcano plot
volcanoPlot(resCol48)

#' a look at the results
# baseMean (linear expression, library size corrected)
# lfcSE (standard error on the log2 FC)
# padj == FDR (multiple testing correction, "P value adjusted")
head(resCol48)

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selCol48 <- resCol48$padj <= padj &
  abs(resCol48$log2FoldChange) >= lfc &
  ! is.na(resCol48$padj)

write.csv(resCol48,file="analysis/kallisto/Group_Col.08hpi_vs_Col.04hpi_batch.csv")
write(rownames(resCol48[selCol48,]),file="analysis/kallisto/Group_Col.08hpi_vs_Col.04hpi_batch.txt")

#' # CONTRAST group C: Col-0 8hpi as reference ======================================
dds$Group <- relevel(dds$Group,"Col-08hpi")
#' New design:
design(dds) <- ~Batch + Group
dds <- DESeq(dds)

#' What contrasts are available
resultsNames(dds)
#[1] "Intercept"                     "Batch_2_vs_1"                  "Group_Col.04hpi_vs_Col.08hpi" 
#[4] "Group_Col.0c_vs_Col.08hpi"     "Group_pae9.24hpi_vs_Col.08hpi" "Group_pae9.28hpi_vs_Col.08hpi"
#[7] "Group_pae9.2c_vs_Col.08hpi" 

#--------------------------------------
#' ## C1. Col-0 8hpi vs pae9-2 8hpi 
resGeno8 <- results(dds,name="Group_pae9.28hpi_vs_Col.08hpi")

#' volcano plot
volcanoPlot(resGeno4)

#' a look at the results
# baseMean (linear expression, library size corrected)
# lfcSE (standard error on the log2 FC)
# padj == FDR (multiple testing correction, "P value adjusted")
head(resGeno8)

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selGeno8 <- resGeno8$padj <= padj &
  abs(resGeno8$log2FoldChange) >= lfc &
  ! is.na(resGeno8$padj)

write.csv(resGeno8,file="analysis/kallisto/Group_pae9.28hpi_vs_Col.08hpi_batch.csv")
write(rownames(resGeno8[selGeno8,]),file="analysis/kallisto/Group_pae9.28hpi_vs_Col.08hpi_batch.txt")


#' # CONTRAST group D: pae9-2 c as reference ======================================
dds$Group <- relevel(dds$Group,"pae9-2c")
#' New design:
design(dds) <- ~Batch + Group
dds <- DESeq(dds)

#' What contrasts are available
resultsNames(dds)
#[1] "Intercept"                   "Batch_2_vs_1"                "Group_Col.08hpi_vs_pae9.2c" 
#[4] "Group_Col.04hpi_vs_pae9.2c"  "Group_Col.0c_vs_pae9.2c"     "Group_pae9.24hpi_vs_pae9.2c"
#[7] "Group_pae9.28hpi_vs_pae9.2c" 

#--------------------------------------
#' ## D1. pae9-2 c vs pae9-2 4hpi 
resPaec4 <- results(dds,name="Group_pae9.24hpi_vs_pae9.2c")

#' volcano plot
volcanoPlot(resPaec4)

#' a look at the results
# baseMean (linear expression, library size corrected)
# lfcSE (standard error on the log2 FC)
# padj == FDR (multiple testing correction, "P value adjusted")
head(resPaec4)

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selPaec4 <- resPaec4$padj <= padj &
  abs(resPaec4$log2FoldChange) >= lfc &
  ! is.na(resPaec4$padj)

write.csv(resPaec4,file="analysis/kallisto/Group_pae9.24hpi_vs_pae9.2c_batch.csv")
write(rownames(resPaec4[selPaec4,]),file="analysis/kallisto/Group_pae9.24hpi_vs_pae9.2c_batch.txt")

#--------------------------------------
#' ## D2. pae9-2 c vs pae9-2 8hpi 
resPaec8 <- results(dds,name="Group_pae9.28hpi_vs_pae9.2c")

#' volcano plot
volcanoPlot(resPaec8)

#' a look at the results
# baseMean (linear expression, library size corrected)
# lfcSE (standard error on the log2 FC)
# padj == FDR (multiple testing correction, "P value adjusted")
head(resPaec8)

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selPaec8 <- resPaec8$padj <= padj &
  abs(resPaec8$log2FoldChange) >= lfc &
  ! is.na(resPaec8$padj)

write.csv(resPaec8,file="analysis/kallisto/Group_pae9.28hpi_vs_pae9.2c_batch.csv")
write(rownames(resPaec8[selPaec8,]),file="analysis/kallisto/Group_pae9.28hpi_vs_pae9.2c_batch.txt")


#' # CONTRAST group E: pae9-2 4hpi as reference ======================================
dds$Group <- relevel(dds$Group,"pae9-24hpi")
#' New design:
design(dds) <- ~Batch + Group
dds <- DESeq(dds)

#' What contrasts are available
resultsNames(dds)
#[1] "Intercept"                      "Batch_2_vs_1"                   "Group_pae9.2c_vs_pae9.24hpi"   
#[4] "Group_Col.08hpi_vs_pae9.24hpi"  "Group_Col.04hpi_vs_pae9.24hpi"  "Group_Col.0c_vs_pae9.24hpi"    
#[7] "Group_pae9.28hpi_vs_pae9.24hpi"

#--------------------------------------
#' ## E1. pae9-2 4 hpi vs pae9-2 8hpi 
resPae48 <- results(dds,name="Group_pae9.28hpi_vs_pae9.24hpi")

#' volcano plot
volcanoPlot(resPae48)

#' a look at the results
# baseMean (linear expression, library size corrected)
# lfcSE (standard error on the log2 FC)
# padj == FDR (multiple testing correction, "P value adjusted")
head(resPae48)

#' Select the genes that are DE
#' Schurch et al.,RNA, 2016
padj <- 0.01
lfc <- 0.5
selPae48 <- resPae48$padj <= padj &
  abs(resPae48$log2FoldChange) >= lfc &
  ! is.na(resPae48$padj)

write.csv(resPae48,file="analysis/kallisto/Group_pae9.28hpi_vs_pae9.24hpi_batch.csv")
write(rownames(resPae48[selPae48,]),file="analysis/kallisto/Group_pae9.28hpi_vs_pae9.24hpi_batch.txt")

#===========================================================================================

#' ## Heatmaps significant genes-------------------------------

#' Attach sample information
colnames(vsd)
samples <- read.delim("/mnt/picea/projects/arabidopsis/balbrectsen/arabidopsis-aphids/analysis/kallisto/Sample_names.txt")
samples <- samples[grep("L001", samples$SciLifeID),]
samples <- samples[order(samples$SampleID),]
samples$lab <- paste(samples$Genotype, samples$Treatment, sep = "_")
colnames(vsd) <- samples$lab

#' Col-0 within (sign genes between treatment)
heatmap.2(vsd[rownames(c(resColc4[selColc4,], resColc8[selColc8,], resCol48[selCol48,])), #significant genes 
              grep("Col-0",colnames(vsd))], #samples
          trace="none",
          col=hpal,scale="row",
          labCol = colnames(vsd)[grep("Col-0", colnames(vsd))],
          labRow = rownames(c(resColc4[selColc4,], resColc8[selColc8,], resCol48[selCol48,])))

#' pae9-2 within (sign genes between treatment)
heatmap.2(vsd[rownames(resPaec8[selPaec8,], c(resPaec4[selPaec4,], resPae48[selPae48,])), #significant genes 
              grep("pae9-2",colnames(vsd))], #samples
          trace="none",
          col=hpal,scale="row",
          labCol = colnames(vsd)[grep("pae9-2", colnames(vsd))],
          labRow = rownames(resPaec8[selPaec8,], c(resPaec4[selPaec4,], resPae48[selPae48,])))

#' Significant genes between genotypes (only c)
heatmap.2(vsd[rownames(resGenoc[selGenoc,]), #significant genes
              -grep("hpi",colnames(vsd))], #samples 
          trace="none",
          col=hpal,scale="row",
          labCol = colnames(vsd)[-grep("hpi", colnames(vsd))],
          labRow = rownames(resGenoc[selGenoc,]))

#' Significant genes between genotypes (only 4)
heatmap.2(vsd[rownames(resGeno4[selGeno4,]), #significant genes
              grep("4hpi",colnames(vsd))], #samples 
          trace="none",
          col=hpal,scale="row",
          labCol = colnames(vsd)[grep("4hpi", colnames(vsd))],
          labRow = rownames(resGeno4[selGeno4,]))


#' Significant genes between genotypes (only 8)
heatmap.2(vsd[rownames(resGeno8[selGeno8,]), #significant genes
              grep("8hpi",colnames(vsd))],  #samples
          trace="none",
          col=hpal,scale="row",
          labCol = colnames(vsd)[grep("8hpi", colnames(vsd))],
          labRow = rownames(resGeno8[selGeno8,]))

#' Significant genes between genotypes (within all treatments)
heatmap.2(vsd[rownames(resGenoc[selGenoc,], resGeno4[selGeno4,], resGeno8[selGeno8,]),], #significant genes 
          trace="none",
          col=hpal,scale="row",
          labCol = colnames(vsd),
          labRow = rownames(resGenoc[selGenoc,], c(resGeno4[selGeno4,], resGeno8[selGeno8,])))

#' ## Gene lists-------------------------------------
#' DE genes overview
#' a. DE between Genotype within Treatment
Genoc <- cbind(rownames(resGenoc[selGenoc,]), data.frame(resGenoc[selGenoc,c(2,6)]))
colnames(Genoc) <- c("GeneID", "LFC_Genoc", "Padj_Genoc")
Geno4 <- cbind(rownames(resGeno4[selGeno4,]),data.frame(resGeno4[selGeno4,c(2,6)]))
colnames(Geno4) <- c("GeneID", "LFC_Geno4", "Padj_Geno4")
Geno8 <- cbind(rownames(resGeno8[selGeno8,]),data.frame(resGeno8[selGeno8,c(2,6)]))
colnames(Geno8) <- c("GeneID", "LFC_Geno8", "Padj_Geno8")
#' b. DE between Treatments within Col-0
Colc4 <- cbind(rownames(resColc4[selColc4,]), data.frame(resColc4[selColc4,c(2,6)]))
colnames(Colc4) <- c("GeneID", "LFC_Colc4", "Padj_Colc4")
Colc8 <- cbind(rownames(resColc8[selColc8,]), data.frame(resColc8[selColc8,c(2,6)]))
colnames(Colc8) <- c("GeneID", "LFC_Colc8", "Padj_Colc8")
#' c. DE between Treatments within pae9-2
Paec4 <- cbind(rownames(resPaec4[selPaec4,]), data.frame(resPaec4[selPaec4,c(2,6)]))
colnames(Paec4) <- c("GeneID", "LFC_Paec4", "Padj_Paec4")
Paec8 <- cbind(rownames(resPaec8[selPaec8,]), data.frame(resPaec8[selPaec8,c(2,6)]))
colnames(Paec8) <- c("GeneID", "LFC_Paec8", "Padj_Paec8")
#' Merge
deg <- merge(Genoc, Geno4, by = "GeneID", all = TRUE)
deg <- merge(deg, Geno8, by = "GeneID", all = TRUE)
deg <- merge(deg, Colc4, by = "GeneID", all = TRUE)
deg <- merge(deg, Colc8, by = "GeneID", all = TRUE)
deg <- merge(deg, Paec4, by = "GeneID", all = TRUE)
deg <- merge(deg, Paec8, by = "GeneID", all = TRUE)
#' Get GO annotations TAIR
suppressPackageStartupMessages(library(org.At.tair.db))
deg2 <- data.frame(GeneID=as.character(deg$GeneID),
                   SYMBOL=mapIds(org.At.tair.db,keys=as.character(deg$GeneID),"SYMBOL","TAIR"),
                   GENENAME=mapIds(org.At.tair.db,keys=as.character(deg$GeneID),"GENENAME","TAIR"))
deg3 <- merge(deg, deg2, by = "GeneID")
write.csv(deg3, "analysis/kallisto/DEgenes_allconstrast_raw.csv")

#' # 2. DE genes within Col-0 versus Geno4 and Geno8
#' Are the aphid-induced genes also DE in pae9-2? Does pae9-2 have the full wildtype aphid response?
Colc4 <- cbind(rownames(resColc4[selColc4,]), data.frame(resColc4[selColc4,c(2,6)]))
colnames(Colc4) <- c("Genes", "LFC_Colc4", "Padj")
Colc8 <- cbind(rownames(resColc8[selColc8,]), data.frame(resColc8[selColc8,c(2,6)]))
colnames(Colc8) <- c("Genes", "LFC_Colc8", "Padj")
Col48 <- cbind(rownames(resCol48[selCol48,]), data.frame(resCol48[selCol48,c(2,6)]))
colnames(Col48) <- c("Genes", "LFC_Col48", "Padj")
colnames(Genoc) <- c("Genes", "LFC_c", "Padj")
colnames(Geno4) <- c("Genes", "LFC_4", "Padj")
colnames(Geno8) <- c("Genes", "LFC_8", "Padj")
Colx <- merge(Colc4, Colc8, by = "Genes", all = TRUE)
Colx <- merge(Colx, Col48, by = "Genes", all = TRUE)
Colx <- merge(Colx, Genoc, by = "Genes", all = TRUE)
Colx <- merge(Colx, Geno4, by = "Genes", all.x = TRUE, all.y = FALSE)
Colx <- merge(Colx, Geno8, by = "Genes", all.x = TRUE, all.y = FALSE)
Colx <- merge(Colx, tair, by = "Genes", all.x = TRUE, all.y = FALSE)
write.csv(Colx, "analysis/kallisto/DEgenes_Col-0_and_Geno-effect.csv")


#' # 3. DE genes within Pae-0 versus within Col-0
#' What is the pae9 * aphid interaction? What genes are aphid-affected in pae9-2 but not in Col-0?
Paec4 <- cbind(rownames(resPaec4[selPaec4,]), data.frame(resPaec4[selPaec4,c(2,6)]))
Paec4 <- Paec4[,c(1,2)]
colnames(Paec4) <- c("Genes", "LFC_Paec4")
Paec8 <- cbind(rownames(resPaec8[selPaec8,]), data.frame(resPaec8[selPaec8,c(2,6)]))
Paec8 <- Paec8[,c(1,2)]
colnames(Paec8) <- c("Genes", "LFC_Paec8")
Pae48 <- cbind(rownames(resPae48[selPae48,]), data.frame(resPae48[selPae48,c(2,6)]))
Pae48 <- Pae48[,c(1,2)]
colnames(Pae48) <- c("Genes", "LFC_Pae48")
Paex <- merge(Paec4, Paec8, by = "Genes", all = TRUE)
Paex <- merge(Paex, Pae48, by = "Genes", all = TRUE)
Paex <- merge(Paex, Colc4, by = "Genes", all.x = TRUE, all.y = FALSE)
Paex <- merge(Paex, Colc8, by = "Genes", all.x = TRUE, all.y = FALSE)
Paex <- merge(Paex, Col48, by = "Genes", all.x = TRUE, all.y = FALSE)
Paex <- merge(Paex, tair, by = "Genes", all.x = TRUE, all.y = FALSE)
write.csv(Paex, "analysis/kallisto/DEgenes_within_Pae9-2_Col-0.csv")

#' Number of DE genes
numDEgenes <- 
  data.frame(cbind("Contrast" =
                     c("Col-0_c4", "Col-0_c8", "Col-0_48",
                       "Pae9-2_c4", "Pae9-2_c8", "Pae9-2_48",
                       "Geno_c", "Geno_4", "Geno_8"),
                   "numDE" = 
                     c(nrow(resColc4[selColc4,]), nrow(resColc8[selColc8,]), nrow(resCol48[selCol48,]),  
                       nrow(resPaec4[selPaec4,]), nrow(resPaec8[selPaec8,]), nrow(resPae48[selPae48,]),
                       nrow(resGenoc[selGenoc,]), nrow(resGeno4[selGeno4,]),nrow(resGeno8[selGeno8,]))))
numDEgenes
write.csv(numDEgenes, "analysis/kallisto/numDEgenes.csv", row.names = FALSE)

#' Number of DE genes split in UP and DOWN
upColc4 <- sum(ifelse(resColc4[selColc4,2] > 0, 1, 0))
downColc4 <- sum(ifelse(resColc4[selColc4,2] < 0, 1, 0))
upColc8 <- sum(ifelse(resColc8[selColc8,2] > 0, 1, 0))
downColc8 <- sum(ifelse(resColc8[selColc8,2] < 0, 1, 0))
upPaec4 <- sum(ifelse(resPaec4[selPaec4,2] > 0, 1, 0))
downPaec4 <- sum(ifelse(resPaec4[selPaec4,2] < 0, 1, 0))
upPaec8 <- sum(ifelse(resPaec8[selPaec8,2] > 0, 1, 0))
downPaec8 <- sum(ifelse(resPaec8[selPaec8,2] < 0, 1, 0))

numUpDown <- data.frame(cbind(rbind(upColc4, upColc8, upPaec4, upPaec8),
                              rbind(downColc4, downColc8, downPaec4, downPaec8), 
                              c(rep("Col-0", 2), rep("pae9-2",2))))
colnames(numUpDown) <- c("up", "down", "Genotype")

pdf("analysis/kallisto/numGenes_upanddown.pdf", width = 3, height = 6)
barplot(as.numeric(as.character(numUpDown[,1])),beside=T,ylim=c(-1000,600),
        yaxt="n", xaxt="n", col="red")
barplot(-as.numeric(as.character(numUpDown[,2])),beside=T,ylim=c(-1000,600),
        yaxt="n",col="blue", las = 2,
        names.arg=c("4 hpi","8 hpi","4 hpi", "8 hpi"), add=T)
axis(2,ylim=c(-1000,500),labels=c(-1000, -500, 0 , 500),at=c(-1000, -500, 0 , 500),col="black")
mtext(2,text="Number of DE genes",line=3,font=2)
mtext(3,text="Col-0",line=0.2,font=1, adj = 0.2)
mtext(3,text=expression(paste(italic("pae9-2"))),line=0,font=2, adj = 0.9)
segments(x0 = 0.4, y0 = 550, x1 = 2.2, y1 = 550, lwd = 2)
segments(x0 = 3, y0 = 550, x1 = 4.8, y1 = 550, lwd = 2)
dev.off()

# VennDiagram sign genes between genotypes
plot.new()

grid.draw(venn.diagram(list(rownames(resGenoc[selGenoc,]),
                            rownames(resGeno4[selGeno4,]),
                            rownames(resGeno8[selGeno8,])),
                       filename=NULL,category.names = c("c","4 hpi","8 hpi"),
                       fill=pal[1:3],
                       resolution=500))


plot.new()
grid.draw(venn.diagram(list(rownames(resColc4[selColc4,]),
                            rownames(resColc8[selColc8,]),
                            rownames(resPaec4[selPaec4,]),
                            rownames(resPaec8[selPaec8,])),
                       filename=NULL,
                       category.names = c("Col-0_4hpi","Col-0_8hpi", "pae9-2_4hpi", "pae9-2_8hpi"),
                       fill=pal[1:4],
                       resolution=500))


# list common DE genes
common.Genoc4 <- intersect(rownames(resGenoc[selGenoc,]),rownames(resGeno4[selGeno4,]))
common.Genoc4
common.Genoc8 <- intersect(rownames(resGenoc[selGenoc,]),rownames(resGeno8[selGeno8,]))
common.Genoc8
common.Geno48 <- intersect(rownames(resGeno4[selGeno4,]),rownames(resGeno8[selGeno8,]))
common.Geno48
c4 <- data.frame(cbind("Genes" = common.Genoc4, "c_4hpi" = rep(TRUE, length(common.Genoc4))))
c8 <- data.frame(cbind("Genes" = common.Genoc8, "c_8hpi" = rep(TRUE, length(common.Genoc8))))
c48 <- data.frame(cbind("Genes" = common.Geno48, "4hpi_8hpi" = rep(TRUE, length(common.Geno48))))
commongeno <- merge(c4, c8, by = "Genes", all = TRUE)
commongeno <- merge(commongeno, c48, by = "Genes", all = TRUE)
commongeno[is.na(commongeno),]
write.csv(commongeno,file="analysis/kallisto/CommonDE_Genotype.csv")

# list all unique DE genes between treatments 
diff.Genoc4 <- union(rownames(resGenoc[selGenoc,]),rownames(resGeno4[selGeno4,]))
diff.Genoc4
diff.Genoc8 <- union(rownames(resGenoc[selGenoc,]),rownames(resGeno8[selGeno8,]))
diff.Genoc8
diff.Geno48 <- union(rownames(resGeno4[selGeno4,]),rownames(resGeno8[selGeno8,]))
diff.Geno48
c4 <- data.frame(cbind("Genes" = diff.Genoc4, "c_4hpi" = rep(TRUE, length(diff.Genoc4))))
c8 <- data.frame(cbind("Genes" = diff.Genoc8, "c_8hpi" = rep(TRUE, length(diff.Genoc8))))
c48 <- data.frame(cbind("Genes" = diff.Geno48, "4hpi_8hpi" = rep(TRUE, length(diff.Geno48))))
diffgeno <- merge(c4, c8, by = "Genes", all = TRUE)
diffgeno <- merge(diffgeno, c48, by = "Genes", all = TRUE)
write.csv(diffgeno,file="analysis/kallisto/diffDE_Genotype.csv")

