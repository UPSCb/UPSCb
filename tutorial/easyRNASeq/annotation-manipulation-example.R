#' ---
#' title: "Annotation manipulation example"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#'
#' # Aim
#' The aim is to get the sequence names in the annotation and BAM files' header
#' in sync to overcome the easyRNASeq error
#' ```{r easyRNASeq error, echo=FALSE}
#'  message("There is no common genomic references between your BAM files and the provided annotation. Fix one or the other.")
#' ```
#'
#' # Setup
#' Load the libraries
library(easyRNASeq)
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(genomeIntervals))
library(pander)

#' Source an helper file
source("https://microasp.upsc.se/root/upscb-public/raw/master/src/R/createSyntheticTranscripts.R")

#' Download the annotation file (gtf) EnsEMBL:
#'
#' ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.80.gtf.gz
#'
#' Note that here, I use the local copy I have downloaded and copied into my
#' current directory (where this script is).

#' # Overview
#'
#' 1. First, we will process the annotation file to create synthetic transcripts
#'
#' 2. Second, we will edit the genomic reference names to match those in the
#' BAM files

#' # Process
#' ## Synthetic transcripts creation
#' This function takes a gtf or gff3 _filename_ as input.
#'
#' The _input_ parameter defines the file format (default to gff3).
#'
#' The _feature_ parameter defines which feature to look for in the provided
#' file. Commonly mRNA for gff3 and transcript for gtf. It defaults to mRNA.
#' Several parameter can ge given as argument.
#'
#' The _output_ paramter defines the type of object that is returned.
#' It can generate a **Genome_intervals** or a __GRanges__ class of objects.
#' The former can be saved as a gff3 using the writeGff3 function from the
#' genomeIntervals package (loaded). The latter can be saved as an RData object
#' and/or be used directly in the construction of an AnnotParam.
gAnnot <- createSyntheticTranscripts(
  filename="Homo_sapiens.GRCh38.80.gtf.gz",
  input="gtf",
  feature="transcript",
  output="GRanges")

#' This created the synthetic transcripts. We now change the sequence names
#' to match those in the BAM files (which are prepended with a 'chr' string).
seqlevels(gAnnot) <- paste("chr",seqlevels(gAnnot),sep="")

#' We also 'rescue' the mitochondria. The other sequence names (haplotypes)
#' would require a manual curation to be added. Right now, we should have the
#' 22 autosomes, the sexual chromosomes and the mitochondria.
seqlevels(gAnnot)[seqlevels(gAnnot)=="chrMT"] <- "chrM"

#' ## Export
#' Save the object for later re-use
save(gAnnot, file="Homo_sapiens.GRCh38.80-synthetic-transcripts.rda")

#' ## Summarization
#' ### Set the params
#' Here it has been tricky. The flag in the BAM files says that the
#' record are Paired-end reads, but a closer look revealed that they
#' are not. An example record is:
#'
#' FCC64Y1ACXX:1:2311:16630:3201#GCCAATAT  81      chr10   60021   0       49M     *       0       0       GCATCGGGGTGCTCTGGTTTTGTTGTTGTTATTTCTGAATGACATTTAC       hiihhiiiiiiiiiiiiihdhiiiihhfiiihihhheggggeeeeebb_       NM:i:0  MD:Z:49
#'
#' where the 7 to 9th columns represent the mate alignment, which is in all
#' occurences empty (* for the sequence name, 0 for the pos and 0 for the
#' insert size)
#'
#' To remediate that, we set the argument paired to 'FALSE' and we need
#' to call 'simpleRNASeq' with the 'override' argument set to 'TRUE' to
#' avoid the parameter auto-detection.
param <- RnaSeqParam(annotParam=AnnotParam(datasource=gAnnot),
                     bamParam=BamParam(paired=FALSE),
                     countBy="gene")

#' ### Get the BAM files
#' Here again, I use the provided BAM files, which are now in my script
#' directory
bamFiles <- getBamFileList(
  filenames=dir(".",pattern=".*.bam$",full.names=TRUE))

#' ### Run
sexp <- simpleRNASeq(bamFiles=bamFiles,param=param,override=TRUE,verbose=TRUE)

#' ### Check
pander(colSums(assay(sexp)))

#' # Session Info
sessionInfo()
