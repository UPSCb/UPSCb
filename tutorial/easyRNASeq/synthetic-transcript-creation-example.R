#' ---
#' title: "Synthetic transcripts generation example"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#'
#' # Setup
#' Load the libraries
library(easyRNASeq)
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(genomeIntervals))
library(pander)

#' Source an helper file
source("https://microasp.upsc.se/root/upscb-public/raw/master/src/R/createSyntheticTranscripts.R")

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
  filename="~/Box Sync/Projects/easyRNASeq/Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz",
  input="gtf",
  feature="transcript",
  output="GRanges")

#' ## Export
#' Save the object for later re-use
save(gAnnot, file="Drosophila_melanogaster.BDGP5.77-synthetic-transcripts.rda")

#' ## Summarization
#' ### Set the params
param <- RnaSeqParam(annotParam=AnnotParam(datasource=gAnnot),
                     bamParam=BamParam(paired=FALSE))

#' ### Get the BAM files
bamFiles <- getBamFileList(filenames=
                             dir(system.file(package="RnaSeqTutorial","extdata"),
                                 pattern="[A,T].*.bam$",
                                 full.names=TRUE))

#' ### Run
sexp <- simpleRNASeq(bamFiles=bamFiles,param=param,verbose=TRUE)

#' ### Check
pander(colSums(assay(sexp)))

#' # Session Info
sessionInfo()
