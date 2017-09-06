#' ---
#' title: "P. tremula genome self-self Blast results parsing"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---

#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/aspseq/tremula_tremuloides_comp/BLAST-subset")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/tremula_tremuloides_comp/BLAST-subset")
#' ```
# #' Load libraries
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(VennDiagram))

#' Source helpers
source("~/Git/UPSCb/src/R/blastUtilities.R")
source('~/Git/UPSCb/src/R/plot.multidensity.R')

#' Create a palette
pal <- brewer.pal(5,"Dark2")

#' # Initial attempt
#' To enhance the speed, I have split the original file in subset of 
#' 10000 lines. The issue is though, that some scaffold HSPs will be
#' in separate files and need to be recalculated afterwards
PotraBLAST <- mclapply(dir(".",pattern=".*\\.blt"),
  readBlast,ignoreSelf = TRUE, 
                        format=c("query.id",
                                 "subject.id",
                                 "percent.identity",
                                 "alignment.length",
                                 "mismatches",
                                 "gap.opening",
                                 "query.start",
                                 "query.end",
                                 "subject.start",
                                 "subject.end",
                                 "e.value",
                                 "bit.score",
                                 "query.length",
                                 "subject.length"),
                        plot = FALSE,mc.cores=4L)

#' # Process all the chunks
#' And re-calculate the scaffold values which occured when the original
#' file was split
dat <- do.call(rbind,mclapply((1:length(PotraBLAST)),function(i,bl){
  
  ## get the scaffolds in the previous chunk
  if(i==1){
    p.scfs = ""
  } else {
    p.scfs <- unique(bl[[i-1]][["df"]]$query.id)
  }
  
  ## get the scaffolds in the chunk
  scfs <- unique(bl[[i]][["df"]]$query.id)
  
  ## get the scaffolds in the next chunk
  if(i==length(bl)){
    n.scfs <- ""
  } else {
    n.scfs <- unique(bl[[i+1]][["df"]]$query.id)
  }
  
  ## keep the scf we want
  scf <- c(setdiff(scfs,p.scfs),intersect(scfs,n.scfs))
  
  ## get the data
  if(i==length(bl)){
    dat <- bl[[i]][["df"]][bl[[i]][["df"]]$query.id %in% scf,]
  } else {
    dat <- rbind(bl[[i]][["df"]][bl[[i]][["df"]]$query.id %in% scf,],
                 bl[[i+1]][["df"]][bl[[i+1]][["df"]]$query.id %in% scf,])
  }
  
  ## identify the doublon
  pos <- which(duplicated(dat[,c("query.id","subject.id")]))
  if(length(pos)){
    q <- dat$query.id[pos]
    s <- dat$subject.id[pos]
    dat <- dat[- which(dat$query.id== q & dat$subject.id==s), ]
    
    ## recalculate the cumulative coverage
    hsp <- rbind(bl[[i]][["blf"]],bl[[i+1]][["blf"]])
    hsp <- hsp[hsp$query.id==q & hsp$subject.id==s,]
    
    dat <- rbind(dat,data.frame(query.id=q,
                     subject.id=s,
                     query.cum.cov=sum(width(reduce(
                       IRanges(
                         start=ifelse(hsp$query.start>hsp$query.end,hsp$query.end,hsp$query.start),
                         end=ifelse(hsp$query.start<hsp$query.end,hsp$query.end,hsp$query.start))
                     )))/hsp$query.length[1],
                     subject.cum.cov=sum(width(reduce(
                       IRanges(
                         start=ifelse(hsp$subject.start>hsp$subject.end,hsp$subject.end,hsp$subject.start),
                         end=ifelse(hsp$subject.start<hsp$subject.end,hsp$subject.end,hsp$subject.start))
                       )))/hsp$subject.length[1],stringsAsFactors=FALSE))
  }
  
  return(dat)  
},PotraBLAST,mc.cores=4))

#' The obtained object contains all reported hits per scaffold (its cumulative
#' coverage).
#' We sort it first by by query.coverage - this minimally matters as hit are 
#' probably duplicated, e.g. it is likely that 
#' scf1, scf2, 1, 0.5 will also appear as scf2, scf1, 0.5, 1 if scaffolds are
#' very similar.
s.dat <- dat[order(dat$query.cum.cov,decreasing=TRUE),]
f.dat <- s.dat[match(unique(s.dat$query.id),s.dat$query.id),]

#' However, to avoid loosing information if some scaffolds are not the best hit
#' or are not reported, we check the assumption above.
seq.annot <- read.delim(
  "/mnt/picea/storage/reference/Populus-tremula/v1.0/fasta/Potra01-genome.fa.fai",
  header=FALSE)[,1:2]
colnames(seq.annot) <- c("scf","len")

#' What are the combination of subject-query?
sprintf("Out of %s unique scaffolds, %s are present as either query or subject",
  length(union(s.dat$query.id,s.dat$subject.id)),
  length(intersect(s.dat$query.id,s.dat$subject.id)))

qs <- paste(s.dat$query.id,s.dat$subject.id,sep="-")
sq <- paste(s.dat$subject.id,s.dat$query.id,sep="-")

#' Plotting takes a long time and the results are symmetric
#' anyway
sprintf("There are %s common query-subject and %s unique to either set",
        sum(qs %in% sq),
        sum(!qs %in% sq))

#' Get them all, re-order, and keep the best hit.
#' We need to create a tmp object for reverting the data.frame as
#' rbind on data.frames uses the column names to do the binding, not their 
#' position, unlike applying rbind to a matrix.
r.dat <- s.dat[!sq %in% qs,c(2,1,4,3)]
colnames(r.dat) <- colnames(s.dat)
f.dat <- rbind(s.dat[qs %in% sq,],
               s.dat[!qs %in% sq,],
               r.dat)

f.dat <- f.dat[order(f.dat$query.cum.cov,decreasing=TRUE),]
f.dat <- f.dat[match(unique(f.dat$query.id),f.dat$query.id),]

#' Now have a look at the query cumulative coverage distribution
plot(density(f.dat$query.cum.cov),col=pal[1],lwd=2,
     main="Query Cumulative coverage")

comparisonplot(f.dat$query.cum.cov,f.dat$subject.cum.cov,
               xlab="query cum. cov.",ylab="subject cum. cov.",
               main="")

heatscatter(f.dat$query.cum.cov,f.dat$subject.cum.cov,
            xlab="query cum. cov.",ylab="subject cum. cov.",
            main="")

abline(h=0.97,v=0.97,col=2,lty=2,lwd=2)

#' This is more as expected, a lot of small sequences are almost fully covered 
#' in big sequences, but there's otherwise a wide distribution. 
sel <- f.dat$query.cum.cov == 1 & f.dat$subject.cum.cov == 1
sprintf("There are %s scaffolds that appear to be fully percent redundant",
        sum(sel))

#' Let's look at these in a pairwise alignment fashion. Subset the necessary data
red <- f.dat[sel,]

#' Then load the sequences
seq <- readDNAStringSet("/mnt/picea/storage/reference/Populus-tremula/v1.0/fasta/Potra01-genome.fa")

#' And perform a set if pairwise alignments. This shows the flaw of this approach
#' which has ignored the percentage of identity of the HSPs. It nevertheless
#' revealed that a lot of sequences are contained within other sequence with
#' an high level of identity (>80%). This is a good sign that the assembly
#' managed to integrate a reasonable number of somewhat repetitive elements.
pA <- pairwiseAlignment(reverseComplement(seq[[red[1,2]]]),seq[[red[1,1]]])
length(pattern(pA))
as.character(pattern(pA))
as.character(subject(pA))
nchar(pA)
length(mismatch(subject(pA))[[1]])
length(mismatch(pattern(pA))[[1]])
indel(pA)

#' # Second attempt
#' This attempt is based on filtering the HSPs based on their percentage of 
#' identity and then construct the cumulative coverage as done in the first
#' attempt.
blt <- do.call(rbind,mclapply(dir(".",pattern=".*\\.blt"),
                function(fil){
                  readBlast(fil,ignoreSelf = TRUE, 
                            format=c("query.id",
                                     "subject.id",
                                     "percent.identity",
                                     "alignment.length",
                                     "mismatches",
                                     "gap.opening",
                                     "query.start",
                                     "query.end",
                                     "subject.start",
                                     "subject.end",
                                     "e.value",
                                     "bit.score",
                                     "query.length",
                                     "subject.length"),
                            plot = FALSE)$blf},mc.cores=4L))

#' Reading the chunks is fast but combining the final results takes a while,
#' it might be worth to integrate the percent identity filtering in the 
#' blastUtility.R helper.
#' ```{r devnull, echo=FALSE, eval=FALSE}
#' ```
#' The overall density of the percentage identity of all the HSPs is as
#' follows. It is really interesting to observe these very defined peaks
#' towards the right end of the graph. We have a number of perfect hits, some
#' more hits around 97% identity (possible haplotypes), a tinier peak at 95% 
#' and a shoulder around 93%. The bulge of the remaining hits centers around 
#' 88%. The peak intervals are surprisingly constant and agrees well with the 
#' estimated heterozygosity rate. Concerning the peaks lower than 97%, this could
#' lead to wild hypotheses :-)
plot(density(blt$percent.identity),main="HSPs percentage identity",
     col=pal[1],lwd=2)
abline(v=c(95:100),lty=2,col="grey")

#' Next we define a function (which should be integrated in the blastUtility,R) 
#' that filters HSPs based on percent identity and calculate and sort the
#' obtained cumulative coverage
getCumulativeCoverage <- function(blt,perc.ident=95){
  blf <- blt[blt$percent.identity>=perc.ident,]
  
  ids <- paste(blf$query.id,blf$subject.id,sep="+")
  suids <- sort(unique(ids))
  
  df <- data.frame(query.id=sub("\\+.*","",suids),
                   subject.id=sub(".*\\+","",suids),
                   query.cum.cov=sum(width(reduce(split(
                     IRanges(
                       start=ifelse(blf$query.start>blf$query.end,blf$query.end,blf$query.start),
                       end=ifelse(blf$query.start<blf$query.end,blf$query.end,blf$query.start)),
                     ids))))/blf$query.length[match(sub("\\+.*","",suids),blf$query.id)],
                   subject.cum.cov=sum(width(reduce(split(
                     IRanges(
                       start=ifelse(blf$subject.start>blf$subject.end,blf$subject.end,blf$subject.start),
                       end=ifelse(blf$subject.start<blf$subject.end,blf$subject.end,blf$subject.start)),
                     ids))))/blf$subject.length[match(sub(".*\\+","",suids),blf$subject.id)],
                   stringsAsFactors=FALSE)
  
  return(df[order(df$query.cum.cov,df$subject.cum.cov,decreasing=TRUE),])
}

#' Now iteratively get the cumulative coverage
res <- mclapply(seq(95,100,1),function(p,blt){
  return(getCumulativeCoverage(blt,p))  
},blt,mc.cores=6L)
names(res) <- paste("Perc","Ident",seq(95,100,1),sep=".")

#' And have a look at the number of scaffold pairs linked by HSPs
barplot(sapply(res,nrow),main="# of linked scaffold pairs")

#' Have a look at the relationship query - subject coverage
dev.null <- sapply(1:length(res),function(i,res){
  re <- res[[i]]
  comparisonplot(re$query.cum.cov,
                 re$subject.cum.cov,
                 xlab="query cumulative coverage",
                 ylab="subject cumulative coverage",
                 main=names(res)[i])
},res)

#' Have a look at number of unique scaffolds involved
scfs <- lapply(1:length(res),function(i,res){
  re <- res[[i]]
  unique(sort(c(re$query.id,re$subject.id)))
},res)
names(scfs) <- names(res)

#' And how do they overlap?
#' As expected the lower percent identity contains all the others and the amount
#' of scaffolds decreases with increasing identity. Nevertheless the vast 
#' majority of scaffolds is present at a 100% identity.
plot.new()
grid.draw(venn.diagram(scfs[2:6],
                       filename=NULL,
                       col=pal[1:5],
                       category.names=names(scfs)[2:6])
)

#' Let us use the subset of 100% identity to identify redundant and artefactual
#' scaffolds; i.e. those having a 100% query cumulative coverage (redundant). If
#' the subject cumulative coverage is also a 100%, then there are considered 
#' artifacts.
sel <- res[["Perc.Ident.100"]]$query.cum.cov == 1
sprintf("There are %s scaffolds that are redundant",sum(sel))

#' Most of them are contained
plot(density(res[["Perc.Ident.100"]]$subject.cum.cov[sel]),
     main="subject coverage of the redundant scaffolds")

#' And a few are redundant
sel <- sel & res[["Perc.Ident.100"]]$subject.cum.cov == 1
sprintf("%s of which are artefactual",sum(sel))

#' Extend the annotation
annot <- read.delim("../../asp201/genome/annotation/Potra01-meta-matrix.tsv")
annot <- annot[,- grep("redund",colnames(annot))]
annot$redundant <- annot$ID %in% 
  res[["Perc.Ident.100"]]$query.id[res[["Perc.Ident.100"]]$query.cum.cov == 1]
annot$artefactual <- annot$ID %in% 
  res[["Perc.Ident.100"]]$query.id[res[["Perc.Ident.100"]]$query.cum.cov == 1
                                   & res[["Perc.Ident.100"]]$subject.cum.cov == 1]

#' Use a value of 97% identity and 97% coverage
#' Note that some scaffolds will be duplicated - i.e. have several hits.
#' We only report the best one here.
dat <- res[["Perc.Ident.97"]]
sel <- match(annot$ID,dat$query.id)
annot$haplotype.ID <- dat$subject.id[sel]
annot$haplotype.query.cum.cov <- dat$query.cum.cov[sel]
annot$haplotype.subject.cum.cov <- dat$subject.cum.cov[sel]
annot$putative.haplotype <- annot$haplotype.query.cum.cov >= 0.97
sel <- annot$putative.haplotype & ! is.na(annot$putative.haplotype) & ! annot$redundant

plot.multidensity(list(all=na.omit(annot$haplotype.subject.cum.cov),
                       hap=annot$haplotype.subject.cum.cov[sel]),
     main="subject coverage of the putative haplotype scaffolds",lwd=2)

plot.multidensity(list(all=na.omit(annot$haplotype.query.cum.cov),
                       hap=annot$haplotype.query.cum.cov[sel]),
     main="query coverage of the putative haplotype scaffolds",lwd=2)

sprintf("There are %s putative haplotype scaffolds",sum(sel))

#' # Save annotation
write.table(annot,"row.names"=FALSE,quote=FALSE,
            sep="\t",file="../../asp201/genome/annotation/Potra01-meta-matrix.tsv")

