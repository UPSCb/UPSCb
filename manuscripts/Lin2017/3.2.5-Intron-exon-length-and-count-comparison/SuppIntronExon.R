#'---
#'title: Intron/Exon summary stats
#'author: Bastian Schiffthaler
#'output:
#'  html_document
#'---
#'
#' #Setting up and loading all the annotations:
suppressPackageStartupMessages({
library(GenomicFeatures);
library(ggplot2);
library(RColorBrewer);
library(genomeIntervals);
})
options(warn = -1)
suppressMessages(expr = {
Potri.gff <- makeTranscriptDbFromGFF("/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gff/Ptrichocarpa_v3.0_210_gene_exons.gff3");
Potra.gff <- makeTranscriptDbFromGFF("/mnt/picea/storage/reference/Populus-tremula/v1.0/gff3/Potra01-gene-complete.gff3");
Potrs.gff <- makeTranscriptDbFromGFF("/mnt/picea/storage/reference/Populus-tremuloides/v1.0/gff3/Potrs01-genome.gff3");
})
cols <- brewer.pal(10,"Paired")

#' # Counting the number of introns and exons per gene
#' 
#' `GenomicFeatures` does not define a "intronsBy()" type of function except for `intronsByTranscript`. I would like 
#' a function to count introns per gene, therefore I define a small funcction that first gets `elementLengths` from 
#' the `intronsByTranscript` function (therby counting the number of introns in each transcripts). I then get the
#' `elementLengths` for `transcriptsBy(x,"gene")` and repeat each element's position by its value, which yields an
#' index of equal length to the `elementLengths(intronsByTranscript())` vector that classifies each entry by the 
#' gene the transcript belongs to. Then I split the `elementLengths(intronsByTranscript())` by the aforementioned 
#' index and apply the statistical summarisation method (defaults to `max()`).
count.introns <- function(TxDB,count=max){
  iByT <- elementLengths(intronsByTranscript(TxDB))
  tByG <- elementLengths(transcriptsBy(TxDB,"gene"))
  tByG.ser <- unlist(lapply(seq_along(tByG), function(f){rep(f,unname(tByG[f]))}))
  iByG <- split(unname(iByT),tByG.ser)
  iByG.max <- lapply(iByG,count)
  return(unlist(iByG.max))
}

#' All the data is prepared as a `data.frame` so I can then easily `ggplot`.

Potri.intron.n <- data.frame("N" = count.introns(Potri.gff,max), 
                             "Species"="P.trichocarpa","Feature" = "intron")
Potra.intron.n <- data.frame("N" = count.introns(Potra.gff,max),
                             "Species"="P.tremula","Feature" = "intron")
Potrs.intron.n <- data.frame("N" = count.introns(Potrs.gff,max), 
                             "Species"="P.tremuloides","Feature" = "intron")
#' Doing the same transformation for exons is a lot simpler, as we only need to count `elementLengths` of
#' `exonsBy(x,"gene")`.
Potri.exon.n <- data.frame("N" = elementLengths(exonsBy(Potri.gff,"gene")),
                           "Species"="P.trichocarpa","Feature" = "exon")
Potra.exon.n <- data.frame("N" = elementLengths(exonsBy(Potra.gff,"gene")),
                            "Species"="P.tremula","Feature" = "exon")
Potrs.exon.n <- data.frame("N" = elementLengths(exonsBy(Potrs.gff,"gene")),
                           "Species"="P.tremuloides","Feature" = "exon")

#' To avoid typing all the objects into `rbind` I `lapply(function=get)` searching the envirnoment via `ls()`
#' and `grep()` and wrap the list into a `do.call(rbind,list)`.
df <- as.data.frame(do.call(rbind,
                            lapply(grep("intron\\.n|exon\\.n",
                                        ls(),value = TRUE),get)))

#+ fig.width=16, fig.height=9
ggplot(df, aes(y = N, x = Feature, fill = Species)) + geom_boxplot() + 
  scale_y_log10() + scale_fill_manual(values = cols[c(2,10,4,2,10,4)]) + 
  ylab("Number of exons or introns per gene")

#' The intron and exon medians are extremely close between all species, with _P. tremuloides_ lacking behind in genes
#' with >100 exons, while introns seem to be equivalent to the other two species. I will for now refrain from interpreting
#' this result, as _P. tremuloides_ still needs to have its annotation updated with Pasa and also lacks the inclusion
#' of RNA-Seq data, which both are very likely to result in a higher number of annotated exons.

```{r, empty}

#' # Getting the length of introns and exons per gene
#' 
#' ## Summarizing by the longest intron and exon
#' 
#' This function does work very similar to `count.introns()` which was defined before. The only big difference
#' is that instead of getting `elementLengths` of `intronsByTranscript` we can already run `width` and `max` on
#' the `intronsByTranscript` object, as it supports these operations. We still need to split later on to get the
#' longest intron of each gene. As `max` of an empty `integer()` returns `-Inf` we need to remove based on 
#' `is.infinite` twice.
length.introns <- function(TxDB,count=max){
  iByT <- max(width(intronsByTranscript(TxDB)))
  tByG <- elementLengths(transcriptsBy(TxDB,"gene"))
  tByG.ser <- unlist(lapply(seq_along(tByG), function(f){rep(f,unname(tByG[f]))}))
  iByG <- split(iByT, tByG.ser)
  iByG.max <- lapply(iByG,function(f){count(f[!is.infinite(f)])})
  return(unlist(iByG.max)[!is.infinite(unlist(iByG.max))])
}

Potri.intron.l <- data.frame("N" = length.introns(Potri.gff,max), 
                             "Species"="P.trichocarpa","Feature" = "intron")
Potra.intron.l <- data.frame("N" = length.introns(Potra.gff,max),
                             "Species"="P.tremula","Feature" = "intron")
Potrs.intron.l <- data.frame("N" = length.introns(Potrs.gff,max), 
                             "Species"="P.tremuloides","Feature" = "intron")

Potri.exon.l <- data.frame("N" = max(width(exonsBy(Potri.gff,"gene"))),
                           "Species"="P.trichocarpa","Feature" = "exon")
Potra.exon.l <- data.frame("N" = max(width(exonsBy(Potra.gff,"gene"))),
                           "Species"="P.tremula","Feature" = "exon")
Potrs.exon.l <- data.frame("N" = max(width(exonsBy(Potrs.gff,"gene"))),
                           "Species"="P.tremuloides","Feature" = "exon")

df <- as.data.frame(do.call(rbind,
                            lapply(grep("intron\\.l|exon\\.l",
                                        ls(),value = TRUE),get)))

#+ fig.width=16, fig.height=9
ggplot(df, aes(y = N, x = Feature, fill = Species)) + geom_boxplot() + 
  scale_y_log10() + scale_fill_manual(values = cols[c(2,10,4,2,10,4)]) + 
  ylab("Length of longest exon or intron per gene")

#' When considering the longest intron all three species again are quite similar. What is noteworthy in this case is 
#' that _P. trichocarpa_ seems to have many introns which are < 10bp, which seems highly suspicious. So I check the 
#' gff3

#' Here I start off by creating `GRanges` objects from the gff.
tmp <- readGff3("/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gff/Ptrichocarpa_v3.0_210_gene_exons.gff3")
tmp.gene <- GRanges(seq_name(tmp[tmp$type=="gene"]),
                    IRanges(start = tmp[tmp$type=="gene"]@.Data[,1], 
                            end=tmp[tmp$type=="gene"]@.Data[,2]))
tmp.exon <- GRanges(seq_name(tmp[tmp$type=="exon"]),
                    IRanges(start = tmp[tmp$type=="exon"]@.Data[,1], 
                            end=tmp[tmp$type=="exon"]@.Data[,2]))

diff.g.e <- setdiff(tmp.gene,tmp.exon)
length(diff.g.e[which(width(diff.g.e) <= 10 )])

#' This naïve approach seems to come to the same conclusion as TxDB, which is that there are 388 introns which are 
#' less than 10 bp in the _P. trichocarpa_ annotation. To dig deeper, I ran `gt gff3 -tidy --addintrons 
#' Ptrichocarpa_v3.0_210_gene_exons.gff3 > Ptrichocarpa_v3.0_210_gene_exons_introns.gff3` from the `genomeTools`
#' software suite. 

Potri.int <- readGff3("/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gff/Ptrichocarpa_v3.0_210_gene_exons_introns.gff3")
Potri.int.intron <- GRanges(seq_name(Potri.int[Potri.int$type=="intron"]),
                    IRanges(start = Potri.int[Potri.int$type=="intron"]@.Data[,1], 
                            end=Potri.int[Potri.int$type=="intron"]@.Data[,2]))
Potri.int.intron[width(Potri.int.intron) <= 10 ]

#' `genomeTools` comes to the same conclusion, but detects less (129) introns with length smaller or equal than 10
#' basepairs. This is due to the fact that `genomeTools` better accounts for gff3 hierarchy when inferring introns
#' (adding introns only between exons that are part of the same transcript), but `GenomicFeatures` accounts for the
#' same:

sum(unlist(width(intronsByTranscript(Potri.gff))) <= 10)

#' This at least assures that the previous results were not due an error in the determination of introns. Yet introns
#' this short do not reflect the biology accurately and we should have a look into the annotation.
```{r, empty}
#' ## Summarizing by the longest intron and exon
#' 
#' Instead of getting the longest introns, this time I am simply interested in overall intron and exon length.
#' Therefore, there is no need for factor splitting and we can use `GenomicFeatures` built-ins:

Potri.intron.lf <- data.frame("N" = unlist(width(intronsByTranscript(Potri.gff))), 
                             "Species"="P.trichocarpa","Feature" = "intron")
Potra.intron.lf <- data.frame("N" = unlist(width(intronsByTranscript(Potra.gff))),
                             "Species"="P.tremula","Feature" = "intron")
Potrs.intron.lf <- data.frame("N" = unlist(width(intronsByTranscript(Potrs.gff))), 
                             "Species"="P.tremuloides","Feature" = "intron")

Potri.exon.lf <- data.frame("N" = unlist(width(exonsBy(Potri.gff,"gene"))),
                           "Species"="P.trichocarpa","Feature" = "exon")
Potra.exon.lf <- data.frame("N" = unlist(width(exonsBy(Potra.gff,"gene"))),
                           "Species"="P.tremula","Feature" = "exon")
Potrs.exon.lf <- data.frame("N" = unlist(width(exonsBy(Potrs.gff,"gene"))),
                           "Species"="P.tremuloides","Feature" = "exon")

df <- as.data.frame(do.call(rbind,
                            lapply(grep("intron\\.lf|exon\\.lf",
                                        ls(),value = TRUE),get)))
#+ fig.width=16, fig.height=9
ggplot(df, aes(y = N, x = Feature, fill = Species)) + geom_boxplot() + 
  scale_y_log10() + scale_fill_manual(values = cols[c(2,10,4,2,10,4)]) + 
  ylab("Length of all exons or introns")
#' The trend is very similar to before, with _P. tremula_ having slightly longer introns as well as exons and
#' _P. tremuloides_ having slightly shorter ones. As mentioned before, this needs re-evaluating once the PASA
#' integration is complete for _P. tremuloides_. All in all, the three species seem to be very similar in terms
#' of exon and intron length.
#' 
sessionInfo()