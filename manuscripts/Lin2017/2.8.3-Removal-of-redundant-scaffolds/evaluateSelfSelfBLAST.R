#'---
#'title: Self self BLAST downstream processing Potra
#'author: Bastian Schiffthaler
#'---

setwd("/mnt/picea/projects/aspseq/tremula_tremuloides_comp/BLAST/")
source("~/Git/UPSCb/src/R/blastUtilities.R")
source("~/Git/private/misc/shitty_progress_bar.R")
source("~/Git/UPSCb/src/R/ggComparisonplot.R")
library(ggplot2)
library(Biostrings)
library(igraph)
###########
### Reading the BLAST output
###########

PotraBLAST <- readBlast(file = "Potra_self.txt", ignoreSelf = TRUE, 
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
                        plot = FALSE)
              
###########
### Quick plots to check cumulative coverage
###########
head(PotraBLAST$df)
s <- sample(1:nrow(PotraBLAST$df),2e5,FALSE)
ggplot(PotraBLAST$df[s,],aes(x=query.cum.cov,y=subject.cum.cov)) + 
  stat_density2d(geom = "tile", contour = FALSE, aes(fill=..density..)) + 
  scale_fill_gradient2()

comparisonplot(PotraBLAST$df$query.cum.cov[s], PotraBLAST$df$subject.cum.cov[s])

ggComparisonPlot(PotraBLAST$df$query.cum.cov[s], PotraBLAST$df$subject.cum.cov[s])

#In order to get a scaffold for each query with associated stats, we first
#split into a list
PotraBLAST.blf.split <- split(PotraBLAST$blf,PotraBLAST$blf$query.id)
#We deselect alignments shorter than 100bp, as it often involves garbage scaffolds
#which should be removed in the final assembly
PotraBLAST.blf.split <-mclapply(PotraBLAST.blf.split, mc.cores = 8, function(f){
  f[f$alignment.length>=100,]
})
#Next, we sort each list item descendingly by percent identity
PotraBLAST.blf.split <-mclapply(PotraBLAST.blf.split, mc.cores = 8, function(f){
  f[order(f$percent.identity, decreasing = TRUE),]
})
#We keep the top hit [TODO test if `head` is faster then `[` ]
PotraBLAST.blf.split.top <- do.call(rbind,lapply(PotraBLAST.blf.split,function(f){f[1,]}))
#Create a match vector to match against rownames in the $df
query.match.vec <- paste(PotraBLAST.blf.split.top[,1], PotraBLAST.blf.split.top[,2], sep="+")
#Match the vector and rownames
query.match <- match(query.match.vec, row.names(PotraBLAST$df))
#Result is used as selector
PotraBLAST.df.top <- PotraBLAST$df[query.match,]

#Fail if $df and $blf are not the same length
if(!nrow(PotraBLAST.blf.split.top)==nrow(PotraBLAST.df.top)){
  stop("Error when finding cumulative overlap for some scaffolds")
}

meta.df <- data.frame(cbind(PotraBLAST.blf.split.top[,1:4],PotraBLAST.df.top[,3:4]))

#fix class errors
meta.df <- as.data.frame(apply(meta.df,2,as.character), stringsAsFactors = FALSE)
meta.df[,3] <- as.numeric(meta.df[,3])
meta.df[,4] <- as.integer(meta.df[,4])
meta.df[,5] <- as.numeric(meta.df[,5])
meta.df[,6] <- as.numeric(meta.df[,6])

###########
### Average percent identity of scaffolds = 100% cum. cov. of subject and query
###########
Sca.index <- PotraBLAST$df[,3]==1 & PotraBLAST$df[,4]==1

# How many unique subjects are >95% cum.cov.
length(unique(PotraBLAST$df$subject.id[Sca.index]))

qID <- PotraBLAST$df[Sca.index,1]
sID <- PotraBLAST$df[Sca.index,2]

# Split BLAST results by query ID
byQuery <- split(PotraBLAST$blf, factor(PotraBLAST$blf$query.id))

# Define function to calculate mean perc. id. given an index
getAvgPId <- function(ind){
  query.id <- qID[ind]
  subject.id <- sID[ind]
  sel <- byQuery[[query.id]]
  return(c(query.id, subject.id,
    as.numeric(mean(sel[subject.id %in% sel$subject.id,3]))
  ))
} 

# Genome to naormalize scaffold factor indices
Genome <- readDNAStringSet("/mnt/picea/storage/reference/Populus-tremula/v1.0/fasta/Potra01-genome.fa")

tmp <- mclapply(seq_along(qID),mc.cores = 8, function(f){
  getAvgPId(f)
})

tmp <- as.data.frame(do.call(rbind,tmp), stringsAsFactors = FALSE)

tmp[,3] <- as.numeric(tmp[,3])
tmp[,1] <- factor(tmp[,1],levels = names(Genome))
tmp[,2] <- factor(tmp[,2],levels = names(Genome))

colnames(tmp) <- c("query","subject","avg.identity")

# Also define query and subject scaffold length
tmp$query.width <- width(Genome[as.character(tmp[,1])])
tmp$subject.width <- width(Genome[as.character(tmp[,2])])

#############
### Plots of Avg. Perc. Identity
#############
ggplot(tmp, aes(x = query.width, y = subject.width)) + geom_hex() + 
  ggtitle("Scaffold width relationship(s) where cumulative coverage of query and subject ==1")

ggplot(tmp, aes(x = log10(query.width), y = log10(subject.width))) + 
  geom_point(aes(color = avg.identity)) + geom_density2d(color = "firebrick2", size = 1.5) + 
  ggtitle("Scaffold width relationship(s) where cum.cov.=1, points shaded by avg. % ident")

ggplot(tmp[tmp$avg.identity>99,], aes(x = query.width, y = subject.width)) + 
  geom_point(aes(color = avg.identity)) + geom_smooth(method = "lm") + 
  ggtitle("Scaffold width relationship(s) where cum.cov.=1, points shaded by avg. % ident")

annotFile <- "../../asp201//genome/annotation/Potra01-meta-matrix.tsv"
annot <- read.delim(annotFile)
tmp2 <- tmp[tmp$avg.identity>99,]
tmp2$query.cov.mean <- annot[annot$ID %in% as.character(tmp2[,1]),"cov.mean"]
tmp2$subject.cov.mean <- annot[annot$ID %in% as.character(tmp2[,2]),"cov.mean"]

dup <- duplicated(apply(tmp2[,1:2],1,function(f){
  paste(sort(f),collapse = "_")
}))

scaSel <- apply(tmp2[,4:7],1,function(f){
  if(f[1]!=f[2]){
    return(which.min(f[1:2]))
  } else {
    return(which.max(f[3:4]))
  }
})

res <- unique(sapply(1:nrow(tmp2),function(f){
  if(!dup[f]){
    as.character(tmp2[f,scaSel[f]])} else {NA}
}))

gr <- graph.edgelist(as.matrix(tmp2[,1:2]),directed = FALSE)
gr <- set.vertex.attribute(gr,"size", value = width(Genome[get.vertex.attribute(gr,"name")])/max(width(Genome[get.vertex.attribute(gr,"name")]))*15)
cols <- sapply(get.vertex.attribute(gr,"name") %in% res, function(tes){
  if(tes){
    return("firebrick2")
  } else {
    return("skyblue")
  }
})
gr <- set.vertex.attribute(gr,"color",value = cols)
gr <- set.vertex.attribute(gr,"width",value = width(Genome[get.vertex.attribute(gr,"name")]))
plot.igraph(gr,vertex.size = get.vertex.attribute(gr,"size"), 
            vertex.label = get.vertex.attribute(gr,"width"),
            vertex.color = get.vertex.attribute(gr,"color"),
            layout = layout.grid)

res <- res[!is.na(res)]
res
resMatch <- sapply(seq_along(scaSel), function(f){
  se <- ifelse(scaSel[f]==1,2,1)
  sca <- as.character(tmp2[f,se])
  sca
})[!dup]
annot$putative.redundant <- NA
annot$putative.redundant[annot$ID%in%res] <- resMatch

write.table(annot, file=annotFile, quote = FALSE, col.names = TRUE, row.names = FALSE)
