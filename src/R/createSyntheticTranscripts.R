"createSyntheticTranscripts" <- function(filename,
                                         input = c("gff3","gtf"),
                                         features = c("mRNA", "miRNA", "tRNA", "transcript"),
                                         output = c("Genome_intervals","GRanges")) {

  #' load libraries
  stopifnot(require(genomeIntervals))
  stopifnot(require(S4Vectors))
  stopifnot(require(IRanges))
  stopifnot(require(easyRNASeq))

  #' first check
  stopifnot(file.exists(filename))

  #' get the values
  input <- match.arg(input)
  features <- match.arg(features, several.ok = TRUE)
  output <- match.arg(output)

  #' define some global variables
  relation <- switch(input,
                     "gff3"=list(ID="ID",Parent="Parent"),
                     "gtf"=list(ID="transcript_id",Parent="gene_id"))

  #' read the gff3/gtf file
  dat <- readGff3(filename)

  #' If gtf, reformat the attributes and drop the double quotes
  if(input=="gtf"){
    dat$gffAttributes <- gsub("\"","",easyRNASeq:::.convertGffToGtfAttributes(dat$gffAttributes))
  }

  ## get the gene <-> mRNA/transcript map
  # This is mRNA IDs and their parents (genes)
  sel <- dat$type %in% features

  # That step would not necessary for gtf, but it is easier to implement in a
  # similar way for both format
  idMap <- data.frame(type = dat[sel]$type,
                      getGffAttribute(dat[sel],relation$ID),
                      getGffAttribute(dat[sel],relation$Parent))

  ## extract the exons and group by gene ID
  sel <- dat$type == "exon"

  ## we can drop multiple Parents (i.e. comma separated Parent values as we are
  ## collapsing them anyway)
  mRnaID <- sub(",.*","",getGffAttribute(dat[sel],switch(input,
                                                         "gff3"=relation$Parent,
                                                         "gtf"=relation$ID)))

  ## avoid unwanted features
  rngs <- IRanges(start = dat[sel, 1],
                  end = dat[sel, 2])[mRnaID %in% idMap[,relation$ID]]

  ## create a set of synthetic exons
  rngList <- IRanges::reduce(
    IRanges::split(rngs,
                   idMap[match(mRnaID[mRnaID %in% idMap[,relation$ID]],
                               idMap[,relation$ID]),relation$Parent]))

  ## export the gene, exon and features as gff3
  ## create the new gff object
  ## select the gene
  sel <- dat$type == "gene"

  ## create the gene gff
  geneID <- getGffAttribute(dat[sel],switch(input,
                                            "gff3"=relation$ID,
                                            "gtf"=relation$Parent))
  geneGff <- dat[sel][geneID %in% idMap[,relation$Parent]]
  if (input == "gtf"){
    geneGff$gffAttributes <- sub(relation$Parent,"ID",geneGff$gffAttributes)
  }

  ## create gffs for each feature
  featureGff <- Reduce(c, lapply(features, function(f) {
    f.sel <- geneID %in% idMap[,relation$Parent][idMap$type == f]
    fGff <- dat[sel][f.sel]
    fGff$type <- f
    fGff$gffAttributes <- paste("ID=",
                                getGffAttribute(fGff,relation$Parent),
                                ".0;Parent=",
                                getGffAttribute(fGff,relation$Parent),
                                sep="")
    fGff
  }))

  ## create the exon gff
  rngList <- rngList[match(geneID[geneID %in% idMap[,relation$Parent]], names(rngList))]
  exonNumber <- elementLengths(rngList)
  exonGff <- dat[rep(which(sel)[geneID %in% idMap[,relation$Parent]], exonNumber)]
  exonGff[,1] <- IRanges::unlist(start(rngList))
  exonGff[,2] <- IRanges::unlist(end(rngList))

  exonID <- sapply(exonNumber, ":", 1)
  sel <- geneGff$strand == "+"
  exonID[sel] <- sapply(exonID[sel], rev)
  ID <- getGffAttribute(exonGff, switch(input,"gff3"=relation$ID,"gtf"=relation$Parent))
  exonGff$gffAttributes <- paste0("ID=", paste(ID, "exon", unlist(exonID, use.names=FALSE), sep="."),
                                 ";Name=", paste(ID, "exon", unlist(exonID, use.names=FALSE), sep="."),
                                 ";Parent=", paste(ID,"0",sep = "."))
  exonGff$type <- "exon"

  ## combine
  newgff <- c(geneGff, featureGff, exonGff)

  ## change the source
  newgff$source <- "UPSC"

  ## sort
  newgff <- newgff[order(seq_name(newgff), newgff[, 1],
                         factor(as.character(newgff$type),
                                labels = seq_len(2 + length(features)),
                                levels = c("gene", features, "exon"))), ]

  return(switch(output,
                "Genome_intervals" = newgff,
                "GRanges" = as(newgff, "GRanges"),
                stop(paste("Cannot generate output of type", output))))
}
