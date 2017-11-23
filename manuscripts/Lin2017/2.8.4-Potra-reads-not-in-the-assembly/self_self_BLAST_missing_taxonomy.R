#'---
#'author: Bastian Schiffthaler
#'title: Unaligned reads Potra to Potra. Taxonomy
#'output: html_document
#'---

knitr::opts_knit$set(root.dir = "/mnt/picea/projects/aspseq/tremula_tremuloides_comp/BLAST/")

library(data.table)
library(stringr)
library(RSQLite)
library(ggplot2)
setwd("/mnt/picea/projects/aspseq/tremula_tremuloides_comp/BLAST/")

#' Read output from self self BLAST
dat <- fread("asp201_unmapped_1.blt")

#' We care only/mostly about the best hit, therefore we select the first
#' BLAST hit (table is sorted) by sequence ID
dat.u <- unique(dat, by="V1")
#' Then, we extract the gi ID which corresponds to NCBI identifiers
gi <- as.numeric(str_extract(dat.u$V2,"\\d+"))
#' Connect to the local copy of the NCBI taxonomy SQLite database
con <- dbConnect(dbDriver("SQLite"),paste("/mnt/picea/storage/reference",
                                          "/Taxonomy/20150407/taxonomy.sqlite",
                                          sep = ""))
#' Formulate the SQL call for each gid
dbCall <- paste("SELECT gid,taxonomy.tid,nam,div_nam from gi_taxid",
                "left join taxonomy on gi_taxid.tid = taxonomy.tid",
                "left join node on taxonomy.tid = node.tax_id",
                "left join division on node.div_id = division.div_id",
                "WHERE class like 'scientific name%' AND gid IN (",
                paste(sort(gi),collapse=","),");")
#' Send the query and fetch all output
tax <- dbSendQuery(con, dbCall)
tax.dt <- data.table(fetch(tax, n = -1))
#' Plot the division as bars
ggplot(tax.dt, aes(x = div_nam)) + geom_bar()

#' Finally we want to select the genus, so we get the first "word" and
#' replace all the lower case statements
taxd <- str_extract(sapply(strsplit(tax.dt$nam, " "), "[[", 1), "^[A-Z].+")
taxd <- as.data.frame(table(taxd[!is.na(taxd)]))
#' Plot the top 50 hits
ggplot(head(taxd[order(taxd$Freq, decreasing = TRUE),],50),
       aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  xlab("Genus") + ylab("Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
