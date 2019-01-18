#' ---
#' title: "Feature selection"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # What is this about?
#' 
#' This file is only a source file containing a function to select feature
#' provided they are expressed at an `exp` cutoff (exp) in at least `nrep` 
#' replicates of any `condition`.
#' 
#' To use it in your script do:
#' ```{r, eval=FALSE}
#' source("~/Git/UPSCb/src/R/featureSelection.R")
#' ```
#' 
#' In the following, you will first find the S4 generic of the function
#'  and then its implementation, as well as some description of its arguments.
#' 
#' # Generic
setGeneric(name="featureSelect",
           def=function(counts=matrix(),
                        conditions=factor(),
                        exp=1L,
                        nrep=2L
           ){
             standardGeneric("featureSelect")
           })

#' # Implementation
#' 
#' ## Arguments
#' The 'featureSelect' function expect the following arguments:
#' 
#' 1. `counts` the matrix of expression counts (normalised or not)
#' 
#' 2. `conditions` a factor that describes the different condition in the 
#' experiment
#'
#' It has 2 optional arguments
#' 
#' 1. `exp` the minimal value of expression to be used as a cutoff
#' 
#' 2. `nrep` the minimal number of replicates that need to have an expression
#' over the `exp` cutoff in any given condition
#'
#' ## Value
#' 
#' The function returns a boolean vector where TRUE indicates that the corresponding
#' feature passed the selected thresholds. The features are ordered as the 
#' rows of the input matrix

setMethod(f = "featureSelect", 
          signature = c("matrix","factor"),
          definition = function(
            counts=matrix(),
            conditions=factor(),
            exp=1L,
            nrep=2L
          ){
            
            ## validation
            stopifnot(all(dim(counts) > 1))
            stopifnot(nlevels(conditions)>=2)
            
            ## compute
            rowSums(sapply(lapply(
              split.data.frame(t(counts >= exp),conditions)
              ,colSums), ">=", nrep)) >= 1
                        
          })
