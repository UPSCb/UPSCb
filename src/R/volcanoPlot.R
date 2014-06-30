setGeneric(name="volcanoPlot",def=function(object,alpha=0.001){
  standardGeneric("volcanoPlot")
})

setMethod(f="volcanoPlot",
          signature="DataFrame",
          definition=function(object,alpha=0.001){

            ## lib
            require(LSD)
            
            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha 

            ## plot
            heatscatter(object$log2FoldChange[sel],
                        -log10(object$padj[sel]),
                        main="Volcano",xlab="Log2 Fold Change", 
                        ylab="- log(10) adj. p-value")

            ## legend
            legend("topleft",bty="n",paste("cutoff @",alpha),lty=2,col="gray")
            
            ## points
            points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col="lightblue",pch=19)
            points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col="dodgerblue3",pch=19,cex=0.5)

            ## circle the points for the dot plot
            abline(h=-log10(alpha),lty=2,col="gray")
            })
