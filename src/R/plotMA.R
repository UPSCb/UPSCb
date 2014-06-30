setMethod(f="plotMA",
          signature="DataFrame",
          definition=function(object,alpha=0.001){
            
            ## lib
            require(LSD)
            
            ## check
            if(!existsFunction("densityPlot")){
              stop("Load the densityPlot function prior to using this function.")
            }
            
            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha 
            
            ## graphic params
            orig.par <- par(no.readonly=TRUE)
            par(mfrow=c(2,1))
            
            ## plots
            densityPlot(log10(object$baseMean[sel]),
                      object$log2FoldChange[sel],
                      grid=250,ncol=30,nlevels=10,
                      main="MA density estimation"
                      )
            mtext("log10 mean expression",side=2,line=2)
            mtext("log2 FC",side=1,line=2)
            
            heatscatter(log10(object$baseMean[sel]),
                        object$log2FoldChange[sel],
                        add.contour=TRUE,main="MA",
                        xlab="log10 mean expression",
                        ylab="log2 FC",sub=paste(sum(sel2),
                                                 "sign. feats. @",
                                                 alpha,"cutoff"))
            
            points(log10(object$baseMean[sel][sel2]),
                   object$log2FoldChange[sel][sel2],
                   col="darkred",pch=19,cex=.5)
            
            legend("topright",pch=19,col="darkred","sign. feats.")
            
            par(orig.par,no.readonly=TRUE)
            invisible(TRUE)
          })
