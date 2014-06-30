"plot.multidensity" <- function(x,xlab="x",col=brewer.pal(8,"Dark2"),
                                legend.x="top",xlim=NULL,ylim=NULL,
                                lty=1,legend.cex=1,...){

  ## pckg
  library(RColorBrewer)
  
  ## check
  stopifnot(is.list(x))
  
  ## densities
  dens <- lapply(x,density)
  if(is.null(xlim)){
    xlim <- range(sapply(dens,"[[","x"))
  }
  if(is.null(ylim)){
    ylim <- range(sapply(dens,"[[","y"))
  }

  ## lty
  if(length(lty)==1){
    lty <- rep(lty,length(dens))
  }
  
  ## plot
  plot(0,0,xlim=xlim,ylim=ylim,ylab="density",type="n",xlab=xlab,...)
  
  ## lines
  sapply(1:length(dens),function(i,dens,col,...){lines(dens[[i]],col=col[i],lty=lty[i],...)},dens,col,...)
  
  ## legend
  legend(legend.x,col=col[1:length(x)],bty="n",legend=names(x),lty=lty,cex=legend.cex)
}

