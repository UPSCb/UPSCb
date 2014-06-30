## this function is originally from the lSD package
## modified to suppress the change of margins

setGeneric(name="densityPlot",
           def=function(x, y, grid = 100, 
                        ncol = 30, 
                        nlevels = 10, ...){
  standardGeneric("densityPlot")
})

setMethod(f="densityPlot",
          signature=c("numeric","numeric"),
          definition=function(x, y, grid = 100, ncol = 30, 
                              nlevels = 10, ...){
            if (!is.vector(x) | !is.vector(y)) 
              stop("First two argument must be vectors !")
            if (length(x) != length(y)) 
              stop("Data vectors must be of the same length !")
            d = LSD:::kde2d.adj(x, y, n = grid)
            z <- d$z
            nrz <- nrow(z)
            ncz <- ncol(z)
            couleurs <- tail(topo.colors(trunc(1.4 * ncol)), ncol)
            image(d, col = couleurs, ...)
            contour(d, add = TRUE, nlevels = nlevels)
            box()
          })