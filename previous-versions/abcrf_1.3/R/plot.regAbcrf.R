plot.regAbcrf <-
function(x, n.var=min(30, nrow(x$model.rf$importance)), ... ){
  if (!inherits(x, "regAbcrf")) 
    stop("First argument not of class regAbcrf")
  varImpPlot(x$model.rf, n.var=n.var, ...)
}
