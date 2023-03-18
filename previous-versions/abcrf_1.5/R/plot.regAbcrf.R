plot.regAbcrf <- function(x, n.var=min(30, length(x$model.rf$variable.importance)), ...){
  if (!inherits(x, "regAbcrf")) 
    stop("First argument not of class regAbcrf")
  variableImpPlot(x, n.var=n.var)
}
