variableImpPlot <- function(object, n.var=min(30, length(object$model.rf$variable.importance)))
{
  if (!inherits(object, "regAbcrf") && !inherits(object, "abcrf") )
    stop("object not of class abcrf or regAbcrf")
  imp <- object$model.rf$variable.importance
  ord <- rev(order(imp, decreasing = TRUE)[1:n.var])
  xmin <- 0
  dotchart(imp[ord], xlab = 'Impurity', ylab = "", main = "Variable Importance", xlim = c(xmin, max(imp)))
  invisible(imp)
}