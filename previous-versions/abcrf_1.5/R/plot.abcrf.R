plot.abcrf <- function(x, training, obs=NULL, n.var=20, pdf=FALSE, ...)
{
  
  if (!inherits(obs, "data.frame") && !is.null(obs) ) 
    stop("obs needs to be a data.frame object or NULL")
  if (!inherits(training, "data.frame"))
    stop("training needs to be a data.frame object")
  
	old.par <- par(no.readonly = TRUE)
	if (length(x$model.rf$variable.importance)<20) n.var <- length(x$model.rf$variable.importance)

	mf <- match.call(expand.dots=FALSE)
	mf <- mf[1]
	mf$formula <- x$formula
	mf$data <- training
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame() )
	mt <- attr(mf, "terms")
	modindex <- model.response(mf)
	
 	if (x$lda) {
 	  if (pdf) { 
 	    pdf("graph_varImpPlot.pdf")
		  variableImpPlot(x, n.var=n.var)
		  dev.off()
 	  }
 	  variableImpPlot(x, n.var=n.var)
		nmod <- length(x$model.rf$forest$levels)
		nstat <- x$model.rf$num.independent.variables
		projections <- predict(x$model.lda, training)$x
		if  (!is.null(obs)) projobs <- predict(x$model.lda, obs)$x
		coloris <- rainbow(nmod)
		colo <- coloris[modindex]
		readline("Press <ENTER> to Continue")
    if (nmod > 2) {
      if (pdf)
        {
        pdf("graph_lda.pdf")
        plot(projections[,1:2], col=colo, pch=3)
        legend("topleft", legend = as.character(x$model.rf$forest$levels), col = coloris, 
               pch = 15, bty = "o", pt.cex = 2, cex = .8, horiz = TRUE, 
               inset = c(.01, .01), title = "Models", bg = "white")
	  	  if  (!is.null(obs)) points(projobs[1],projobs[2],pch="*",cex=5.3)
	  	  dev.off()
		    }
      plot(projections[,1:2], col=colo, pch=3)
      legend("topleft", legend = as.character(x$model.rf$forest$levels), col = coloris, 
             pch = 15, bty = "o", pt.cex = 2, cex = .8, horiz = TRUE, 
             inset = c(.01, .01), title = "Models", bg = "white")
      if  (!is.null(obs)) points(projobs[1],projobs[2],pch="*",cex=5.3)
    } else {
      l1 <- x$model.rf$forest$levels[1]
      l2 <- x$model.rf$forest$levels[2]
      d1 <- density(projections[modindex == l1])
      d2 <- density(projections[modindex == l2])
      coloris <- c("blue", "orange")
      xrange <- range(c(d1$x, d2$x))
      yrange <- c(0, 1.2*max(c(d1$y, d2$y)))
      if (pdf)
        {
        pdf("graph_lda.pdf")
        plot(d1, xlim = xrange, ylim = yrange,
             col=coloris[1], main="", xlab="")
        lines(d2, col=coloris[2])
        legend("topleft", legend = as.character(x$model.rf$forest$levels), col = coloris, 
                cex = .8, horiz = TRUE, lty=1, bty="o",
               inset = c(.01, .01), title = "Models", bg = "white")
      	if  (!is.null(obs)) abline(v=projobs)
        dev.off()
      }
      plot(d1, xlim = xrange, ylim = yrange,
           col=coloris[1], main="", xlab="")
      lines(d2, col=coloris[2])
      legend("topleft", legend = as.character(x$model.rf$forest$levels), col = coloris, 
              cex = .8, horiz = TRUE, lty=1, bty="o",
             inset = c(.01, .01), title = "Models", bg = "white")
      if  (!is.null(obs)) abline(v=projobs)
    }
	} else {
	  if (pdf)
	    {
	    pdf("graph_varImpPlot.pdf")
	    variableImpPlot(x , n.var=n.var)
	    dev.off()
	    }
	  variableImpPlot(x , n.var=n.var)
	}
	par(old.par)
}