abcrf.formula <- function(formula, data, lda=TRUE, ntree=500, sampsize=min(1e5, nrow(data)), paral=FALSE, 
                          ncores= if(paral) max(detectCores()-1,1) else 1, ...) 
{
  
  # formula and data.frame check
  
  if (!inherits(formula, "formula"))
    stop("abcrf.formula is only for formula objects")
  if (!inherits(data, "data.frame"))
    stop("data needs to be a data.frame object")
  if(ncores > detectCores() || ncores < 1)
    stop("incorrect number of CPU cores")
  if ( (!is.logical(paral)) && (length(paral) != 1L) )
    stop("paral should be TRUE or FALSE")
  
  # modindex and sumsta recovery
  
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data"), names(mf))
  mf <- mf[c(1L,m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame() )
  
  if (!is.factor(model.response(mf)))
    stop("response should be a factor containing the model indexes")
  if (nrow(data) == 0L || is.null(nrow(data)))
    stop("no simulation in the reference table (response, sumstat)")
  if ( (!is.logical(lda)) && (length(lda) != 1L) )
    stop("lda should be TRUE or FALSE")
  
  if (lda) {
    model.lda <- lda(formula, data)
    data <- cbind(data, as.matrix(predict(model.lda, data)$x ))
    } else { 
    model.lda <- NULL
  }
  
  m <- names(match.call(expand.dots=TRUE))
  if ((!"sampsize" %in% m) && (nrow(data) <= 15)) 
    sampsize <- as.integer(sampsize / 10)
  if(sampsize > nrow(data))
    stop("sampsize too large")
  
  model.rf <- ranger(formula, data, num.trees=ntree, sample.fraction=sampsize/nrow(data), 
                     num.threads = ncores, keep.inbag = TRUE, importance = 'impurity', ...)
  
  # class error computation
  class.error = vector()
  for (i in 1:nrow(model.rf$confusion.matrix)) {
    rowSum <- sum(model.rf$confusion.matrix[i,])
    accurate <- diag(model.rf$confusion.matrix)[i]
    error <- rowSum - accurate
    class.error[i] <- error / rowSum
  }
  model.rf$confusion.matrix <- cbind(model.rf$confusion.matrix, class.error)
  colnames(model.rf$confusion.matrix) <- c(paste(model.rf$forest$levels),"class.error")
  
  model.rf$model.rf
  
  cl <- match.call()
  cl[[1]] <- as.name("abcrf")

  x <- list(call=cl, lda=lda, formula=formula, model.rf=model.rf, model.lda=model.lda, prior.err=model.rf$prediction.error)
  class(x) <- "abcrf"
  x
}

abcrf.default <- function(...) {
  cl <- match.call()
  cl[[1]] <- as.name("abcrf")
  cat("call:\n")
  print(cl)
  stop("the first argument should be a formula")
}
  
abcrf <- function(...) UseMethod("abcrf")

print.abcrf <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  if (x$lda) cat("includes the axes of a preliminary LDA\n\n")  
  cat("Number of simulations: ", length(x$model.rf$num.samples), "\n", sep="")
  cat("     Prior error rate: ", round(x$prior.err * 100, digits = 4), "%\n\n", sep = "")
  cat("Confusion matrix:\n")
  print(x$model.rf$confusion, ...)  
}