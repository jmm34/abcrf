predictOOB.regAbcrf <- function(object, training, quantiles=c(0.025,0.975),
                             paral = FALSE, ncores = if(paral) max(detectCores()-1,1) else 1, rf.weights = FALSE,...)
{
  ### Checking arguments
  
  if (!inherits(training, "data.frame"))
    stop("training needs to be a data.frame object")
  if (nrow(training) == 0L || is.null(nrow(training)))
    stop("no simulation in the training reference table (response, sumstat)")
  if ( (!is.logical(paral)) && (length(paral) != 1L) )
    stop("paral should be TRUE or FALSE")
  if ( (!is.logical(rf.weights)) && (length(rf.weights) != 1L) )
    stop("paral should be TRUE or FALSE")
  if( ncores > detectCores() || ncores < 1 )
    stop("incorrect number of CPU cores")
  if(min(quantiles)<0 | max(quantiles)>1 )
    stop("quantiles must be in [0,1]")
  
  # modindex and sumsta recovery
  
  mf <- match.call(expand.dots=FALSE)
  mf <- mf[1]
  mf$formula <- object$formula
  mf$data <- training
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame() )
  mt <- attr(mf, "terms")
  
  obj <- object$model.rf
  inbag <- matrix(unlist(obj$inbag.counts, use.names=FALSE), ncol=obj$num.trees, byrow=FALSE)
  
  obj[["origNodes"]] <- predict(object$model.rf, training, predict.all=TRUE, num.threads=ncores)$predictions
  obj[["origObs"]] <- model.response(mf)
  
  ### prediction
  
  origObs <- obj$origObs
  origNodes <- obj$origNodes
  
  nobs <- object$model.rf$num.samples
  
  quant <- matrix(nrow=nobs,ncol=length(quantiles))
  mediane <- matrix(nrow=nobs, ncol=1)
  
  ntree <- obj$num.trees
  
  result <- predictOob_cpp(origNodes = as.matrix(origNodes), inbag = as.matrix(inbag), nobs = as.integer(nobs), ntree = as.integer(ntree))
  
  weights <- matrix(result, nrow= nobs)
  weights.std <- weights
  
  esper <- sapply(1:nobs, function(x) weights.std[,x]%*%origObs)
  
  # Out of bag expectations
  
  predict.oob <- object$model.rf$predictions
  
  # squared residuals
  
  residus.oob.sq <- (origObs - predict.oob)^2
  
  # variance estimation
  
  variance <- sapply(1:nobs, function(x) weights.std[,x] %*% residus.oob.sq)
  
  ## Variance obtained using cdf
  
  variance.cdf <- sapply(1:nobs, function(x) weights.std[,x] %*% ( ( origObs - esper[x] )^2 ) ) 
  
  # Quantiles calculation
  
  ord <- order(origObs)
  origObs <- origObs[ord]
  weights <- weights[ord,,drop=FALSE]
  cumweights <- colCumsums(weights)
  cumweights <- sweep(cumweights,2,as.numeric(cumweights[nobs,]),FUN="/")
  
  # quantiles (from Meins)
  
  for (qc in 1:length(quantiles)){
    larg <- cumweights<quantiles[qc]
    wc <- colSums(larg)+1
    ind1 <- which(wc<1.1)
    indn1 <- which(wc>1.1)
    quant[ind1,qc] <- rep(origObs[1],length(ind1))
    quantmax <- origObs[wc[indn1]]
    quantmin <- origObs[wc[indn1]-1]
    weightmax <- cumweights[cbind(wc[indn1],indn1)]
    weightmin <- cumweights[cbind(wc[indn1]-1,indn1)]
    factor <- numeric(length(indn1))
    indz <- weightmax-weightmin<10^(-10)
    factor[indz] <- 0.5
    factor[!indz] <- (quantiles[qc]-weightmin[!indz])/(weightmax[!indz]-weightmin[!indz])
    quant[indn1,qc] <- quantmin + factor* (quantmax-quantmin)
  }
  
  colnames(quant) <- paste("quantile=",quantiles)
  
  # mediane estimation
  
  larg <- cumweights< 0.5
  wc <- colSums(larg)+1
  ind1 <- which(wc<1.1)
  indn1 <- which(wc>1.1)
  mediane[ind1,1] <- rep(origObs[1],length(ind1))
  quantmax <- origObs[wc[indn1]]
  quantmin <- origObs[wc[indn1]-1]
  weightmax <- cumweights[cbind(wc[indn1],indn1)]
  weightmin <- cumweights[cbind(wc[indn1]-1,indn1)]
  factor <- numeric(length(indn1))
  indz <- weightmax-weightmin<10^(-10)
  factor[indz] <- 0.5
  factor[!indz] <- (0.5-weightmin[!indz])/(weightmax[!indz]-weightmin[!indz])
  mediane[indn1,1] <- quantmin + factor* (quantmax-quantmin)
  
  MSE = mean( (obj$origObs - esper)^2)
  NMAE = mean( abs((obj$origObs - esper)/obj$origObs) )
  
  coverage = NULL
  mean.q.range = NULL
  if(length(quantiles)==2){
    if(quantiles[1] < quantiles[2]){
      coverage = mean( (quant[,1] <= obj$origObs) && (obj$origObs <= quant[,2]) )
    } else{
      coverage = mean( (quant[,2] <= obj$origObs) && (obj$origObs <= quant[,1]) )
    }
    mean.q.range = mean(abs((quant[,2] - quant[,1])/obj$origObs))
  }

  if(rf.weights == TRUE){
    tmp <- list(expectation = esper, med = mediane, variance = variance, variance.cdf = variance.cdf, quantiles = quant, weights=weights.std, MSE = MSE, NMAE = NMAE, coverage = coverage, mean.q.range = mean.q.range)
  } else{
    tmp <- list(expectation = esper, med = mediane, variance = variance, variance.cdf = variance.cdf, quantiles = quant, MSE = MSE, NMAE = NMAE, coverage = coverage, mean.q.range = mean.q.range)
  }
  class(tmp) <- "regAbcrfOOBpredict"
  tmp
}

predictOOB <-
  function(...) UseMethod("predictOOB")

print.regAbcrfOOBpredict <-
  function(x, ...){
    cat("\nOut-of-bag mean squared error: ", x$MSE, "\n")
    cat("Out-of-bag normalized mean absolute error: ", x$NMAE, "\n")
    if(!is.null(x$coverage)) cat("Out-of-bag credible interval coverage: ", x$coverage, "\n")
    if(!is.null(x$mean.q.range)) cat("Out-of-bag credible interval relative range: ", x$mean.q.range)
  }

as.list.regAbcrfOOBpredict <-
  function(x, ...){
    if(is.null(x$weights)){
      list(expectation = x$expectation, med = x$med , variance = x$variance, x$variance.cdf, quantiles=x$quantiles,
           MSE = x$MSE, NMAE = x$NMAE, coverage = x$coverage, mean.q.range = x$mean.q.range, ...)
    }else{
      list(expectation = x$expectation, med = x$med , variance = x$variance, x$variance.cdf, quantiles=x$quantiles,
           weights = x$weights, MSE = x$MSE, NMAE = x$NMAE, coverage = x$coverage, mean.q.range = x$mean.q.range, ...) 
    }
  }