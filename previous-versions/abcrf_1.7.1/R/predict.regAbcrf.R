predict.regAbcrf <- function(object, obs, training, quantiles=c(0.025,0.975),
                             paral = FALSE, ncores = if(paral) max(detectCores()-1,1) else 1, rf.weights = FALSE,...)
{
  ### Checking arguments

  if (!inherits(obs, "data.frame"))
    stop("obs needs to be a data.frame object")
  if (!inherits(training, "data.frame"))
    stop("training needs to be a data.frame object")
  if (nrow(training) == 0L || is.null(nrow(training)))
    stop("no simulation in the training reference table (response, sumstat)")
  if ( (!is.logical(paral)) || (length(paral) != 1L) )
    stop("paral should be TRUE or FALSE")
  if ( (!is.logical(rf.weights)) || (length(rf.weights) != 1L) )
    stop("paral should be TRUE or FALSE")
  if(is.na(ncores)){
    warning("Unable to automatically detect the number of CPU cores, \n1 CPU core will be used or please specify ncores.")
    ncores <- 1
  }
  
  if(min(quantiles)<0 | max(quantiles)>1 )
    stop("quantiles must be in [0,1]")
  
  # modindex and sumsta recovery
  
  mf <- match.call(expand.dots=FALSE)
  mf <- mf[1]
  mf$formula <- object$formula
  

  mf$data <- training
  
  training <- mf$data
  
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame() )
  mt <- attr(mf, "terms")
  
  obj <- object$model.rf
  inbag <- matrix(unlist(obj$inbag.counts, use.names=FALSE), ncol=obj$num.trees, byrow=FALSE)
  
  obj[["origNodes"]] <- predict(object$model.rf, training, predict.all=TRUE, num.threads=ncores)$predictions
  obj[["origObs"]] <- model.response(mf)
  
  x <- obs
  if(!is.null(x)){
    if(is.vector(x)){
      x <- matrix(x,ncol=1)
    }
    if (nrow(x) == 0)
      stop("obs has 0 rows")
    if (any(is.na(x)))
      stop("missing values in obs")
  }

  
  ### prediction
  
  origObs <- obj$origObs
  origNodes <- obj$origNodes
  
  nnew <- nrow(x)
  
  quant <- matrix(nrow=nnew,ncol=length(quantiles))
  mediane <- matrix(nrow=nnew, ncol=1)
  
  nodes <- predict(object$model.rf, x, predict.all=TRUE, num.threads=ncores)$predictions
  if(is.null(dim(nodes))) nodes <- matrix(nodes, nrow=1)
  ntree <- obj$num.trees
  
  nobs <- object$model.rf$num.samples

  weights <- findweights(as.matrix(origNodes), as.matrix(inbag), as.matrix(nodes), as.integer(nobs), as.integer(nnew), as.integer(ntree)) # cpp function call

  weights.std <- weights/ntree
  
  esper <- sapply(1:nnew, function(x) weights.std[,x]%*%origObs)
  
  # Out of bag expectations
  
  predict.oob <- object$model.rf$predictions
  
  # squared residuals
  
  residus.oob.sq <- (origObs - predict.oob)^2
  
  # variance estimation
  
  variance <- sapply(1:nnew, function(x) weights.std[,x] %*% residus.oob.sq)
  
  ## Variance obtained using cdf
  
  variance.cdf <- sapply(1:nnew, function(x) weights.std[,x] %*% ( ( origObs - esper[x] )^2 ) ) 
  
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
  if(rf.weights == TRUE){
    tmp <- list(expectation = esper, med = mediane, variance = variance, variance.cdf = variance.cdf, quantiles = quant, weights=weights.std)
  } else{
    tmp <- list(expectation = esper, med = mediane, variance = variance, variance.cdf = variance.cdf, quantiles = quant)
  }
  class(tmp) <- "regAbcrfpredict"
  tmp
}


print.regAbcrfpredict <-
  function(x, ...){
    ret <- cbind(x$expectation, x$med, x$variance, x$variance.cdf, x$quantiles)
    colnames(ret) <- c("expectation", "median", "variance", "variance.cdf" , colnames(x$quantiles) )
    print(ret, ...)
  }

as.data.frame.regAbcrfpredict <-
  function(x, ...) {
    ret <- cbind(x$expectation, x$med, x$variance, x$variance.cdf, x$quantiles)
    colnames(ret) <- c("expectation", "median", "variance", "variance.cdf" , colnames(x$quantiles) )
    as.data.frame(ret,  row.names=NULL, optional=FALSE, ...)
  }

as.matrix.regAbcrfpredict <-
  function(x, ...){
    ret <- cbind(x$expectation, x$med, x$variance, x$variance.cdf, x$quantiles)
    colnames(ret) <- c("expectation", "median", "variance", "variance.cdf" , colnames(x$quantiles) )
    ret
  }

as.list.regAbcrfpredict <-
  function(x, ...){
    if(is.null(x$weights)){
      list(expectation = x$expectation, med = x$med , variance = x$variance, x$variance.cdf, quantiles=x$quantiles, ...)
    } else{
      list(expectation = x$expectation, med = x$med , variance = x$variance, x$variance.cdf, quantiles=x$quantiles, weights=x$weights, ...)
    }
  }