predict.regAbcrf <-
function(object, newdata, quantiles=c(0.025,0.975) , ...){

  train.data <- object$sumsta
  obj <- object$model.rf
  inbag <- obj$inbag
  x <- newdata

  class(obj) <- c("quantregForest", "randomForest")
  obj[["origNodes"]] <- predict(object$model.rf, train.data, predict.all=TRUE)$individual
  obj[["origObs"]] <- object$model.rf$y
  obj[["importance"]] <- object$model.rf$importance[,-1]
  obj[["quantiles"]] <- NULL

  ### Checking arguments
  if (!inherits(object, "regAbcrf"))
    stop("object not of class regAbcrf")
  if(min(quantiles)<0 | max(quantiles)>1 )
    stop("quantiles must be in [0,1]")

  x <- newdata
  if(!is.null(x)){
    if(is.vector(x)){
      x <- matrix(x,ncol=1)
    }
    if (nrow(x) == 0)
      stop("newdata has 0 rows")
    if (any(is.na(x)))
      stop("missing values in newdata")
  }

  ### prediction

  origObs <- obj$origObs
  origNodes <- obj$origNodes

  quant <- matrix(nrow=nrow(x),ncol=length(quantiles))
  mediane <- matrix(nrow=nrow(x), ncol=1)
  nodes <- predict(object$model.rf, x, predict.all=TRUE)$individual
  ntree <- obj$ntree

  nobs <- length(origObs)
  nnew <- nrow(x)
  normalise <- 1
  weightvec <- rep(0,nobs*nnew)
  counti <- rep(0,nobs)
  thres <- 5*.Machine$double.eps
  result <- .C("findweightsAmelioree",
               as.double(as.vector(origNodes)),
               as.integer(as.vector(inbag)),
               as.double(as.vector(nodes)),
               weightvec=as.double(weightvec),
               as.integer(nobs),
               as.integer(nnew),
               as.integer(ntree),
               as.double(thres),
               as.integer(counti),
               as.integer(normalise),
               PACKAGE="abcrf")

  weights <- matrix(result$weightvec,nrow= nobs)
  weights.std <- sapply(1:nrow(x),function(x) weights[,x]/sum(weights[,x])) # weights std

  # esper <- sapply(1:nrow(x), function(x) weights.std[,x]%*%origObs)
  esperRf <- predict(object$model.rf, newdata=x, ...)

  # Out of bag expectations

  predict.oob <- predict(object$model.rf)

  # squared residuals

  residus.oob.sq <- (origObs - predict.oob)^2

  # variance estimation

  variance <- sapply(1:nrow(x), function(x) weights.std[,x] %*% residus.oob.sq)

  # Quantiles calculation

  ord <- order(origObs)
  origObs <- origObs[ord]
  weights <- weights[ord,,drop=FALSE]
  cumweights <- apply(weights,2,cumsum)
  cumweights <- sweep(cumweights,2,as.numeric(cumweights[nobs,]),FUN="/")

  # quantiles (code Meins)

  for (qc in 1:length(quantiles)){
    larg <- cumweights<quantiles[qc]
    wc <- apply(larg,2,sum)+1
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
  wc <- apply(larg,2,sum)+1
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

  tmp <- list(expectation = esperRf, med = mediane, variance = variance, quantiles = quant)
  class(tmp) <- "regAbcrfpredict"
  tmp

}


print.regAbcrfpredict <-
function(x, ...){

  ret <- cbind(x$expectation, x$med, x$variance, x$quantiles)
  colnames(ret) <- c("expectation", "median", "variance", colnames(x$quantiles) )
  print(ret, ...)

}

as.data.frame.regAbcrfpredict <-
function(x, ...) {

  ret <- cbind(x$expectation, x$med, x$variance, x$quantiles)
  colnames(ret) <- c("expectation", "median", "variance", colnames(x$quantiles) )
  as.data.frame(ret,  row.names=NULL, optional=FALSE, ...)

}

as.matrix.regAbcrfpredict <-
function(x, ...){

  ret <- cbind(x$expectation, x$med, x$variance, x$quantiles)
  colnames(ret) <- c("expectation", "median", "variance", colnames(x$quantiles) )
  ret

}

as.list.regAbcrfpredict <-
function(x, ...){

    list(expectation = x$expectation, med = x$med , variance = x$variance, quantiles=x$quantiles, ...)

}
