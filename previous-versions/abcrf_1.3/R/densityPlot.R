densityPlot.regAbcrf <- 
function(object, newdata, main="Posterior density", ...){
    
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
    
    #####################
    
    origObs <- obj$origObs
    origNodes <- obj$origNodes
    
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
    
    plot(density(obj$y, weights=weights.std, ...), main=main )
    
}

densityPlot <-
  function(...) UseMethod("densityPlot")