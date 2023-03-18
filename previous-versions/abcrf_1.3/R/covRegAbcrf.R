covRegAbcrf.regAbcrf <-
function(regForest1, regForest2, newdata, ntree=500, sampsize=min(1e5, nrow(regForest1$sumsta)), paral=FALSE, ... ){
  
  if (!inherits(regForest1, "regAbcrf")) 
    stop("regForest1 not of class regAbcrf")
  if (!inherits(regForest2, "regAbcrf")) 
    stop("regForest2 not of class regAbcrf")

  if(any(dim(regForest1$sumsta) != dim(regForest2$sumsta) ))
    stop("regForest1 and regForest2 training data are not constructed on the same summaries")
  if(any(regForest1$sumsta != regForest2$sumsta) )
    stop("regForest1 and regForest2 training data are not constructed on the same summaries")
  
  train.data <- regForest1$sumsta
  object1 <- regForest1$model.rf
  object2 <- regForest2$model.rf
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
  
  if(!is.null(x)){
    vname <- if (is.null(dim(object2$importance))) {
      names(object2$importance)
    } else {
      rownames(object2$importance)
    }
    
    if (any(colnames(x) != vname))
      stop("names of predictor variables do not match")
  }
  
  
  # residuals
  
  res1 <- regForest1$model.rf$y - predict(regForest1$model.rf)
  res2 <- regForest2$model.rf$y - predict(regForest2$model.rf)
  
  res12 <- res1*res2 # new response varible
  
  # forest construction
  
  if (paral==TRUE) {
    ncores <- max(detectCores()-1,1) 
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    if (trunc(ntree/ncores)==ntree/ncores) ntrees <- rep(ntree/ncores, ncores) else
      ntrees <- c(rep(trunc(ntree/ncores), ncores),ntree-trunc(ntree/ncores)*ncores)
    model.rf <- foreach(ntree=ntrees, .combine= combine, .multicombine=TRUE, .packages='randomForest') %dorng% {
      randomForest(regForest1$sumsta, res12, ntree=ntree, sampsize=sampsize, keep.inbag=TRUE, ...)
    }
    stopCluster(cl)
    pred.noob <- predict(model.rf, newdata=regForest1$sumsta, predict.all=TRUE)
    mat <- pred.noob$individual
    for( j in 1:model.rf$ntree ){
      mat[model.rf$inbag[,j]!=0,j] <- NA
    }
    model.rf$predicted <- sapply(1:nrow(regForest1$sumsta), function(x) mean(mat[x,!is.na(mat[x,])]) )
  } else model.rf <- randomForest(regForest1$sumsta, res12, ntree=ntree, sampsize=sampsize, keep.inbag = TRUE, ...)
  
  result <- predict(model.rf, newdata=x)
  
  return(result)
  
}


covRegAbcrf <-
function(...) UseMethod("covRegAbcrf")
