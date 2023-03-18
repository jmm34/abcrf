regAbcrf.numeric <-
  function(resp, sumsta, ntree=500, sampsize=min(1e5, length(resp)), paral=FALSE, ... ){
    if(!is.numeric(resp))
      stop("resp must be a numeric vector")
    sumsta <- as.matrix(sumsta)
    if(length(resp)!=nrow(sumsta))
      stop("length of resp differs from the number of lines in sumsta")
    if(length(resp)==0L)
      stop("no simulation in the reference table (resp, sumsta)")
    if (is.null(colnames(sumsta))) colnames(sumsta) <- paste("V",1:dim(sumsta)[2],sep="")
    if (paral==TRUE) {
      ncores <- max(detectCores()-1,1) 
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      if (trunc(ntree/ncores)==ntree/ncores) ntrees <- rep(ntree/ncores, ncores) else
        ntrees <- c(rep(trunc(ntree/ncores), ncores),ntree-trunc(ntree/ncores)*ncores)
      model.rf <- foreach(ntree=ntrees, .combine= combine, .multicombine=TRUE, .packages='randomForest') %dorng% {
        randomForest(sumsta, resp, ntree=ntree, sampsize=sampsize, keep.inbag=TRUE, ...)
      }
      stopCluster(cl)
      pred.noob <- predict(model.rf, newdata=sumsta, predict.all=TRUE)
      mat <- pred.noob$individual
      for( j in 1:model.rf$ntree ){
        mat[model.rf$inbag[,j]!=0,j] <- NA
      }
      model.rf$predicted <- sapply(1:nrow(sumsta), function(x) mean(mat[x,!is.na(mat[x,])]) )
    } else model.rf <- randomForest(sumsta, resp, ntree=ntree, sampsize=sampsize, keep.inbag = TRUE, ...)
    cl <- match.call()
    cl[[1]] <- as.name("regAbcrf")
    x <- list(call=cl, model.rf = model.rf, sumsta=sumsta)
    class(x) <- "regAbcrf"
    x
  }

regAbcrf <-
function(...) UseMethod("regAbcrf")

print.regAbcrf <-
function(x, ...){

  cat("\nCall:\n", deparse(x$call), "\n")
  cat("Number of simulations: ", length(x$model.rf$y), "\n", sep="")
  cat("Number of trees: ", x$model.rf$ntree, "\n", sep="")
  cat("No. of variables tried at each split: ", x$model.rf$mtry, "\n", sep="")
      
}