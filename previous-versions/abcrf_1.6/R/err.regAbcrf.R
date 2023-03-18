err.regAbcrf <- function(object, training, paral=FALSE, ncores= if(paral) max(detectCores()-1,1) else 1)
{
  if (!inherits(training, "data.frame"))
    stop("training needs to be a data.frame object")
  if ( (!is.logical(paral)) && (length(paral) != 1L) )
    stop("paral should be TRUE or FALSE")
  if(ncores > detectCores() || ncores < 1)
    stop("incorrect number of CPU cores")
  
  ntrain <- nrow(training)
  ntree <- object$model.rf$num.trees
  
  mf <- match.call(expand.dots=FALSE)
  mf <- mf[1]
  mf$formula <- object$formula
  mf$data <- training
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame() )
  mt <- attr(mf, "terms")
  resp <- model.response(mf)
  
  inbag <- matrix(unlist(object$model.rf$inbag.counts, use.names=FALSE), ncol=object$model.rf$num.trees, byrow=FALSE)
  
  pred <- predict(object$model.rf, training, predict.all=TRUE , num.threads=ncores)$predictions
  
  if (ntree < 40) stop("the number of trees in the forest should be greater than 40")
  sequo <- seq(40,object$model.rf$num.trees, length.out = 20)
  
  res <- oobErrorsReg(as.integer(floor(sequo)), as.integer(ntrain), as.integer(ntree), as.numeric(resp), inbag, pred)
  
  plot(floor(sequo), res, ylab="out-of-bag mean squared error",xlab="Number of trees",type="l")
  cbind(ntree=floor(sequo), oob_mse=res)
}