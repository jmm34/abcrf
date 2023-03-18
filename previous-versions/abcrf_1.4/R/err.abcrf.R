err.abcrf <- function(object, training, paral=FALSE, ncores= if(paral) max(detectCores()-1,1) else 1)
{
  if (!inherits(training, "data.frame"))
    stop("training needs to be a data.frame object")
  if ( (!is.logical(paral)) && (length(paral) != 1L) )
    stop("paral should be TRUE or FALSE")
  if(ncores > detectCores() || ncores < 1)
    stop("incorrect number of CPU cores")
  
  nmod <- length(object$model.rf$forest$levels)
  ntrain <- nrow(training)
   
  mf <- match.call(expand.dots=FALSE)
  mf <- mf[1]
  mf$formula <- object$formula
  mf$data <- training
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame() )
  mt <- attr(mf, "terms")
  modindex <- model.response(mf)
  
  if (object$lda) {
    training <- cbind( training, predict(object$model.lda, training)$x )
  }
  inbag <- matrix(unlist(object$model.rf$inbag.counts, use.names=FALSE), ncol=object$model.rf$num.trees, byrow=FALSE)
  
  mimi <- predict(object$model.rf, training, predict.all=TRUE, num.threads=ncores)$predictions
  if (object$model.rf$num.trees < 10) stop("the number of trees in the forest should be greater than 10")
  sequo <- seq(10,object$model.rf$num.trees, length.out = 20)

  res <- rep(0,length(sequo))
  h <- 0
  pb=txtProgressBar(min=0,max=length(sequo),style=3)
  for (j in sequo) 
  {
    mama <- sample(paste(1:nmod),ntrain,replace=TRUE)
    for (i in 1:ntrain) 
    {
      outbag <- (1:j)[inbag[i,1:j]==0]
      if (length(outbag)>0) mama[i] <- names(which.max(table(as.factor(mimi[i,outbag]))))
    }
  h <- h+1
  res[h] <- mean(as.factor(mama)!=modindex)
  setTxtProgressBar(pb,h)
  }
  close(pb)

  plot(sequo,res,ylab="Prior error rate",xlab="Number of trees",type="l")
  cbind(ntree=floor(sequo), error.rate=res)
}
