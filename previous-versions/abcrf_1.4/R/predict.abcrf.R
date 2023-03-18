predict.abcrf <- function(object, obs, training, ntree = 1000, sampsize = min(1e5,  object$model.rf$num.samples ), 
                          paral = FALSE, ncores = if(paral) max(detectCores()-1,1) else 1,
                          paral.predict = FALSE, ncores.predict = if(paral.predict) max(detectCores()-1,1) else 1 , ...)
{
  if (!inherits(obs, "data.frame")) 
    stop("obs needs to be a data.frame object")
  if (!inherits(training, "data.frame"))
    stop("training needs to be a data.frame object")
  if (nrow(obs) == 0L || is.null(nrow(obs)))
    stop("no data in obs")
  if (nrow(training) == 0L || is.null(nrow(training)))
    stop("no simulation in the training reference table (response, sumstat)")
  if ( (!is.logical(paral)) && (length(paral) != 1L) )
    stop("paral should be TRUE or FALSE")
  if ( (!is.logical(paral.predict)) && (length(paral.predict) != 1L) )
    stop("paral.predict should be TRUE or FALSE")
  if( (ncores > detectCores() || ncores < 1) || (ncores.predict > detectCores() || ncores.predict < 1) )
    stop("incorrect number of CPU cores")
  
  nmod <- length(object$model.rf$forest$levels)
  ntrain <- object$model.rf$num.samples
  nstat <- object$model.rf$num.independent.variables
  
  # modindex and sumsta recovery
  
  mf <- match.call(expand.dots=FALSE)
  mf <- mf[1]
  mf$formula <- object$formula
  mf$data <- training
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame() )
  mt <- attr(mf, "terms")
  modindex <- model.response(mf)
  sumsta <- model.frame(terms(reformulate(attributes(mt)$term.labels)), data.frame(mf))
  
  nobs <- nrow(obs)

	if (is.null(colnames(obs))) {
	  colnames(obs) <- object$model.rf$forest$independent.variable.names[1:(nstat-nmod-1)]
    warning("Columns of obs have no names")
	}
  if(sampsize > ntrain)
    stop("sampsize too large")
  
	old.options <- options(); options(warn=-1)
	if (object$lda) {
	  obs <- cbind(obs, predict(object$model.lda, obs)$x)
		sumsta <- cbind( sumsta, predict(object$model.lda, training)$x )
	}
	
	allocation <- predict(object$model.rf, obs, num.threads=ncores.predict)$predictions
	
	# vote :
	
	vote <- matrix( nrow=nobs, ncol=nmod)
	colnames(vote) <- object$model.rf$forest$levels
	
	pred.all <- predict(object$model.rf, obs, predict.all=TRUE, num.threads=ncores.predict)$predictions
	for(i in object$model.rf$forest$class.values ){
	  vote[,i] <- sapply(1:nobs, function(x) mean(pred.all[x,] == i ) )
	}
	
	local.error <- as.numeric(object$model.rf$predictions == modindex)

		
	# data.frame for ranger 
	data.ranger <- data.frame(local.error, sumsta )
	
  error.rf <- ranger(local.error~., data=data.ranger, num.trees = ntree, 
                       sample.fraction=sampsize/ntrain, num.threads = ncores, ...)

  options(old.options)
	tmp <- list(allocation=allocation, vote=vote, post.prob=predict(error.rf, obs, num.threads=ncores.predict)$predictions )
  class(tmp) <- "abcrfpredict"
  tmp
}

summary.abcrfpredict <- function(object, ...) {
  cat("Number of affectations per model:\n")
  summary(object$allocation, ...)
}

print.abcrfpredict <- function(x, ...) {
ret <- cbind.data.frame(x$allocation, x$vote, x$post.prob)
  colnames(ret) <- c("selected model", paste("votes model",1:dim(x$vote)[2],sep=""), "post.proba")
  print(ret, ...)
}

as.matrix.abcrfpredict <- function(x, ...) {
  ret <- cbind(x$allocation, x$vote, x$post.prob)
  colnames(ret) <- c("selected model", paste("votes model",1:dim(x$vote)[2],sep=""), "post.proba")
  ret
}

as.data.frame.abcrfpredict <- function(x, ...) {
  ret <- cbind(x$allocation, x$vote, x$post.prob)
  colnames(ret) <- c("selected model", paste("votes model",1:dim(x$vote)[2],sep=""), "post.proba")
  as.data.frame(ret,  row.names=NULL, optional=FALSE, ...)
}

as.list.abcrfpredict <- function(x, ...) {
  list(allocation = x$allocation, vote = x$vote, post.prob = x$post.prob, ...)
}