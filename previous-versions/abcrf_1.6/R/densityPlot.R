densityPlot.regAbcrf <- 
function(object, obs, training,  main="Posterior density", log="", paral=FALSE, ncores= if(paral) max(detectCores()-1,1) else 1, ...)
{
    ### Checking arguments
    if (!inherits(object, "regAbcrf")) 
      stop("object not of class regAbcrf")
  
    if (!inherits(training, "data.frame"))
      stop("training needs to be a data.frame object")
  
    if (!inherits(obs, "data.frame")) 
      stop("obs needs to be a data.frame object")
    if (nrow(obs) == 0L || is.null(nrow(obs)))
      stop("no data in obs")
    if (nrow(training) == 0L || is.null(nrow(training)))
      stop("no simulation in the training reference table (response, sumstat)")
    if ( (!is.logical(paral)) && (length(paral) != 1L) )
      stop("paral should be TRUE or FALSE")
    if( ncores > detectCores() || ncores < 1 )
      stop("incorrect number of CPU cores")
    if( !is.character(log) )
      stop("log needs to be a character string")
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

    # resp and sumsta recover
  
    mf <- match.call(expand.dots=FALSE)
    mf <- mf[1]
    mf$formula <- object$formula
    mf$data <- training
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame() )
    mt <- attr(mf, "terms")
    resp <- model.response(mf)
    
    obj <- object$model.rf
    inbag <- matrix(unlist(obj$inbag.counts, use.names=FALSE), ncol=obj$num.trees, byrow=FALSE)
    
    obj[["origNodes"]] <- predict(obj, training, predict.all=TRUE, num.threads=ncores)$predictions
    obj[["origObs"]] <- model.response(mf)
    
    #####################

    origObs <- obj$origObs
    origNodes <- obj$origNodes
    
    nodes <- predict(obj, x, predict.all=TRUE, num.threads=ncores )$predictions
    ntree <- obj$num.trees
    nobs <- object$model.rf$num.samples
    nnew <- nrow(x)

    weights <- findweights(origNodes, inbag, nodes, as.integer(nobs), as.integer(nnew), as.integer(ntree)) # cpp function call
    weights.std <- weights/ntree

    for(i in 1:nnew){
      plot(density(resp, weights=weights.std[,i], ...), main=main, log=log)
      if(nnew>1 && i<nnew) readline("Press <ENTER> to Continue")
    }
}

densityPlot <-
  function(...) UseMethod("densityPlot")