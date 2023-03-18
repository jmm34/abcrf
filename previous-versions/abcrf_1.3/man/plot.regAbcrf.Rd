\name{plot.regAbcrf}
\alias{plot.regAbcrf}
\title{Plot of a reg-ABC-RF object}
\description{
\code{plot.regAbcrf} provides a variable importance plot used to construct the reg-ABC-RF object.
}
\usage{
\method{plot}{regAbcrf}(x, n.var = min(30, nrow(x$model.rf$importance)), ...)
}
\arguments{
  \item{x}{a \code{\link{regAbcrf}} object.}
  \item{n.var}{number of variables in the variable importance representation. The default value is equal to the minimum between 30 and the number of summary statistics.}
 \item{\dots}{additional arguments to be passed on to \code{\link{varImpPlot}}.}
}
\references{
Marin, J.-M., Raynal, L., Pudlo, P., Ribatet, M. and Robert, C.P. (2016)
\emph{ABC random forests for Bayesian parameter inference}
\url{http://arxiv.org/pdf/1605.05537v1.pdf}
}
\seealso{
\code{\link{regAbcrf}},
\code{\link{predict.regAbcrf}},
\code{\link{covRegAbcrf}},
\code{\link{randomForest}},
\code{\link{densityPlot}}
}
\examples{
  data(snp)
  modindex <- snp$modindex
  sumsta <- snp$sumsta[modindex == "3",]
  r <- snp$param$r[modindex == "3"]
  r <- r[1:500]
  sumsta <- sumsta[1:500,]
  model.rf.r <- regAbcrf(r, sumsta, ntree=100)
  plot(model.rf.r)
}
