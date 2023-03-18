\name{densityPlot}
\alias{densityPlot}
\alias{densityPlot.regAbcrf}
\title{Plot the posterior density given a new summary statistic}
\description{
Given a reg-ABC-RF object and a new value of the summary statistics, \code{densityPlot} gives the corresponding posterior density plot of the parameter.
}
\usage{
\method{densityPlot}{regAbcrf}(object, newdata, main="Posterior density", ...)
}
\arguments{
  \item{object}{a \code{\link{regAbcrf}} object.}
  \item{newdata}{a new dataset of summary statistics.}
  \item{main}{main title to be used for the posterior density plot.}
  \item{\dots}{additional arguments to be passed on to \code{\link{density}}.}
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
\code{\link{plot.regAbcrf}}
}
\examples{
  data(snp)
  modindex <- snp$modindex
  sumsta <- snp$sumsta[modindex == "3",]
  r <- snp$param$r[modindex == "3"]
  r <- r[1:500]
  sumsta <- sumsta[1:500,]
  model.rf.r <- regAbcrf(r, sumsta, ntree=100)
  data(snp.obs)
  densityPlot(model.rf.r, snp.obs$sumsta[1,], main = "Posterior density of r")
}
