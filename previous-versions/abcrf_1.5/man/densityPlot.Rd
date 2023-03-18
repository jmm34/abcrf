\name{densityPlot}
\alias{densityPlot}
\alias{densityPlot.regAbcrf}
\title{Plot the posterior density given a new summary statistic}
\description{
Given a reg-ABC-RF object and a new value of the summary statistics, 
\code{densityPlot} gives the corresponding posterior density plot of the parameter.
}
\usage{
\method{densityPlot}{regAbcrf}(object, obs, training,
main="Posterior density", paral=FALSE, 
ncores= if(paral) max(detectCores()-1,1) else 1, ...)
}
\arguments{
  \item{object}{a \code{\link{regAbcrf}} object.}
  \item{obs}{a data frame containing the summary statistics of the observed data sets.}
  \item{training}{the data frame containing the reference table used to train the \code{\link{regAbcrf}} object.}
  \item{main}{main title to be used for the posterior density plot.}
  \item{paral}{a boolean that indicates if random forests predictions should be parallelized.}
  \item{ncores}{the number of CPU cores to use for the regression random forest predictions. If paral=TRUE, it is used the number of CPU cores minus 1.}
  \item{\dots}{additional arguments to be passed on to \code{\link{density}}.}
}
\references{
Marin, J.-M., Raynal, L., Pudlo, P., Ribatet, M. and Robert, C.P. (2016)
\emph{ABC random forests for Bayesian parameter inference}
\url{http://arxiv.org/pdf/1605.05537v2.pdf}
}
\seealso{
\code{\link{regAbcrf}},
\code{\link{predict.regAbcrf}},
\code{\link{err.regAbcrf}},
\code{\link{covRegAbcrf}},
\code{\link{ranger}},
\code{\link{plot.regAbcrf}}
}
\examples{
data(snp)
modindex <- snp$modindex
sumsta <- snp$sumsta[modindex == "3",]
r <- snp$param$r[modindex == "3"]
r <- r[1:500]
sumsta <- sumsta[1:500,]
data2 <- data.frame(r, sumsta)
model.rf.r <- regAbcrf(r~., data2, ntree=100)
data(snp.obs)
densityPlot(model.rf.r, snp.obs, data2, main = "Posterior density of r")
}
