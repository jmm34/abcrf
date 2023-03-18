\name{predict.regAbcrf}
\alias{predict.regAbcrf}
\title{
Predict posterior expectation, median, variance and quantiles
given a new dataset using a reg-ABC-RF object
}
\description{
Based on a reg-ABC-RF object this function predicts the posterior expectation, median, variance, quantiles for the corresponding parameter given new dataset.
}
\usage{
\method{predict}{regAbcrf}(object, newdata, quantiles = c(0.025, 0.975), ...)
}
\arguments{
  \item{object}{a \code{\link{regAbcrf}} object.}
  \item{newdata}{a matrix or data frame containing the summary statistics of the observed data sets.}
  \item{quantiles}{numeric vector of probabilities with values in [0,1]. The default value is equal to \code{c(0.025, 0.975)}.}
  \item{\dots}{optional arguments to be passed on to the function \code{\link{predict.randomForest}}.}
}

\value{
An object of class \code{regAbcrfpredict}, which is a list the following components:
\item{expectation}{predicted posterior expectation for each oberved data set,}
\item{med}{predicted posterior median for each oberved data set,}
\item{variance}{predicted posterior variance for each observed data set, computed by reusing weights,}
\item{quantiles}{predicted posterior quantiles for each observed data set.}
}
\references{
Marin, J.-M., Raynal, L., Pudlo, P., Ribatet, M. and Robert, C.P. (2016)
\emph{ABC random forests for Bayesian parameter inference}
\url{http://arxiv.org/pdf/1605.05537v1.pdf}
}
\seealso{
\code{\link{regAbcrf}},
\code{\link{plot.regAbcrf}},
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
  data(snp.obs)
  predict(model.rf.r, snp.obs$sumsta)
}
