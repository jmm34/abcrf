\name{regAbcrf}
\alias{regAbcrf}
\alias{regAbcrf.numeric}
\title{
Create a reg-ABC-RF object: a regression random forest from a reference table aimed out predicting posterior expectation, variance and quantiles for a parameter}
\description{
\code{regAbcrf} constructs a regression random forest from a reference table towards predicting posterior expectations, variances and quantiles of a parameter.
}
\usage{
\method{regAbcrf}{numeric}(resp, sumsta, ntree=500, sampsize=min(1e5, length(resp)), paral=FALSE, ... )
}
\arguments{
  \item{resp}{a numeric vector containing the response variable of the reference table.}
  \item{sumsta}{a matrix or a data frame containing the summary statistics of the reference table.}
  \item{sampsize}{size of the sample from the reference table used to grow a tree of the regression forest, by default the minimum between the number of elements of the reference table and 100,000.}
  \item{ntree}{number of trees to grow in the forest, by default 500 trees.}
   \item{paral}{a boolean that indicates whether or not the calculations of the regression random forest (forest used to predict a response from the observed dataset) should be parallelized. If yes, it is the number of CPU cores minus 1.}
   \item{\dots}{additional arguments to be passed on to \code{\link{randomForest}} used to construct
the regression random forest that predicts the response variable.}
}
\value{
An object of class \code{regAbcrf}, which is a list with the
following components:
\item{call}{the original call to \code{regAbcrf},}
\item{model.rf}{an object of class \code{randomForest} containing the trained forest with the reference table,}
\item{sumsta}{the summary statistics of the reference table.}
}
\references{
Marin, J.-M., Raynal, L., Pudlo, P., Ribatet, M. and Robert, C.P. (2016)
\emph{ABC random forests for Bayesian parameter inference}
\url{http://arxiv.org/pdf/1605.05537v1.pdf}
}

\seealso{
\code{\link{plot.regAbcrf}},
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
  model.rf.r
}
