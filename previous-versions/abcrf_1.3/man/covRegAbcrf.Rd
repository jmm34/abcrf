\name{covRegAbcrf}
\alias{covRegAbcrf}
\alias{covRegAbcrf.regAbcrf}
\title{
Predict posterior covariance between two parameters for new data using two reg-ABC-RF objects
}
\description{
Using two reg-ABC-RF objects constructed on the same reference table for two different response variables, this function predicts the posterior covariance between those two response variables, given a new dataset of summaries.
}
\usage{
\method{covRegAbcrf}{regAbcrf}(regForest1, regForest2, newdata, ntree=500,
            sampsize=min(1e5, nrow(regForest1$sumsta)), paral=FALSE, ... )
}
\arguments{
  \item{regForest1, regForest2}{\code{\link{regAbcrf}} objects.}
  \item{newdata}{a matrix or data frame containing the summary statistics of the observed data sets.}
  \item{ntree}{number of trees to grow in the forest, by default equal to 500 trees.}
  \item{sampsize}{size of the sample from the reference table used to grow a tree of the regression forest, by default the minimum between the number of elements of the reference table and 100,000.}
  \item{paral}{a boolean that indicates whether or not the calculations of the regression random forest (forest used to predict a response from the observed dataset) should be parallelized. If yes, it is the number of CPU cores minus 1.}
  \item{\dots}{additional arguments to be passed on to \code{\link{randomForest}} used to construct the regression random forest that predicts posterior covariance.}
}
\value{
\code{covRegAbcrf} returns predicted posterior covariances between response variables of two reg-ABC-RF objects, for a new data set.
}
\references{
Marin, J.-M., Raynal, L., Pudlo, P., Ribatet, M. and Robert, C.P. (2016)
\emph{ABC random forests for Bayesian parameter inference}
\url{http://arxiv.org/pdf/1605.05537v1.pdf}
}
\seealso{
\code{\link{regAbcrf}},
\code{\link{predict.regAbcrf}},
\code{\link{plot.regAbcrf}},
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
  N1 <- snp$param$N1[modindex == "3"]
  N1 <- N1[1:500]
  model.rf.N1 <- regAbcrf(N1, sumsta, ntree=100)
  data(snp.obs)
  covRegAbcrf(model.rf.r, model.rf.N1, snp.obs$sumsta, ntree=100)
}
