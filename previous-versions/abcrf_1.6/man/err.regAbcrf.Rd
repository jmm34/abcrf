\name{err.regAbcrf}
\alias{err.regAbcrf}
\title{Calculate and plot for different numbers of tree, the out-of-bag mean squared errors associated with a REG-ABC-RF object}
\description{
\code{err.regAbcrf} returns out-of-bag mean squared errors and plot them.
}
\usage{
err.regAbcrf(object, training, paral=FALSE, 
ncores= if(paral) max(detectCores()-1,1) else 1)
}
\arguments{
  \item{object}{a \code{regAbcrf} object.}
  \item{training}{the data frame containing the reference table used to train the \code{\link{regAbcrf}} object.}
  \item{paral}{a boolean that indicates if random forests predictions should be parallelized.}
  \item{ncores}{the number of CPU cores to use for the random forest predictions. If paral=TRUE, it is used the number of CPU cores minus 1.}
}
\value{
  \item{}{a matrix with 2 columns: the number of trees and the out-of-bag mean squared errors. NAs might be returned if the number of trees is too low. Errors are computed from 40 trees to the total number.}
}
\references{
Marin, J.-M., Raynal, L., Pudlo, P., Ribatet, M. and Robert, C.P. (2016)
\emph{ABC random forests for Bayesian parameter inference}
\url{http://arxiv.org/pdf/1605.05537v2.pdf}
}
\seealso{
  \code{\link{regAbcrf}},
  \code{\link{predict.regAbcrf}},
  \code{\link{plot.regAbcrf}},
  \code{\link{densityPlot}},
  \code{\link{covRegAbcrf}},
  \code{\link{ranger}}
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
  err.regAbcrf(model.rf.r, data2)
}