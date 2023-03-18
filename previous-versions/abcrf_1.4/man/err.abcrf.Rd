\name{err.abcrf}
\alias{err.abcrf}
\title{Calculate and plot for different numbers of tree, the out-of-bag errors associated with an ABC-RF object}
\description{
\code{err.abcrf} returns out-of-bag errors and plot them.
}
\usage{
err.abcrf(object, training, paral=FALSE, ncores= if(paral) max(detectCores()-1,1) else 1)
}
\arguments{
\item{object}{an \code{abcrf} object.}
\item{training}{the data frame containing the reference table used to train the \code{\link{abcrf}} object.}
\item{paral}{a boolean that indicates if random forests predictions should be parallelized.}
\item{ncores}{the number of CPU cores to use for the random forest predictions. If paral=TRUE, it is used the number of CPU cores minus 1.}
}
\value{
\item{}{a matrix with 2 colunms: the number of trees and the out-of-bag errors}
}
\references{
Pudlo, P., Marin, J.-M., Estoup, A., Cornuet, J.-M., Gautier, M. and Robert, C.P. (2016)
\emph{Reliable ABC model choice via random forests} Bioinformatics
\url{http://bioinformatics.oxfordjournals.org/content/32/6/859}
}
\seealso{
\code{\link{abcrf}},
\code{\link{predict.abcrf}},
\code{\link{plot.abcrf}}
}
\examples{
data(snp)
modindex <- snp$modindex[1:500]
sumsta <- snp$sumsta[1:500,]
data1 <- data.frame(modindex, sumsta)
model.rf <- abcrf(modindex~., data1, ntree=100)
err.rf <- err.abcrf(model.rf, data1)
}
