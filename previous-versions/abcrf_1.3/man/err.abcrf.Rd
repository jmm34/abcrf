\name{err.abcrf}
\alias{err.abcrf}
\title{Calculate and plot for different numbers of tree, the out-of-bag errors associated with an ABC-RF object}
\description{
\code{err.abcrf} returns out-of-bag errors and plot them.
}
\usage{
err.abcrf(x)
}
\arguments{
\item{x}{an \code{abcrf} object.}
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
model.rf <- abcrf(modindex, sumsta, ntree=100)
err.rf <- err.abcrf(model.rf)
}
