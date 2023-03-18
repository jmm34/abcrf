\name{abcrf}
\alias{abcrf}
\alias{abcrf.formula}
\alias{abcrf.factor}

\title{Create an ABC-RF object: a classification random forest from a reference table towards performing an ABC model choice}
\description{
\code{abcrf} constructs a random forest from a reference table towards performing
an ABC model choice. Basically, the reference table (i.e. the dataset that will
be treated with the present package) includes a column with the index
of the models to be compared and additional columns corresponding
to the values of the simulated summary statistics.
}
\usage{
\method{abcrf}{factor}(modindex, sumsta, lda=TRUE, ntree=500,
sampsize=min(1e5, length(modindex)), paral=FALSE, ...)
\method{abcrf}{formula}(formula, data=NULL, ...)
}
\arguments{
\item{modindex}{a factor containing the model indices of the reference table.}
\item{sumsta}{a matrix or a data frame containing the summary statistics
of the reference table.}
\item{lda}{should LDA scores be added to the list of summary statistics?}
\item{ntree}{number of trees to grow in the forest, by default 500 trees.}
\item{sampsize}{size of the sample from the reference table to grow a tree
of the classification forest, by default the minimum between the number of elements of the reference table and 100,000.}
\item{paral}{a boolean that indicates if the calculations of the classification random forest
(forest used to assign a model to the observed dataset) should be parallelized.
If yes, it is used the number of CPU cores minus 1.}
\item{formula}{a formula: left of ~, variable representing
the model index; right of ~, summary statistics.}
\item{data}{an optional data frame containing the reference table.}
\item{...}{additional arguments to be passed on to \code{\link{randomForest}} used to construct
the classification random forest that preditcs the selected model.}
}
\value{
An object of class \code{abcrf}, which is a list with the
following components:
\item{call}{the original call to \code{abcrf},}
\item{lda}{a boolean indicating if LDA scores have been added to the list of summary statistics,}
\item{model.rf}{an object of class \code{randomForest} containing the trained forest with the reference table,}
\item{model.lda}{an object of class \code{lda} containing the Linear Discriminant Analysis based on the reference table,}
\item{prior.err}{prior error rates of model selection on the reference table, estimated with the "out-of-bag" error of the forest,}
\item{sumsta}{the reference table including the LDA scores if lda=TRUE.}
}
\references{
Pudlo, P., Marin, J.-M., Estoup, A., Cornuet, J.-M., Gautier, M. and Robert, C.P. (2016)
\emph{Reliable ABC model choice via random forests} Bioinformatics
\url{http://bioinformatics.oxfordjournals.org/content/32/6/859}
}
\seealso{
\code{\link{plot.abcrf}},
\code{\link{predict.abcrf}},
\code{\link{err.abcrf}},
\code{\link{randomForest}}
}
\examples{
data(snp)
modindex <- snp$modindex[1:500]
sumsta <- snp$sumsta[1:500,]
model.rf <- abcrf(modindex, sumsta, ntree=100)
model.rf
}
