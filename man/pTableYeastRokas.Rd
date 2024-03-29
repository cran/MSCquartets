% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSCquartets-package.R
\docType{data}
\name{pTableYeastRokas}
\alias{pTableYeastRokas}
\title{pTable for Yeast dataset}
\format{
an R data file
}
\source{
\url{https://wiki.rice.edu/confluence/download/attachments/8898533/yeast.trees?version=1&modificationDate=1360603275797&api=v2}
}
\usage{
data(pTableYeastRokas)
}
\description{
An .rda file dataset for the "dataYeastRokas" dataset. This is a subset of the data
of \insertCite{Rokas03;textual}{MSCquartets}.
}
\details{
This is provided primarily so that examples of other functions run more quickly. It can be reproduced by the following
example code below.
}
\examples{
\donttest{
gtrees=read.tree(file = system.file("extdata","dataYeastRokas",package="MSCquartets"))
QT=quartetTable(gtrees)
RQT=quartetTableResolved(QT)
pTable=quartetCutTestInd(RQT)
pTable=quartetTreeTestInd(pTable)
pTableYeastRokas=quartetStarTestInd(pTable)
}
}
\references{
\insertRef{Rokas03}{MSCquartets}
}
\keyword{datasets}
