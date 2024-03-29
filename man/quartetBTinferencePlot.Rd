% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TreeOfBlobs.R
\name{quartetBTinferencePlot}
\alias{quartetBTinferencePlot}
\title{Produce simplex plot with  results of B/T-quartet inference}
\usage{
quartetBTinferencePlot(pTable, Bquartets, test, alpha, beta, cex = 1)
}
\arguments{
\item{pTable}{table of quartets and p-values}

\item{Bquartets}{indicator vector for B-quartets (1=B, 0=T), ordered as in pTable}

\item{test}{test model used for tree null hypothesis; options are \code{"cut"}, \code{"T3"}}

\item{alpha}{significance level used by TINNIK for test \code{test}}

\item{beta}{significance level used by TINNIK for star tree test}

\item{cex}{scaling factor for size of plotted symbols}
}
\description{
Plot a 2-d probability simplex, with points for all normalized quartet count
vectors. Colors indicate B- or T-quartets from TINNIK algorithm, at specified
test levels.
}
\details{
The first argument of this function is a table of quartets and p-values. The
plot may show results using either the T3, or 2-cut
test, and a star tree test.
The p-values must be computed by or before previous calls to
\code{TINNIK}. The second argument is the indicator vector for B/T quartets produced by \code{TINNIK}.
}
\examples{
data(pTableYeastRokas)
out=TINNIKdist(pTableYeastRokas,test="T3",alpha=.05,beta=.05)
quartetBTinferencePlot(pTableYeastRokas,out$B,test="T3",alpha=.05,beta=.05)

}
\seealso{
\code{TINNIK}, \code{\link{quartetTestPlot}}
}
