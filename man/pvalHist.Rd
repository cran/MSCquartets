% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pHist.R
\name{pvalHist}
\alias{pvalHist}
\title{Plot histogram of log p-values in table}
\usage{
pvalHist(pTable, model, pmin = 0)
}
\arguments{
\item{pTable}{a quartet table with p-values, such as from \code{NANUQ}, \code{quartetTreeTestInd},
or \code{quartetStarTestInd}}

\item{model}{one of \code{"T1"}, \code{"T3"}, or \code{"star"}, where \code{pTable}
contains a column \code{p_model} of p-values}

\item{pmin}{include only p-values above \code{pmin} in the histogram}
}
\description{
Graphical exploration of extreme p-values from quartet hypothesis tests, to aid in choosing critical values
for hypothesis tests. Log base 10 of p-values exceeding some minimum are plotted, to explore gaps in
the tail of the distribution.
}
\details{
Since logarithms are plotted, p-values close to 0 will appear as negative numbers of large magnitude, putting the tail of the distribution
to the left in the histogram.

When exploring possible critical values for the hypothesis tests in the NANUQ algorithm, use \code{model= "T3"} to
choose values for \code{alpha} which distinguish treelikeness from hybridization, and \code{model= "star"} to
choose values for \code{beta} which distinguish polytomies from resolved trees.
In general, \code{alpha} should be chosen to be small and \code{beta}
to be large so that most quartets are interpreted as resolved trees.
}
\examples{
data(pTableYeastRokas)
pvalHist(pTableYeastRokas,"T3")

}
\seealso{
\code{\link{NANUQ}}, \code{\link{NANUQdist}}
}
