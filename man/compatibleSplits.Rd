% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TreeOfBlobs.R
\name{compatibleSplits}
\alias{compatibleSplits}
\title{Extract compatible splits}
\usage{
compatibleSplits(sp, tol = 0, plot = FALSE)
}
\arguments{
\item{sp}{an object of class splits}

\item{tol}{splits with weights below tol are dropped}

\item{plot}{a logical; if TRUE plots tree displaying remaining spilts}
}
\value{
splits objects containing only those that are compatible and high weight
}
\description{
Given an object of class splits, first discards any with weight less than a tolerance, and
then further removes all remaining splits that are incompatible with any other remaining one.
}
\examples{
data(pTableYeastRokas)
dist=NANUQdist(pTableYeastRokas, alpha=.05, beta=.95,outfile=NULL)
nn=neighborNet(dist)
plot(nn,"2D")
tob=treeFromSplits(compatibleSplits(nn$splits),plot=TRUE) #produce tree of blobs of splits graph

}
\seealso{
\code{\link{treeFromSplits}, \link{TINNIK}}
}
