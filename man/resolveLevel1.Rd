% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ToBtoLevel1.R
\name{resolveLevel1}
\alias{resolveLevel1}
\title{Resolve Tree of Blobs to Level-1 network}
\usage{
resolveLevel1(
  ToB,
  pTable,
  test = "T3",
  alpha,
  beta,
  distance = "NANUQ",
  hdegree = 10,
  plot = 2,
  delta = 10^-6,
  fullResMax = 10
)
}
\arguments{
\item{ToB}{an unrooted tree of blobs (phylo object) as determined by TINNIK or
another method}

\item{pTable}{a table of qcCFs, with columns p_star and p_test}

\item{test}{either "T3" or "cut", indicating test to use for determining what
qcCFs indicate hybridization}

\item{alpha}{test level for p_test}

\item{beta}{test level value for p_star}

\item{distance}{cycle resolution distance to be used ("NANUQ" or "modNANUQ")}

\item{hdegree}{resolve a multifurcation of this degree or larger by a heuristic
method; must be at least 5}

\item{plot}{if 0, no plots; if 1, plot only possible root locations on ToB and
full resolution; if 2, include plots of each individual blob resolution, if 3
include histograms of measure of fit for all hybrid/orders considered in choosing
best}

\item{delta}{cutoff for relative difference in squared residuals and smallest,
(RSS-minRSS)/minRSS, for determining near ties as "best" fit resolutions}

\item{fullResMax}{maximum number of full resolutions (all multifurcations at once)
to form; if the product of the number of resolutions of individual multifurcations
exceeds this, no full resolutions are produced, although \code{\link{combineCycleResolutions}}
can be applied to produce them.}
}
\value{
a list of resolutions and squared residuals:
\itemize{
\item [[1]] is a list of Newick
resolutions of entire network, with all edge lengths 1 (NULL if one cannot be produced or \code{fullResMax} is exceeded),
\item [[2]]-[[n]] are individual resolutions of each  multifurcation on \code{ToB},
each given as a list as output from \code{\link{resolveCycle}}.
}
}
\description{
Given a Tree of Blobs and qcCF information, resolve all multifurcations to
cycles. Resolution is performed by finding a least-squares best-fit of an
empirical distance to an expected distance related to the cycle, as described
in \insertCite{ABRW24;textual}{MSCquartets}.
}
\details{
Possible distances to use are the NANUQ distance and a modified NANUQ distance
that weights quartet trees differently, but from which the cycle structure is
still identifiable.

For multifucations of degree less than a designated cutoff, all possible
circular orders and choices of hybrid nodes are considered in choosing the
best. Above that cutoff, a heuristic method  is used to obtain
a small number of orders likely to be good fits, with the least-squares
fitting applied only to those.
}
\examples{
data(pTableYeastRokas)
out=TINNIK(pTableYeastRokas, alpha=.01, beta=.05)
ToB=labelIntNodes(out$ToB)
resolveLevel1(ToB, pTable=out$pTable, alpha=.01, beta=.05, distance="NANUQ")

}
\references{
\insertRef{ABRW24}{MSCquartets}
}
\seealso{
\code{\link{TINNIK}}, \code{\link{labelIntNodes}}, \code{\link{resolveCycle}},
\code{\link{combineCycleResolutions}}
}
