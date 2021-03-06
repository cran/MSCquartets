% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TreeQuartetDist.R
\name{WQDC}
\alias{WQDC}
\title{Compute Weighted Quartet Distance Consensus tree from gene tree data}
\usage{
WQDC(
  genetreedata,
  taxanames = NULL,
  method = fastme.bal,
  omit = FALSE,
  terminal = 1
)
}
\arguments{
\item{genetreedata}{gene tree data that may be supplied in any of 3 forms:
\enumerate{
\item a character string giving the name of a file containing gene trees in Newick
\item a multiPhylo object containing gene trees
\item a resolved quartet table, as produced by \code{quartetTableResolved}
}}

\item{taxanames}{if \code{genetreedata} is a file or a multiPhylo object, a vector of a subset
of the taxa names on the gene trees 
to be analyzed, if \code{NULL} all taxa on the first gene tree are used; if \code{genetreedata} 
is a quartet table, this argument is ignored and all taxa in the table are used}

\item{method}{a distance-based tree building function, such as \code{fastme.bal} or \code{nj}}

\item{omit}{\code{TRUE} leaves out unresolved quartets, \code{FALSE} treats them as 1/3 of each resolution; ignored if 
\code{genetreedata} is given as  a resolved quartet table}

\item{terminal}{non-negative branch length to supply for terminal branches 
whose length cannot be inferred by \code{WQDC}}
}
\value{
an unrooted metric tree of type phylo
}
\description{
Compute the Weighted Quartet Distance Consensus \insertCite{YR19}{MSCquartets} estimate of a 
species tree from gene tree data. This is a consistent estimator of the unrooted 
species tree topology and all internal branch lengths.
}
\details{
This function is a wrapper which performs the steps of reading in a collection
of gene trees, tallying quartets, estimating quartet internal branch lengths, computing the weighted
quartet distance between taxa, building
a tree, and adjusting edge lengths, to give a consistent estimate of the metric species tree in coalescent units
under the MSC.

If the gene tree data indicates some quartets experienced little to no incomplete lineage 
sorting, this algorithm tends to be less topologically accurate than \code{QDC} 
(which infers no metric information) or \code{WQDCrecursive} (which gives better topologies,
and reasonably accurate lengths for short edges, though long edge lengths may still be unreliable).
}
\examples{
gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
tnames=taxonNames(gtrees)
stree=WQDC(gtrees,tnames[1:6])
write.tree(stree)
plot(stree)

}
\references{
\insertRef{YR19}{MSCquartets}
}
\seealso{
\code{\link{quartetTable}}, 
         \code{\link{quartetTableResolved}}, 
         \code{\link{quartetTableDominant}}, 
         \code{\link{quartetWeightedDist}}, 
         \code{\link{WQDCrecursive}}, 
         \code{\link{WQDS}}, 
         \code{\link{QDC}}
}
