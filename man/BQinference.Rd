% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{BQinference}
\alias{BQinference}
\title{Main loop of B-quartet inference}
\usage{
BQinference(pTable, C, Cn4, n, Bquartets, L1, lenL1, Nrule1, Nrule2, cuttops)
}
\arguments{
\item{pTable}{a quartet table with p-values}

\item{C}{precomputed binomial coefficients}

\item{Cn4}{precomputed binomial coefficient}

\item{n}{number of taxa}

\item{Bquartets}{0/1 vector of initial Bquartets}

\item{L1}{vector of recently inferred B quartets}

\item{lenL1}{lnegth of L1}

\item{Nrule1}{count of inference from rule 1}

\item{Nrule2}{count of inference from rule 2}

\item{cuttops}{inferred cut topologies}
}
\description{
This is a C++ function, called from \code{TINNIKdist}, to
infer B and T quartets.
}
\seealso{
\code{\link{quartetTable}}, \code{\link{quartetTableParallel}}
}
