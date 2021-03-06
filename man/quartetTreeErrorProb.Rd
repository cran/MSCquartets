% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TreeQuartetHypTest.R
\name{quartetTreeErrorProb}
\alias{quartetTreeErrorProb}
\title{Bayesian posterior probability of error in 4-taxon unrooted species tree topology estimate}
\usage{
quartetTreeErrorProb(obs, model = "T3")
}
\arguments{
\item{obs}{vector of counts for 3 topologies}

\item{model}{\code{"T3"} or \code{"T1"}, for the models of \insertCite{MAR19;textual}{MSCquartets} describing an unspecified species
tree topology (\code{"T3"}), or the topology whose count is the first entry of \code{obs} (\code{"T1"})}
}
\value{
\code{(error.prob, top.probs)} where \code{error.prob} is the species tree error probability
and \code{top.probs} is a vector of the three species tree topology probabilities in the order of \code{obs};
for model \code{"T1"} the species tree used is the one
corresponding to the first count; for model \code{"T3"} the species
tree is the one corresponding to the largest count
}
\description{
From a gene quartet count concordance factor (qcCF), computes Bayesian posterior probabilities
of the three 4-taxon species tree topologies and the Bayesian posterior probability that the
assumed topology is incorrect, under the assumption that the counts arise from the MSC on some species tree.
}
\details{
The Jeffreys prior is used for internal branch length, along with the uniform prior
on the resolved topology.
}
\examples{
obs <- c(28,32,30)
quartetTreeErrorProb(obs,model="T1")
quartetTreeErrorProb(obs,model="T3")

}
\references{
\insertRef{MAR19}{MSCquartets}
}
