#' TINNIK algorithm to infer species tree of blobs
#'
#' Apply the TINNIK algorithm of \insertCite{ABMR24;textual}{MSCquartets} (see also \insertCite{ABMR22;textual}{MSCquartets})
#' to infer a tree of blobs for the species network from a collection of gene trees,
#' under the network multispecies coalescent (NMSC) model.
#'
#' @details  This function
#' \enumerate{
#' \item counts displayed quartets across gene trees to form quartet count concordance factors (qcCFs),
#' \item applies appropriate hypothesis tests to judge qcCFs as representing putative hybridization,
#' resolved trees, or unresolved (star) trees using \code{alpha} and \code{beta} as significance levels,
#' \item produces a simplex plot showing results of the hypothesis tests for all qcCFs
#' \item computes the appropriate TINNIK distance table, and infers the tree of blobs from the distance.
#' }
#'
#' A call of \code{TINNIK} with \code{genedata} given as a table previously output from \code{TINNIK} is
#' equivalent to a call of \code{TINNIKdist} followed by tree construction from the distance table.
#' If \code{genedata} is a
#' table previously output from \code{quartetTableResolved}
#' which lacks columns of p-values for hypothesis tests, these will be appended to the table output by \code{TINNIK}.
#' This table must contain a row with quartet counts for every 4 taxon set.
#'
#' If plots are produced, there are 2 simplex plots: The first shows the hypothesis test results,
#' and the second shows
#' inferred B-quartets and T-quartets. In both,
#' each point in the simplex plot corresponds to an empirical quartet concordance factor,
#' color-coded to represent test or inference results.
#'
#' In general, \code{alpha} should be chosen to be small and \code{beta}
#' to be large so that most quartets are interpreted as resolved trees. More quartets judges to have
#' either blob or unresolved relationships will lead to a less resolved blob tree.
#'
#' Usually, an initial call to \code{TINNIK} will not give a good analysis, as values
#' of \code{alpha} and \code{beta} are likely to need some adjustment based on inspecting the data. Saving the returned
#' table of test results from \code{TINNIK} will allow for the results of the time-consuming computation of qcCFs to be
#' saved, along with p-values,
#' for input to further calls of \code{TINNIK} with new choices of \code{alpha} and \code{beta}.
#'
#' See the documentation for \code{\link{TINNIKdist}} for an explanation of a small, rarely noticeable,
#' stochastic element of the algorithm.
#'
#' For data sets of many gene trees, user time may be reduced by using parallel code for
#' counting displayed quartets. See \code{\link{quartetTableParallel}}.
#'
#'
#' @references
#' \insertRef{ABMR22}{MSCquartets}
#'
#' \insertRef{ABMR24}{MSCquartets}
#'
#'
#' @param genedata gene tree data that may be supplied in any of 3 forms:
#' \enumerate{
#' \item as a character string giving the name of a file containing Newick gene trees,
#' \item as a multiPhylo object containing the gene trees, or
#' \item as a table of quartets on the gene trees, as produced by a previous call to
#' \code{TINNIK} or \code{quartetTableResolved}, which has columns only for taxa, resolved quartet counts,
#' and possibly p_T3, p_cut, and p_star
#' }
#' @param omit \code{FALSE} to treat unresolved quartets as 1/3 of each resolution;
#' \code{TRUE} to discard unresolved quartet data; ignored if gene tree data given as quartet table
#' @param epsilon minimum for branch lengths to be treated as non-zero; ignored if gene tree data given as quartet table
#' @param test a hypothesis test to perform, either "cut" or "T3" (default)
#' @param alpha a value or vector of significance levels for judging p-values for test specified by "test";
#' testing a null hypothesis of no hybridization vs. an alternative of hybridization, for each quartet;  a smaller value applies
#' a less conservative test for a tree (more trees), hence a stricter requirement for deciding in favor of hybridization (fewer reticulations)
#' @param beta a value or vector of significance levels for judging p-values testing
#' a null hypothesis of a star tree (polytomy) for each quartet vs. an alternative of anything else; a smaller value applies a less conservative
#' test for a star tree (more polytomies), hence a stricter requirement for deciding in favor of a resolved tree or network;
#' if vectors, \code{alpha} and \code{beta} must have the same length
#' @param treemethod a function implementing a method of tree inference from a distance table,
#' e.g. the ape package's fastme.bal or nj
#' @param delta a minimum edge length to retain in tree of blobs (see \insertCite{ABMR24}{MSCquartets} for related theory); shorter edges are collapsed
#' @param taxanames if \code{genedata} is a file or a multiPhylo object, a vector of a subset
#' of the taxa names on the gene trees
#' to be analyzed, if \code{NULL} all taxa on the first gene tree are used; if \code{genedata}
#' is a quartet table, this argument is ignored and all taxa in the table are used
#' @param plot \code{TRUE} produces simplex plots of hypothesis test results and plots the tree of blobs \code{FALSE} omits plots
#'
#' @return \code{output} (returned invisibly), with \code{output$ToB} the TINNIK tree of blobs, \code{output$pTable}
#' the table of quartets and p-values for judging fit to the MSC on quartet
#' trees, and \code{output$Bquartets} a TRUE/FALSE indicator vector of B-quartets; if \code{alpha, beta} are vectors, \code{output$ToB} is a vector of trees;
#' the table can be used as input to \code{TINNIK} or \code{TINNIKdist} with new choices of \code{alpha, beta}, without re-tallying quartets on
#' gene trees
#'
#' @seealso \code{\link{quartetTable}}, \code{\link{quartetTableParallel}}, \code{\link{quartetTableDominant}}, \code{\link{quartetCutTestInd}},\code{\link{quartetTreeTestInd}},
#' \code{\link{quartetStarTestInd}}, \code{\link{TINNIKdist}}, \code{\link{quartetTestPlot}}, \code{\link{pvalHist}}
#'
#' @examples
#' data(pTableYeastRokas)
#' out=TINNIK(pTableYeastRokas, alpha=.01, beta=.05)
#'
#' @importFrom ape di2multi nodelabels
#' @importFrom Rcpp evalCpp
#'
#' @export
TINNIK = function(genedata,
                  omit = FALSE,
                  epsilon = 0,
                  test = "T3",
                  alpha = .05,
                  beta = .95,
                  treemethod = fastme.bal,
                  delta = 2,
                  taxanames = NULL,
                  plot = TRUE) {
  if (!(test %in% c("cut", "T3"))) {
    stop("Argument 'test' must be 'cut' or 'T3'.")
  }
  pcol = paste0("p_", test)

  if (!(is.numeric(alpha) && is.numeric(beta))) {
    stop("Arguments 'alpha' and 'beta' must be numeric.")
  }
  if (length(alpha)!=length(beta)){
    stop("Arguments 'alpha' and 'beta' must have same number of entries.")
  }

  if ("matrix" %in% class(genedata)) {
    pTable = genedata
    momit=missing(omit)
    mepsilon=missing(epsilon)
    mtaxanames=missing(taxanames)
    if (!momit | !mepsilon | !mtaxanames){
      message(
        "Since genedata supplied as quartet table, ignoring arguments 'omit', 'epsilon', 'taxanames'."
      )}
    taxanames = setdiff(
      colnames(pTable),
      c(
        "12|34",
        "13|24",
        "14|23",
        "1234",
        "p_cut",
        "cutindex",
        "p_star",
        "p_T3",
        "p_T1"
      )
    )
  } else {
    if ("multiPhylo" %in% class(genedata))  {
      genetrees = genedata
    } else {
      if ("character" %in% class(genedata)) {
        genetrees <- read.tree(genedata) #read gene trees
        message(paste("Read", length(genetrees), "gene trees from file."))
      } else {
        stop("Data must be supplied as an object of type multiPhylo, character, or matrix.")
      }
    }

    if (is.null(taxanames)) {
      # if no taxa names specified,
      taxanames = genetrees[[1]]$tip.label   # ... get them from first tree
    }
    taxanames = sort(taxanames)
    if (length(taxanames) <= 25) {
      namelist = paste0(taxanames, collapse = ", ")
    } else {
      namelist = paste0(paste0(taxanames[1:25], collapse = ", "),
                        ",...(see output table for full list)")
    }
    message("Analyzing ", length(taxanames), " taxa: ", namelist)

    pTable = quartetTable(genetrees, taxanames, epsilon = epsilon)   # tally quartets on gene trees
    pTable = quartetTableResolved(pTable, omit)   # treat unresolved quartets
  }



  if (!(pcol %in% colnames(pTable))) {
    if (test == "cut") {
      pTable = quartetCutTestInd(pTable, lambda = 0) #compute p-values for cut test
    }
    if (test == "T3") {
      pTable = quartetTreeTestInd(pTable, model = "T3", lambda = 0)  #compute p-values for cut test
    }
  }



  if (!("p_star" %in% colnames(pTable))) {
    pTable = quartetStarTestInd(pTable)# compute p-values for each quartet CFs fit to star tree
  }

  ToBs=c() #empty list for ToBs

  for  (i in 1:length(alpha) ){
  ToBdist = TINNIKdist(pTable,
                       test = test,
                       alpha = alpha[i],
                       beta = beta[i]) # construct TINNIK distance table
  ToB = treemethod(ToBdist$dist) # infer tree of blobs
  ToB=collapseEdges(ToB,delta) #collapse short edges
  ToB = di2multi(ToB) # suppress edges of length 0

  #ToBdist$Bquartets encodes B-quartets T/F

  if (plot == TRUE) {
    if (test == "cut") {
      quartetTestPlot(pTable, "cut", alpha[i], beta[i]) # plot quartets in simplex, indicating which are interpretted as hybridization
      titletext = bquote("2-cut test," ~ alpha * "=" * .(alpha[i]) * "," ~ beta *
                           "=" * .(beta[i]))
    } else {
      if (test == "T3") {
        quartetTestPlot(pTable, "T3", alpha[i], beta[i]) # plot quartets in simplex, indicating which are interpreted as hybridization
        titletext = bquote("T3 test," ~ alpha * "=" * .(alpha[i]) * "," ~ beta *
                             "=" * .(beta[i]))
      }
    }


    quartetBTinferencePlot(pTable,ToBdist$Bquartets,test,alpha[i], beta[i])

    plot(ToB, type = "unrooted", main = "TINNIK Tree of Blobs")
    mtext(eval(bquote(.(titletext))),
          side = 3,
          line = 0,
          cex = 1) # add rest of title text

    ntaxa = length(taxanames)
    intnodes = (ntaxa + 1):(ntaxa + ToB$Nnode)

    nodelabels(
      node = intnodes,
      pch = 19,
      col = "red",
      cex = 1
    )
  }

  if (is.null(ToBs)) {
      ToBs=ToB
      } else {
    ToBs=c(ToBs,ToB) #append new ToB
      }
  }

  if (length(ToBs)==1) {ToBs=ToBs[[1]]} #if 1-element list, make just the element

  output = list(ToBs, pTable,ToBdist$Bquartets)
  names(output) = c("ToB", "pTable", "Bquartets")
  invisible(output)
}

###########################################################







# Cotangent function
#
#
#@param x argument
#@return function value
#
cot <- function(x) {
  1 / tan(x)
}
#########################################################

# Error function
#
#@param x argument
#@return function value
#
erf <- function(x) {
  2 * pnorm(x * sqrt(2)) - 1
}

#########################################################

# Inverse error function
#
#@param x argument
#@return function value
#
erf.inv <- function(x) {
  qnorm((x + 1) / 2) / sqrt(2)
}


###############################################################

# Log function allowing zero argument
#
# Log function, adjusted for use on zero argument: returns 0 when applied to 0
#
#@param x argument
#@return function value
#
logZ <- function(x) {
  output <- ifelse(x == 0, 0, log(x))
  return(output)
}
######################################################################

#' Probability density function for Cut Model
#'
#' Value of probability density function for Cut Model of \insertCite{ABMR24;textual}{MSCquartets}, Section 3.
#'
#' @references
#' \insertRef{ABMR24}{MSCquartets}
#'
#' @param lambda  statistic value (e.g., likelihood ratio statistic, or other power divergence statistic)
#' @param mu0 parameter of density function
#' @param alpha0 parameter of density function
#' @param beta0 parameter of density function
#'
#' @return  value of density function
#'
#' @seealso \code{\link{T1density},\link{T3density}}
#'
#' @export
cutDensity <- function(lambda, mu0, alpha0, beta0) {
  (1 / (2 * sqrt(2 * pi * lambda))) * (
    exp(-(1 / 2) * lambda) * (2 - erf((
      sqrt(lambda) * tan(alpha0 + beta0) + mu0
    ) / sqrt(2)) - erf((
      sqrt(lambda) * tan(alpha0 + beta0) - mu0
    ) / sqrt(2))) + exp(-(1 / 2) * (sqrt(lambda) + mu0 * cos(alpha0)) ^ 2) *
      (2 - erf((
        sqrt(lambda) * cot(alpha0) - mu0 * sin(alpha0)
      ) / sqrt(2)) - erf((
        sqrt(lambda) * cot(beta0) + mu0 * sin(alpha0)
      ) / sqrt(2))) + exp(-(1 / 2) * (sqrt(lambda) - mu0 * cos(alpha0)) ^ 2) *
      (2 - erf((
        sqrt(lambda) * cot(alpha0) + mu0 * sin(alpha0)
      ) / sqrt(2)) - erf((
        sqrt(lambda) * cot(beta0) - mu0 * sin(alpha0)
      ) / sqrt(2)))
  )
}



######################################################################

#' Maximum likelihood estimate of quartet tree of blobs topology and CF under Cut model
#'
#' Computes Maximum likelihood estimate of quartet tree of blobs topology and CF under the Cut model of
#' \insertCite{ABMR24;textual}{MSCquartets}, Section 3. In case of ties, the topology and CF estimate are chosen randomly among those maximizing
#' the likelihood. In particular, a resolved tree of blobs is always returned.
#'
#'@references
#'\insertRef{ABMR24}{MSCquartets}
#'
#'@param qcCF  a quartet count Concordance Factor (a 3-element vector)
#'
#'@return output with \code{output$topology} = 1, 2, or 3 indicating topology of
#'tree of blobs in accord with ordering of qcCF entries,
#' and \code{output$CFhat} the ML estimate for the CF
#'
#'@examples
#'  quartetCutMLE(c(17,72,11))
#'  quartetCutMLE(c(48,11,41))
#'  quartetCutMLE(c(11,48,41))
#'  quartetCutMLE(c(48,41,11))
#'
#'@export
#'
quartetCutMLE <- function(qcCF) {
  n = sum(qcCF)
  s1 = qcCF[2] + qcCF[3]
  L1 <- qcCF[1] * logZ(2 * qcCF[1]) + s1 * logZ(s1)
  s2 = qcCF[1] + qcCF[3]
  L2 <- qcCF[2] * logZ(2 * qcCF[2]) + s2 * logZ(s2)
  s3 = qcCF[1] + qcCF[2]
  L3 <- qcCF[3] * logZ(2 * qcCF[3]) + s3 * logZ(s3)
  L = c(L1, L2, L3)
  maxs = which(L == max(L))

  if (length(maxs) == 1) {
    top = maxs
  }
  else {
    top = sample(maxs, 1)
  }

  if (top == 1) {
    MLE <- c(qcCF[1], s1 / 2, s1 / 2) / n
  } else if (top == 2) {
    MLE <- c(s2 / 2, qcCF[2], s2 / 2) / n
  } else {
    MLE <- c(s3 / 2, s3 / 2, qcCF[3]) / n
  }
  output = list(top, MLE)
  names(output) = list("topology", "CFhat")
  return(output)
}


#############################################################


#' Hypothesis test for quartet counts fitting a resolved quartet tree of blobs under NMSC
#'
#' Test the hypothesis H_0=Cut model of \insertCite{ABMR24;textual}{MSCquartets},
#' Section 3., vs. H_1= everything else.
#' Returns p-value and estimate of tree of blobs topology.
#'
#' @details
#'
#' The Cut model for quartet CFs is the NMSC combined with the quartet species network having
#' a cut edge separating two of the taxa from the other two.
#'
#' This function implements the test described in \insertCite{ABMR24;textual}{MSCquartets}
#' as well as parametric bootstrapping,
#' with other procedures for when some expected counts are small.
#' These are more accurate tests than, say, a Chi-square with one degree of freedom,
#' which is not theoretically
#' justified near the singularity of the model, nor for small counts.
#'
#' If \code{method="MLtest"}, this uses the test for the Cut model described
#' in Section 3 of \insertCite{ABMR24;textual}{MSCquartets}, using the ML
#' estimate of the generating parameter.
#' As shown in simulations in that paper, the test is conservative when small
#' levels are used for rejection.
#' Although the test generally performs well in practice, it lacks a uniform
#' asymptotic guarantee over the full parameter space.
#'
#' If \code{method="conservative"}, the test uses the Chi-square distribution with 1 degree of freedom
#' (the "least favorable" approach). This is asymptotically guaranteed to reject the null
#' hypothesis at most at a specified level, but at the expense of increased type II errors.
#'
#' If \code{method="bootstrap"}, then parametric bootstrapping is performed,
#' based on ML estimates of the CF. The bootstrap sample size is given by the \code{bootstrap} argument.
#'
#' When some expected topology counts are small, the methods \code{"MLest"} and \code{"conservative"}
#' are not appropriate.
#' The argument \code{smallcounts} determines whether bootstrapping or a faster approximate method is used.
#' These use ML estimates of the CF under the Cut model.
#'
#' If two of the three counts are small (so the estimated CF is near a vertex of the simplex),
#' The approximate approach
#' returns a precomputed p-value, found by replacing the largest observed count
#' with 1e+6 and performing 1e+8 bootstraps. When n is sufficiently large (at least 30) and
#' some expected counts are small, the probability of topological error  is small and the bootstrap p-value is
#' approximately independent of the largest observed count.
#'
#' If one of the three counts is small (so the estimated CF is near an edge of the simplex),
#' a chi-squared test using the binomial model for the larger counts is used, as described
#' by \insertCite{ABMR24;textual}{MSCquartets}.
#'
#' The returned p-value should be taken with caution when there is a small sample size, e.g. less than 30 gene trees.
#'
#' @param obs  vector of 3 counts of resolved quartet frequencies
#' @param lambda  parameter for power-divergence statistic (e.g., 0 for likelihood ratio statistic, 1 for Chi-squared statistic)
#' @param method \code{"MLtest"},\code{"conservative"}, or \code{"bootstrap"}
#' @param smallcounts \code{"bootstrap"} or \code{"approximate"}, method of obtaining p-value when some counts are small
#' @param bootstraps  number of samples for bootstrapping
#'
#' @return output where \code{output$p.value} is a p-value and \code{output$topology} = 1, 2, or 3
#' indicates the ML estimate of the topology of the quartet tree of blobs in accord with ordering of qcCF entries.
#'
#'@examples
#'  quartetCutTest(c(17,72,11))
#'  quartetCutTest(c(48,11,41))
#'  quartetCutTest(c(11,48,41))
#'  quartetCutTest(c(48,41,11))
#'
#' @references
#' \insertRef{ABMR22}{MSCquartets}
#'
#' \insertRef{ABMR24}{MSCquartets}
#'
#' \insertRef{MAR19}{MSCquartets}
#'
#'
#' @seealso \code{\link{quartetCutTestInd}}
#'
#' @importFrom stats pbinom rmultinom integrate pnorm uniroot
#'
#' @export
#'
quartetCutTest <-
  function (obs,
            lambda = 0,
            method = "MLest",
            smallcounts = "approximate",
            bootstraps = 10 ^ 4) {
    if (is.numeric(obs) == FALSE || is.vector(obs) == FALSE ||
        length(obs) != 3 || min(obs) < 0 || sum(obs) < 1 || all.equal(sum(obs),
                                                                      round(sum(obs))) != TRUE || is.numeric(lambda) == FALSE ||
        is.vector(lambda) == FALSE || length(lambda) != 1 ||
        smallcounts %in% c("bootstrap", "approximate") ==
        FALSE || is.numeric(bootstraps) == FALSE || is.vector(bootstraps) ==
        FALSE || length(bootstraps) != 1 || bootstraps <= 0 ||
        bootstraps%%1 != 0 || method %in% c("MLest", "conservative",
                                            "bootstrap") == FALSE) {
      stop("Invalid arguments: obs must be a numeric non-negative vector of length 3 summing to a positive integer;\n      lambda must be a real number;\n      method must be \"MLest\", \"conservative\", or \"bootstrap\";\n      smallcounts must be \"bootstrap\" or \"approximate\";\n      bootstraps must be a positive integer.")
    }
    else {
      n <- sum(obs)
      if (n < 30) {
        warning("The number of gene quartets is <30; p-value may be inaccurate.")
      }
    }
    obs <- as.numeric(obs)
    cutMLE = quartetCutMLE(obs)
    top = cutMLE[[1]]
    expd <- n * (cutMLE[[2]])
    if (top == 2) {
      obs = obs[c(2, 1, 3)]
      expd = expd[c(2, 1, 3)]
    }
    else {
      if (top == 3) {
        obs = obs[c(3, 2, 1)]
        expd = expd[c(3, 2, 1)]
      }
    }
    if (all.equal(sum(abs(obs - expd)), 0) == TRUE) {
      p <- 1
    }
    else {
      stat <- powerDivStat(obs, expd, lambda)
      if (method == "bootstrap" || (method != "bootstrap" &&
                                    min(expd) < 5 && smallcounts == "bootstrap") ||
          (method != "bootstrap" && min(expd) < 5 &&
           expd[1] > expd[2] && smallcounts != "bootstrap" &&
           ((all.equal(sum(abs(obs * 3 - round(obs * 3))),
                       0) != TRUE) || (all.equal(sum(abs(obs * 3 -
                                                         round(obs * 3))), 0) == TRUE && var(obs - trunc(obs)) !=
                                       0)))) {
        message("Bootstrapping has been selected. p-values are approximate. Bootstrapping is only necessary when some expectations are small and approximate bootstrap p-values not available.")
        if (method == "bootstrap" && min(expd) < 5 &&
            smallcounts == "approximate") {
          message("Method is prioritized over smallcounts. Bootstrap has been chosen.")
        }
        if (method != "bootstrap" && min(expd) < 5 &&
            max(expd) > n - 10 && smallcounts != "bootstrap" &&
            ((all.equal(sum(abs(obs * 3 - round(obs * 3))),
                        0) != TRUE) || (all.equal(sum(abs(obs * 3 -
                                                          round(obs * 3))), 0) == TRUE && var(obs - trunc(obs)) !=
                                        0))) {
          warning("Approximate bootstrap p-values not available when two expected counts are below 5 and observed counts are not all integers, all integers + 1/3 or all integers + 2/3.")
        }
        if (bootstraps == 0) {
          bootstraps <- 10^4
        }
        count <- 0
        sims <- rmultinom(bootstraps, n, prob = expd/n)
        for (i in 1:bootstraps) {
          dat <- sims[, i]
          expected <- n * quartetCutMLE(dat)[[2]]
          simstat <- powerDivStat(dat, expected, lambda)
          count <- ifelse(simstat >= stat, count + 1, count)
        }
        p <- count/bootstraps
      }
      else if (method != "bootstrap" && min(expd) < 5 &&
               expd[1] >= expd[2] && smallcounts == "approximate" &&
               (all.equal(sum(abs(obs * 3 - round(obs * 3))), 0) ==
                TRUE) && var(obs - trunc(obs)) == 0) {
        sortobs23times3 = sort(round(3 * obs[2:3]))
        if (sortobs23times3[1]%%3 == 0) {
          if (identical(sortobs23times3, c(0, 27))) {
            p <- 0.000919
          }
          else if (identical(sortobs23times3, c(3, 24))) {
            p <- 0.0209
          }
          else if (identical(sortobs23times3, c(6, 21))) {
            p <- 0.103
          }
          else if (identical(sortobs23times3, c(9, 18))) {
            p <- 0.359
          }
          else if (identical(sortobs23times3, c(12, 15))) {
            p <- 0.79
          }
          else if (identical(sortobs23times3, c(0, 24))) {
            p <- 0.00191
          }
          else if (identical(sortobs23times3, c(3, 21))) {
            p <- 0.0397
          }
          else if (identical(sortobs23times3, c(6, 18))) {
            p <- 0.175
          }
          else if (identical(sortobs23times3, c(9, 15))) {
            p <- 0.523
          }
          else if (identical(sortobs23times3, c(0, 21))) {
            p <- 0.00416
          }
          else if (identical(sortobs23times3, c(3, 18))) {
            p <- 0.0763
          }
          else if (identical(sortobs23times3, c(6, 15))) {
            p <- 0.297
          }
          else if (identical(sortobs23times3, c(9, 12))) {
            p <- 0.768
          }
          else if (identical(sortobs23times3, c(0, 18))) {
            p <- 0.00871
          }
          else if (identical(sortobs23times3, c(3, 15))) {
            p <- 0.129
          }
          else if (identical(sortobs23times3, c(6, 12))) {
            p <- 0.476
          }
          else if (identical(sortobs23times3, c(0, 15))) {
            p <- 0.0183
          }
          else if (identical(sortobs23times3, c(3, 12))) {
            p <- 0.24
          }
          else if (identical(sortobs23times3, c(6, 9))) {
            p <- 0.737
          }
          else if (identical(sortobs23times3, c(0, 12))) {
            p <- 0.0393
          }
          else if (identical(sortobs23times3, c(3, 9))) {
            p <- 0.439
          }
          else if (identical(sortobs23times3, c(0, 9))) {
            p <- 0.086
          }
          else if (identical(sortobs23times3, c(3, 6))) {
            p <- 0.681
          }
          else if (identical(sortobs23times3, c(0, 6))) {
            p <- 0.197
          }
          else if (identical(sortobs23times3, c(0, 3))) {
            p <- 0.478
          }
        }
        else if (sortobs23times3[1]%%3 == 1) {
          if (identical(sortobs23times3, c(1, 28))) {
            p <- 0.00222
          }
          else if (identical(sortobs23times3, c(4, 25))) {
            p <- 0.0237
          }
          else if (identical(sortobs23times3, c(7, 22))) {
            p <- 0.119
          }
          else if (identical(sortobs23times3, c(10, 19))) {
            p <- 0.358
          }
          else if (identical(sortobs23times3, c(13, 16))) {
            p <- 0.777
          }
          else if (identical(sortobs23times3, c(1, 25))) {
            p <- 0.00452
          }
          else if (identical(sortobs23times3, c(4, 22))) {
            p <- 0.0436
          }
          else if (identical(sortobs23times3, c(7, 19))) {
            p <- 0.203
          }
          else if (identical(sortobs23times3, c(10, 16))) {
            p <- 0.52
          }
          else if (identical(sortobs23times3, c(1, 22))) {
            p <- 0.0094
          }
          else if (identical(sortobs23times3, c(4, 19))) {
            p <- 0.0794
          }
          else if (identical(sortobs23times3, c(7, 16))) {
            p <- 0.295
          }
          else if (identical(sortobs23times3, c(10, 13))) {
            p <- 0.754
          }
          else if (identical(sortobs23times3, c(1, 19))) {
            p <- 0.0194
          }
          else if (identical(sortobs23times3, c(4, 16))) {
            p <- 0.143
          }
          else if (identical(sortobs23times3, c(7, 13))) {
            p <- 0.47
          }
          else if (identical(sortobs23times3, c(1, 16))) {
            p <- 0.0404
          }
          else if (identical(sortobs23times3, c(4, 13))) {
            p <- 0.233
          }
          else if (identical(sortobs23times3, c(7, 10))) {
            p <- 0.72
          }
          else if (identical(sortobs23times3, c(1, 13))) {
            p <- 0.085
          }
          else if (identical(sortobs23times3, c(4, 10))) {
            p <- 0.422
          }
          else if (identical(sortobs23times3, c(1, 10))) {
            p <- 0.113
          }
          else if (identical(sortobs23times3, c(4, 7))) {
            p <- 0.665
          }
          else if (identical(sortobs23times3, c(1, 7))) {
            p <- 0.237
          }
          else if (identical(sortobs23times3, c(1, 4))) {
            p <- 0.532
          }
        }
        else {
          if (identical(sortobs23times3, c(2, 26))) {
            p <- 0.00827
          }
          else if (identical(sortobs23times3, c(5, 23))) {
            p <- 0.0436
          }
          else if (identical(sortobs23times3, c(8, 20))) {
            p <- 0.199
          }
          else if (identical(sortobs23times3, c(11, 17))) {
            p <- 0.516
          }
          else if (identical(sortobs23times3, c(2, 23))) {
            p <- 0.0164
          }
          else if (identical(sortobs23times3, c(5, 20))) {
            p <- 0.0777
          }
          else if (identical(sortobs23times3, c(8, 17))) {
            p <- 0.298
          }
          else if (identical(sortobs23times3, c(11, 14))) {
            p <- 0.739
          }
          else if (identical(sortobs23times3, c(2, 20))) {
            p <- 0.0237
          }
          else if (identical(sortobs23times3, c(5, 17))) {
            p <- 0.146
          }
          else if (identical(sortobs23times3, c(8, 14))) {
            p <- 0.466
          }
          else if (identical(sortobs23times3, c(2, 17))) {
            p <- 0.0468
          }
          else if (identical(sortobs23times3, c(5, 14))) {
            p <- 0.238
          }
          else if (identical(sortobs23times3, c(8, 11))) {
            p <- 0.703
          }
          else if (identical(sortobs23times3, c(2, 14))) {
            p <- 0.0926
          }
          else if (identical(sortobs23times3, c(5, 11))) {
            p <- 0.409
          }
          else if (identical(sortobs23times3, c(2, 11))) {
            p <- 0.185
          }
          else if (identical(sortobs23times3, c(5, 8))) {
            p <- 0.645
          }
          else if (identical(sortobs23times3, c(2, 8))) {
            p <- 0.377
          }
          else if (identical(sortobs23times3, c(2, 5))) {
            p <- 0.525
          }
        }
      }
      else if (method != "bootstrap" && min(expd) < 5 &&
               expd[1] < expd[2] && smallcounts == "approximate" &&
               all.equal(sum(abs(obs - round(obs))), 0) == TRUE) {
        p <- 2 * pbinom(min(obs[2:3]), sum(obs[2:3]), 1/2)
      }
      else if (method != "bootstrap" && min(expd) < 5 &&
               expd[1] < expd[2] && smallcounts == "approximate" &&
               all.equal(sum(abs(obs - round(obs))), 0) != TRUE) {
        p <- pchisq(stat, 1, lower.tail = FALSE)
      }
      else if (method != "bootstrap" && min(expd) >=
               5) {
        if (method == "MLest") {
          phinhat <- 3 * (n - expd[1])/(2 * n)
          munhat <- sqrt(2 * n) * (1 - phinhat)/sqrt(phinhat *
                                                       (3 - 2 * phinhat))
          alphanhat <- atan(1/sqrt(3 * (3 - 2 * phinhat)))
          betanhat <- (1/2) * (pi/2 - phinhat)
          p <- ifelse(stat == 0, 1, min(1, max(0, integrate(cutDensity,
                                                            lower = stat, upper = Inf, mu0 = munhat, alpha0 = alphanhat,
                                                            beta0 = betanhat)$value)))
        }
        else if (method == "conservative") {
          p <- pchisq(stat, 1, lower.tail = FALSE)
        }
      }
    }
    ptop = list(p, top)
    names(ptop) = c("p.value", "topology")
    return(ptop)
  }

#############################################################

#' Multiple independent hypothesis tests for quartet counts fitting the Cut model under the NMSC
#'
#' Perform a hypothesis test for all quartet counts in an input table, as if the counts for different choices of 4 taxa
#' are independent.
#'
#' @details This function assumes all quartets are resolved.  The test performed and the arguments
#' are described more fully in \code{quartetCutTest}.
#'
#' @references
#' \insertRef{ABMR24}{MSCquartets}
#'
#' @param rqt  table of resolved quartet counts, as produced by \code{quartetTableResolved}, or \code{quartetStarTestInd}
#' @param lambda power divergence statistic parameter (e.g., 0 for likelihood ratio statistic, 1 for Chi-squared statistic)
#' @param method \code{"MLest"}, \code{"conservative"}, or \code{"bootstrap"}; see \code{quartetCutTest} for explanation
#' @param smallcounts \code{"bootstrap"} or \code{"approximate"}, method of obtaining p-value when some counts are small, so
#' the chosen \code{method} is inappropriate
#' @param bootstraps  number of samples for bootstrapping
#'
#' @return
#'    a copy of \code{rqt} with two columns appended: \code{"p_cut"} with p-values and \code{"cutindex"}
#'    giving index 1,2, or 3 of ML estimate of quartet tree of blobs (1 if 12|34, 2 if 13|24, 3 if 14|23)
#'    under Cut model.
#'
#' @seealso \code{\link{quartetCutTest}}, \code{\link{quartetTestPlot}}, \code{\link{quartetStarTestInd}}, \code{\link{quartetTableResolved}}
#'
#' @examples
#' data(pTableYeastRokas)
#' pT=pTableYeastRokas[1:10,1:11]
#' pTable=quartetCutTestInd(pT)
#'
#' @importFrom ape unroot read.tree cophenetic.phylo
#' @export
quartetCutTestInd <-
  function(rqt,
           lambda = 0,
           method = "MLest",
           smallcounts = "approximate",
           bootstraps = 10 ^ 4) {
    if (is.vector(lambda) == FALSE ||
        length(lambda) != 1 ||
        smallcounts %in% c("bootstrap", "approximate") == FALSE ||
        bootstraps <= 0 || bootstraps %% 1 != 0 ||
        method %in% c("MLest", "conservative", "bootstrap") == FALSE) {
      # Return an error and stop if arguments are misspecified.
      stop(
        'Invalid arguments: lambda must be real; smallcounts must be "bootstrap" or "approximate";
      bootstraps must be a positive integer; method must be "MLest","conservative", or "bootstrap".'
      )
    }
    colnames = colnames(rqt)
    if (("p_cut" %in% colnames) | ("cutindex" %in% colnames)) {
      stop('Input table already has column to be appended.')
    }
    taxanames = setdiff(
      colnames,
      c(
        "12|34",
        "13|24",
        "14|23",
        "1234",
        "p_cut",
        "cutindex",
        "p_star",
        "p_T3",
        "p_T1"
      )
    )

    M = dim(rqt)[1] # number of quartets
    n = length(taxanames) # number of taxa

    pTable = cbind(rqt, p_cut = 0, cutindex = 0) # tack on columns for p-values, and cutindices
    message("Applying hypothesis test for Cut model to ", M, " quartets.")
    for (m in 1:M) {
      # consider each quartet
      obs = unname(rqt[m, c("12|34", "13|24", "14|23")])
      ptop = quartetCutTest(
        obs,
        lambda = lambda,
        smallcounts = smallcounts,
        bootstraps = bootstraps,
        method = method
      )  #  compute p-value
      pTable[m, "p_cut"] = ptop$p.value  # store p-value
      pTable[m, "cutindex"] = ptop$topology # and topology
    }
    return(pTable)

  }


#################################################################
#' Compute TINNIK distance from quartets and hypothesis test p-values
#'
#' Apply the B-quartet inference algorithm of \insertCite{ABMR22;textual}{MSCquartets}, \insertCite{ABMR24;textual}{MSCquartets} to
#' infer all B-quartets from results of hypothesis tests, and then compute an estimate of an intertaxon distance
#' fitting the topological tree of blobs of the species network.
#'
#' @details This function assumes \code{pTable} has columns for taxa and resolved
#' quartet counts as originally produced by \code{quartetTable},
#' and hypothesis test results as produced by
#' \code{quartetStarTestInd}, and either \code{quartetTreeTestInd} for the \code{T3} test or \code{quartetCutTestInd}.
#' Rows must be present for every 4-taxon subset.
#' (Note: Of functions in this package, only \code{HolmBonferroni} might modify the row order from the required one.)
#'
#' This function uses the Rcpp package for significant speed up in computation time.
#'
#' @references
#' \insertRef{ABMR22}{MSCquartets}
#'
#' \insertRef{ABMR24}{MSCquartets}
#'
#' @param pTable  table of resolved quartet counts, as produced by
#' \code{quartetTableResolved}, with extra columns from
#' both star hypothesis test, and either cut or T3 hypothesis tests
#' @param test either "cut" or "T3"
#' @param alpha level for cut or T3 test
#' @param beta level for star test
#'
#' @return a distance table \code{output$dist} and
#'    a vector \code{output$Bquartets} with TRUE/FALSE entries indicating B-quartets
#'    ordered as rows of \code{pTable}.
#'
#' @seealso \code{\link{quartetTable},\link{quartetTableResolved},\link{quartetStarTest}},
#' \code{\link{quartetCutTest}}, \code{\link{quartetStarTestInd}}, \code{\link{quartetCutTestInd}}
#'
#' @examples
#' data(pTableYeastRokas)
#' out=TINNIKdist(pTableYeastRokas,test="T3",alpha=.05,beta=.05)
#'
#' @export
TINNIKdist <- function(pTable,
                       test = "T3",
                       alpha = .05,
                       beta = .05) {
  if (!(test %in% c("cut", "T3")) |
      ((test == "cut") &
       !all(c("p_cut", "cutindex", "p_star") %in% colnames(pTable))) |
      ((test == "T3") &
       !all(c("p_T3", "p_star") %in% colnames(pTable)))) {
    stop('Invalid test argument, or input table missing columns from tests.')
  }


  taxanames = setdiff(
    colnames(pTable),
    c(
      "12|34",
      "13|24",
      "14|23",
      "1234",
      "p_cut",
      "cutindex",
      "p_star",
      "p_T3",
      "p_T1"
    )
  )
  n = length(taxanames)# number of taxa
  m = dim(pTable)[1]# number of 4-taxon sets in table
  if (m != choose(n, 4)) {
    stop("Dimensions of pTable not correct; must have row for each choice of 4 taxa.")# check number of rows is correct
  }

  L1 = c() #empty vectors for storing indices (i.e., row numbers) of B-quartets in pTable
  L2 = L1

  if (test == "cut") {
    cuttops = pTable[, "cutindex"] #already have topology for MLE for tree of blobs
  } else {
    #here if test="T3"
    cuttops = max.col(pTable[, c("12|34", "13|24", "14|23")]) # Use dominant topology (MLE) appropriate for T3 model
  }

  # Initialization of B-quartets from hypothesis tests
  Bquartets = rep(FALSE, m) #create initial indicator and index vectors for B-quartets inferred from CFs
  p_test = paste0("p_", test) #create p-value column names for test
  colptest= which(colnames(pTable) ==p_test)
  colpstar= which(colnames(pTable) =="p_star")
  ret=initBquartets(pTable,m,alpha,beta,colptest,colpstar,Bquartets)
  L1=(ret$Lout)[(ret$Lout)!=0]
  Bquartets=as.logical(ret$Bout)

  lenL1 = length(L1) # number of new B-quartets found
  Binit = lenL1  # save initial number of B-quartets

  if (Binit == 0) {
    message(
      c(
        "No support for hybridizations/gene flow in any quartet using ",
        test,
        " test with  alpha=",
        alpha,
        ", beta=",
        beta,
        "."
      )
    )
  }

  np1 = n + 1 #number of taxa + 1
  C = matrix(nrow = n, ncol = 4) #cache binomial coefficients for efficiency in computing indices
  for (i in 0:(n - 1)) {
    for (j in 1:4) {
      C[i + 1, j] = choose(i, j)
    }
  }
  Cn4 = choose(n, 4) # need one extra value

  Nrule1 = 0 #counters of uses of inference rules
  Nrule2 = 0

  # Perform B-quartet inference procedure using C++ code
  Bquartets=BQinference(pTable,C,Cn4,n,Bquartets,L1,lenL1,Nrule1,Nrule2,cuttops)

  #After all B-quartets are inferred, compute TINNiK distance
  taxa = colnames(pTable)[1:n]# get names of taxa
  D = matrix(0, n, n)# create blank distance matrix
  colnames(D) = taxa
  rownames(D) = taxa

  for (i in 1:m) {
    taxa = which(pTable[i, 1:n] == 1)
    #if not a B-quartet
    if (Bquartets[i] == FALSE) {
      ones = taxa[c(1, cuttops[i] + 1)]
      nones = setdiff(taxa, ones)
      D[ones, nones] = D[ones, nones] + 1 # increment 4 entries of D
      D[nones, ones] = D[nones, ones] + 1 # and another 4
    } else {
      #if a B-quartet
      D[taxa, taxa] = D[taxa, taxa] + 1 #increment 16 entries of D (diagonal will be reset to 0 later) to treat as unresolved
    }
  }

  D = 2 * D + 2 * n - 4 # adjust to appropriate distance formula
  for (i in 1:n) {
    D[i, i] = 0 # and set diagonal back to 0
  }

  #print(paste("B-quartets from initial hypothesis tests:", Binit))
  #print(paste("B-quartets from inference rules:", sum(Bquartets) - Binit))
  #print(paste("Rule 1 used", Nrule1, "times."))
  #print(paste("Rule 2 used", Nrule2, "times."))

  output = list(dist = D, Bquartets = Bquartets) # return B-quartet indicators and distance
  return(output)
}



##############################################

#' Groups taxa by deleting a node in a tree
#'
#' Finds groups of taxa determined by the connected components of the graph resulting from deleting an internal node in a tree.
#'
#' @details
#' When applied to a rooted tree, the last group returned is the set of tips that are
#' non-descendants of the node (provided any exist).
#'
#'@param tree  a tree, of class "phylo"
#'@param nodeNum a node number, representing an internal node in the phylo representation
#'
#'@return a list of lists of tree tip numbers for each group.
#'The union of the groups is the set of all tips.
#'
#'@examples
#'  tree=read.tree(text="((a,b),((c,d,e),(f,g)));")
#'  nodeGroups(tree,8)
#'  nodeGroups(tree,10)
#'  nodeGroups(tree,11)
#'
#'@export
#'
nodeGroups <- function(tree, nodeNum) {
  ntaxa = length(tree$tip.label)

  if ((nodeNum <= ntaxa) || (nodeNum > (ntaxa + tree$Nnode))) {
    stop("Argument 'nodeNum' not in valid range.")
  }

  children = Descendants(tree, nodeNum, "children") #find children of node
  groups = list() #begin list
  alld = c() #and set of all descendents
  for (c in children) {
    d = Descendants(tree, c, 'tips')
    groups = append(groups, d)
    alld = union(alld, d[[1]])
  }

  nonDesc = setdiff(1:ntaxa, alld)

  if (length(nonDesc) > 0) {
    groups = append(groups, list(nonDesc))
  }

  return(groups)
}



################################

#' Produce simplex plot with  results of B/T-quartet inference
#'
#' Plot a 2-d probability simplex, with points for all normalized quartet count
#' vectors. Colors indicate B- or T-quartets from TINNIK algorithm, at specified
#' test levels.
#'
#' @details The first argument of this function is a table of quartets and p-values. The
#' plot may show results using either the T3, or 2-cut
#' test, and a star tree test.
#' The p-values must be computed by or before previous calls to
#' \code{TINNIK}. The second argument is the indicator vector for B/T quartets produced by \code{TINNIK}.
#'
#' @param pTable table of quartets and p-values
#' @param Bquartets indicator vector for B-quartets (1=B, 0=T), ordered as in pTable
#' @param test  test model used for tree null hypothesis; options are \code{"cut"}, \code{"T3"}
#' @param alpha  significance level used by TINNIK for test \code{test}
#' @param beta  significance level used by TINNIK for star tree test
#' @param cex scaling factor for size of plotted symbols
#'
#' @return NULL
#'
#' @seealso \code{TINNIK}, \code{\link{quartetTestPlot}}
#'
#' @examples
#' data(pTableYeastRokas)
#' out=TINNIKdist(pTableYeastRokas,test="T3",alpha=.05,beta=.05)
#' quartetBTinferencePlot(pTableYeastRokas,out$B,test="T3",alpha=.05,beta=.05)
#'
#' @export
quartetBTinferencePlot <- function(pTable,
                            Bquartets,
                            test,
                            alpha,
                            beta,
                            cex=1) {
  if (!(is.numeric(alpha) && is.numeric(beta))) {
    stop("Test levels alpha and beta must be numeric.")
  }

  if ( !(test %in% c("T3","cut") ) ){ stop("Invalid test name.")}

  lineWidth = 1.5
  redColor = "goldenrod"
  blueColor = "darkgreen"

  counts = c("12|34", "13|24", "14|23")
  M = dim(pTable)[1]


  if (test == "T3")  {
    #T3 and star
    titletext = bquote("T3," ~ alpha * "=" * .(alpha) * "," ~ beta *
                         "=" * .(beta))
  } else {
      #cut and star
      titletext = bquote("2-cut," ~ alpha * "=" * .(alpha) * "," ~ beta *
                           "=" * .(beta))
  }

    legtext = c("B-quartet",
                "T-quartet")
    legpch = c(5, 1)
    legcol = c(redColor, blueColor)
    model = test

simplexPrepare(model, maintitle = "TINNIK Inferred B-,T-quartets ", titletext =
                   titletext) # draw outline of simplex and lines for model

  reds = Bquartets
  blues = !Bquartets

  bluepoints = pTable[which(blues), counts, drop = FALSE] # points to plot blue
  nblue = dim(bluepoints)[1]
  redpoints = pTable[which(reds), counts, drop = FALSE] #      red
  nred = dim(redpoints)[1]

  if (nblue > 0) {
    for (i in 1:nblue) {
      simplexPoint(
        bluepoints[i,],
        type = "o",
        pch = 1,
        col = blueColor,
        lwd = lineWidth,
        cex = cex
      )
    }
  }
  if (nred > 0) {
    for (i in 1:nred) {
      simplexPoint(
        redpoints[i,],
        type = "o",
        pch = 5,
        col = redColor,
        lwd = lineWidth,
        cex = cex
      )
    }
  }

  legend(
    'topleft',
    legend = legtext,
    pch = legpch,
    col = legcol,
    cex = 0.8,
    inset = 0.01,
    pt.cex = 1.2,
    lty = 0,
    lwd = 2
  )
}

##################################################################
#' Tree of blobs for a network
#'
#' Given extended newick, an evonet object, or an igraph object for a network, return its reduced, unrooted
#' tree of blobs. Here 'reduced' means all nodes resulting from 2-blobs are suppressed, as are edges above
#' the network's LSA.
#'
#' @param net A network, supplied as an extended Newick string, an evonet object, or an igraph object
#' @param plot if TRUE (default), plot the tree of blobs.
#'
#'
#'@return An object of class \code{phylo} containing the unrooted topological tree
#' derived from the network by contracting all blobs.  All edge lengths are 1.
#'
#'@seealso \code{\link{TINNIK}}
#'
#'@examples network = "(((a:1,d:1):1,(b:1)#H1:1):1,(#H1,c:1):2);"
#' plot(read.evonet(text=network))
#' treeOfBlobs(network, plot=TRUE)
#'
#'@importFrom igraph as.igraph is_directed degree V E all_simple_paths difference intersection
#'@importFrom igraph edge_attr distances
#'@importFrom ape read.evonet .PlotPhyloEnv
#'@importFrom utils combn
#'
#' @export
treeOfBlobs = function(net, plot = FALSE) {
  # convert network to igraph
  if (is(net, "character")) {
    if (length(unlist(gregexpr('#', net)))<2  ) stop("Invalid network argument") # check for a hybrid node
    net = read.evonet(text = net)
    net = igraph::as.igraph(net)
  } else {
    if (is(net, "evonet")) {
      net = igraph::as.igraph(net)
    } else {
      if (!is(net, "igraph")) {
        stop("Network argument must be a newick string, evonet, igraph.")
      }
    }
  }

  if (!(igraph::is_directed(net))) {
    stop("Network must be rooted.")
  }

  leaves = which(igraph::degree(net, v = igraph::V(net), mode = "out")==0, useNames = TRUE)
  root = which(igraph::degree(net, v = igraph::V(net), mode = "in")==0, useNames = TRUE)

  blob_edges = c()

  # obtain blob edges, by getting all paths from the root to taxa
  for (cur_leaf in 1:length(leaves)) {
    S = igraph::all_simple_paths(net, from = root, to = leaves[cur_leaf])
    if (length(S)>1) {

      # get number of pairs of paths in S
      nc2 = combn(1:length(S), 2, simplify = FALSE)
      for (j in 1:length(nc2)) {
        v_path1 = S[[nc2[[j]][1]]]
        v_path2 = S[[nc2[[j]][2]]]
        path1 = igraph::E(net,path=v_path1)
        path2 = igraph::E(net,path=v_path2)
        common_edges = igraph::intersection(path1,path2)
        blob_edges = unique(c(igraph::difference(path1,common_edges),igraph::difference(path2,common_edges),blob_edges))
      }
    }
  }

  igraph::edge_attr(net,"weight") = 1
  igraph::edge_attr(net,"weight",index = blob_edges) <- 0

  dTable = igraph::distances(net,leaves,leaves)
  tob = nj(dTable) #construct tree
  tob = di2multi(unroot(tob), tol = 1e-10) #unroot and suppress degree 2 nodes
  tob$edge.length=rep(1,length(tob$edge.length)) # again make all edges length 1

  # plot tob
  if (plot) {
    plot(tob,main="Tree of blobs",type="unrooted",cex=0.7,lab4ut="axial")
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
    XX <- lastPP$xx[internal_nodes]
    YY <- lastPP$yy[internal_nodes]
    points(XX,YY,pch=16,col="blue",cex=1)
  }
  return(tob)
}

##################################################################
#' Extract compatible splits
#'
#' Given an object of class splits, first discards any with weight less than a tolerance, and
#' then further removes all remaining splits that are incompatible with any other remaining one.
#
#'
#'@param sp an object of class splits
#'@param tol splits with weights below tol are dropped
#'@param plot a logical; if TRUE plots tree displaying remaining spilts
#'
#'@return splits objects containing only those that are compatible and high weight
#'
#'@examples
#' data(pTableYeastRokas)
#' dist=NANUQdist(pTableYeastRokas, alpha=.05, beta=.95,outfile=NULL)
#' nn=neighborNet(dist)
#' plot(nn,"2D")
#' tob=treeFromSplits(compatibleSplits(nn$splits),plot=TRUE) #produce tree of blobs of splits graph
#'
#'@seealso \code{\link{treeFromSplits}, \link{TINNIK}}
#'
#'@export
compatibleSplits = function(sp, tol = 0, plot = FALSE) {
  if (!is(sp,"splits")) {
    stop("First argument must be of class 'splits'.")
  }

  # compute number of splits
  num_splits <- length(attr(sp,"weights"))

  # filter low support splits using value of tol by setting
  # any split of weight < tol to have weight 0
  low_support_splits <- which(attr(sp,"weights") < tol)
  high_support_splits <- setdiff(1:num_splits,low_support_splits)
  sp <- sp[high_support_splits]

  # determine which remaining splits are compatible
  cmat = as.matrix(compatible(sp))

  # get indices of incompatible splits
  nrows = dim(cmat)[1]
  incompatible_splits = c()
  for (rnum in 1:nrows) {
    if (any(cmat[rnum,]==1)) {
      incompatible_splits <- c(incompatible_splits,rnum)
    }
  }

  # eliminate incompatible splits
  num_splits <- length(attr(sp,"weights"))
  comp_splits <- setdiff(1:num_splits,incompatible_splits)
  sp <- sp[comp_splits]

  if (plot) {
    plot(as.networx(sp))
  }

  # return the splits system with only the compatible splits
  invisible(sp)
}

################################################
#'Produce tree from compatible splits
#'
#'@param sp a compatible split system, as produced by compatibleSplits
#'@param plot, a logical, if TRUE, plot tree
#'
#'@return a phylo object for tree displaying splits
#'
#'@examples
#' data(pTableYeastRokas)
#' dist=NANUQdist(pTableYeastRokas, alpha=.05, beta=.95,outfile=NULL)
#' nn=neighborNet(dist)
#' plot(nn,"2D")
#' tob=treeFromSplits(compatibleSplits(nn$splits),plot=TRUE) #produce tree of blobs of splits graph
#'
#' @seealso \code{\link{compatibleSplits}, \link{TINNIK}}
#'
#'@importFrom phangorn compatible as.networx
#'@importFrom ape nj
#'@importFrom stats cophenetic
#'
#'@export
treeFromSplits = function(sp, plot = FALSE) {

  if (length(compatibleSplits(sp))<length(sp) ) {
    stop("Splits are not all compatible.")
  }

  nn = as.networx(sp)

  D = cophenetic(nn)
  tree = nj(D)
  tree = di2multi(unroot(tree), tol = 1e-10)

  if (plot) {
    plot(tree,type = "unrooted",
         main = "Tree from compatible splits system")
  }
  return(tree)
}

############################################################
#' Collapse short tree edges
#'
#' Set all edges of a tree with short lengths to be zero.
#'
#'@param tree a phylo object
#'@param delta minimum edge length to retain
#'
#'@return a phylo object
#'
#'@examples
#' tree=read.tree(text="((a:1,b:1):.5,(c:.5,d:2):1);")
#' newtree=collapseEdges(tree,delta=1)
#' write.tree(newtree)
#'@export

collapseEdges=function(tree,delta) {
  ntree=tree
  ntree$edge.length=ntree$edge.length*(ntree$edge.length>=delta)
  ntree = di2multi(ntree) # suppress edges of length 0
  return(ntree)
}
