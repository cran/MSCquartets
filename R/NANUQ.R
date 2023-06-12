#' Apply NANUQ network inference algorithm to gene tree data
#'
#' Apply the NANUQ algorithm of \insertCite{ABR19;textual}{MSCquartets} to infer a hybridization network from a collection of gene trees,
#' under the level-1 network multispecies coalescent (NMSC) model.
#'
#' @details  This function
#' \enumerate{
#' \item counts displayed quartets across gene trees to form quartet count concordance factors (qcCFs),
#' \item applies appropriate hypothesis tests to judge qcCFs as representing putative hybridization,
#' resolved trees, or unresolved (star) trees using \code{alpha} and \code{beta} as significance levels,
#' \item produces a simplex plot showing results of the hypothesis tests for all qcCFs
#' \item computes the appropriate NANUQ distance table, writing it to a file.
#' }
#' The distance table file
#' can then be opened in the external software SplitsTree \insertCite{SplitsTree}{MSCquartets} (recommended) or within R using the package \code{phangorn} to
#' obtain a circular split system under the Neighbor-Net algorithm, which is then depicted as a splits graph.
#' The splits graph should be interpreted via
#' the theory of \insertCite{ABR19;textual}{MSCquartets} to infer the level-1 species network, or to conclude the data does
#' not arise from the NMSC on such a network.
#'
#' If \code{alpha} and \code{beta} are vectors, they must have the same length k. Then the i-th entries are paired to
#' produce k plots and k output files. This is equivalent to k calls to \code{NANUQ} with scalar values of \code{alpha} and \code{beta}.
#'
#' A call of \code{NANUQ} with \code{genedata} given as a table previously output from \code{NANUQ} is
#' equivalent to a call of \code{NANUQdist}. If \code{genedata} is a table previously output from \code{quartetTableResolved}
#' which lacks columns of p-values for hypothesis tests, these will be appended to the table output by \code{NANUQ}.
#'
#' If plots are produced, each point represents an empirical quartet concordance factor,
#' color-coded to represent test results.
#'
#' In general, \code{alpha} should be chosen to be small and \code{beta}
#' to be large so that most quartets are interpreted as resolved trees.
#'
#' Usually, an initial call to \code{NANUQ} will not give a good analysis, as values
#' of \code{alpha} and \code{beta} are likely to need some adjustment based on inspecting the data. Saving the returned
#' table from \code{NANUQ} will allow for the results of the time-consuming computation of qcCFs to be
#' saved, along with p-values,
#' for input to further calls of \code{NANUQ} with new choices of \code{alpha} and \code{beta}.
#'
#' See the documentation for \code{\link{quartetNetworkDist}} for an explanation of a small, rarely noticeable,
#' stochastic element of the algorithm.
#'
#' For data sets of many gene trees, user time may be reduced by using parallel code for
#' counting displayed quartets. See \code{\link{quartetTableParallel}}, where example commands are given.
#'
#'
#' @references
#' \insertRef{ABR19}{MSCquartets}
#'
#' \insertRef{SplitsTree}{MSCquartets}
#'
#' @param genedata gene tree data that may be supplied in any of 3 forms:
#' \enumerate{
#' \item as a character string giving the name of a file containing Newick gene trees,
#' \item as a multiPhylo object containing the gene trees, or
#' \item as a table of quartets on the gene trees, as produced by a previous call to
#' \code{NANUQ} or \code{quartetTableResolved}, which has columns only for taxa, resolved quartet counts,
#' and possibly p_T3 and p_star
#' }
#' @param outfile  a character string giving an output file name stub for
#' saving a \code{NANUQ} distance matrix in nexus format; to the stub \code{outfile}
#' will be appended an \code{alpha} and \code{beta} value and ".nex";
#' if \code{NULL} then then no file is written
#' @param omit \code{FALSE} to treat unresolved quartets as 1/3 of each resolution;
#' \code{TRUE} to discard unresolved quartet data; ignored if gene tree data given as quartet table
#' @param epsilon minimum for branch lengths to be treated as non-zero; ignored if gene tree data given as quartet table
#' @param alpha a value or vector of significance levels for judging p-values
#' testing a null hypothesis of no hybridization vs. an alternative of hybridization, for each quartet;  a smaller value applies
#' a less conservative test for a tree (more trees), hence a stricter requirement for desciding in favor of hybridization (fewer reticulations)
#' @param beta a value or vector of significance levels for judging p-values testing
#' a null hypothesis of a star tree (polytomy) for each quartet vs. an alternative of anything else; a smaller value applies a less conservative
#' test for a star tree (more polytomies), hence a stricter requirement for deciding in favor of a resolved tree or network;
#' if vectors, \code{alpha} and \code{beta} must have the same length
#' @param taxanames if \code{genedata} is a file or a multiPhylo object, a vector of a subset
#' of the taxa names on the gene trees
#' to be analyzed, if \code{NULL} all taxa on the first gene tree are used; if \code{genedata}
#' is a quartet table, this argument is ignored and all taxa in the table are used
#' @param plot \code{TRUE} produces simplex plots of hypothesis test results, \code{FALSE} omits plots
#'
#' @return a table \code{$pTable} of quartets and p-values for judging fit to the MSC on quartet
#' trees, and a distance table \code{$dist}, or list of distance tables, giving NANUQ distance (returned invisibly);
#' the table can be used as input to \code{NANUQ} or \code{NANUQdist} with new choices of alpha and beta, without re-tallying quartets on
#' gene trees; the distance table is to be used as input to NeighborNet.
#'
#'
#' @seealso \code{\link{quartetTable}}, \code{\link{quartetTableParallel}}, \code{\link{quartetTableDominant}}, \code{\link{quartetTreeTestInd}},
#' \code{\link{quartetStarTestInd}}, \code{\link{NANUQdist}}, \code{\link{quartetTestPlot}}, \code{\link{pvalHist}},
#' \code{\link{quartetNetworkDist}}
#'
#' @examples
#' data(pTableYeastRokas)
#' out=NANUQ(pTableYeastRokas, alpha=.05, beta=.95, outfile = NULL)
#' # Specifying an outfile would write the distance table to it for opening in SplitsTree.
#' # Alternately, to use the phangorn implementation of NeighborNet
#' # within R, enter the following additional lines:
#' nn=neighborNet(out$dist)
#' plot(nn,"2D")
#'
#' @export
NANUQ = function( genedata,
                  outfile = "NANUQdist",
                  omit = FALSE,
                  epsilon=0,
                  alpha = .05,
                  beta = .95,
                  taxanames = NULL,
                  plot = TRUE) {

  if (!(is.numeric(alpha) && is.numeric(beta))) {
    stop("Arguments alpha and beta must be numeric.")
  }


  if ("matrix" %in% class(genedata)) {
    pTable = genedata
    momit=missing(omit)
    mepsilon=missing(epsilon)
    mtaxanames=missing(taxanames)
    if (!momit | !mepsilon | !mtaxanames){
      warning(
        "Since genedata supplied as quartet table, ignoring arguments 'omit', 'epsilon', 'taxanames'."
      )}
      taxanames=NULL
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

  if (!("p_T3"%in% colnames(pTable))){
    pTable = quartetTreeTestInd(pTable, model = "T3", lambda = 0) # compute p-values for quartet CFs fit to MSC on a tree
  }

  if (!("p_star"%in% colnames(pTable))) {
  pTable = quartetStarTestInd(pTable)# compute p-values for each quartet CFs fit to star tree
  }

  dist=NANUQdist(pTable, outfile, alpha, beta, plot) # construct distance tables and put in files

  invisible(list(pTable=pTable, dist=dist))
}

####################################################################################################

#' Compute NANUQ distance and write to file
#'
#' Computes the quartet distance tables for the NANUQ algorithm of \insertCite{ABR19;textual}{MSCquartets}, using precomputed p-values for quartets,
#' for each of several levels specified. Distance tables are written to files, in nexus format.
#'
#' @details
#' If plots are produced, each point represents an empirical quartet concordance factor,
#' color-coded to represent test results giving interpretation as network, resolved tree, or star tree.
#'
#' If \code{alpha} and \code{beta} are vectors, they must be of the same length k. Then the i-th entries are
#' paired to produce k plots and k distance tables/output files. This is equivalent to k
#' calls to \code{NANUQdist} with paired scalar values from the vectors of \code{alpha} and \code{beta}.
#'
#' See the documentation for \code{\link{quartetNetworkDist}} for an explanation of a small, rarely noticeable,
#' stochastic element of the algorithm.
#'
#' @references
#' \insertRef{ABR19}{MSCquartets}
#'
#' @param pTable a table of resolved quartets and p-values, as previously computed by \code{NANUQ}, or by both \code{quartetTreeTestInd} and
#' \code{quartetStarTestInd}, with columns \code{"p_T3"} and \code{"p_star"}
#' @param outfile  a character string giving an output file name stub for
#' saving a \code{NANUQ} distance matrix in nexus format; to the stub \code{outfile}
#' will be appended an \code{alpha} and \code{beta} value and ".nex";
#' if \code{NULL} then not written to file
#' @param alpha a value or vector of significance levels for judging p-values
#' testing a null hypothesis of no hybridization for each quartet;  a smaller value applies
#' a more liberal test for a tree (more trees), hence a stricter requirement for suspecting hybridization (fewer reticulations)
#' @param beta a value or vector of significance levels for judging p-values testing
#' a null hypothesis of a star tree for each quartet; a smaller value applies a more liberal
#' test for a star tree (more polytomies), hence a stricter requirment for suspecting a resolved tree;
#' if vectors, \code{alpha} and \code{beta} must have the same length
#'
#' @param plot \code{TRUE} produces simplex plots of hypothesis tests, \code{FALSE} omits plots
#'
#' @return a NANUQ distance table, or a list of such tables if \code{alpha} and \code{beta}
#' are vectors (returned invisibly)
#'
#' @seealso \code{\link{NANUQ}}, \code{\link{quartetTreeTestInd}}, \code{\link{quartetStarTestInd}}
#'
#' @examples
#' data(pTableYeastRokas)
#' dist=NANUQdist(pTableYeastRokas, alpha=.05, beta=.95, outfile = NULL)
#'
#' @export
NANUQdist = function (pTable,
                      outfile = "NANUQdist",
                      alpha=.05,
                      beta=.95,
                      plot = TRUE) {
  if (!(is.numeric(alpha) && is.numeric(beta))) {
    stop("Critical values alpha and beta must be numeric.")
  }
  n = length(alpha)
  if (length(beta) != n)
    stop("Arguments alpha and beta must have same length")

  allQD=vector("list",n)

  for (i in 1:n) {
    if (plot == TRUE) {
      quartetTestPlot(pTable, "NANUQ", alpha[i], beta[i]) # plot quartets in simplex, indicating which are interpretted as hybridization
    }

    QD = quartetNetworkDist(pTable, alpha[i], beta[i])           # compute the quartet distance table

    if (!is.null(outfile)) {
    outfile1 = paste0(outfile, "_alpha", alpha[i], "_beta", beta[i], ".nex")  # name of the nexus output file
    nexusDist(QD, outfile1)           # Produce the file to be read by SplitsTree
    }
    allQD[[i]]=QD
  }
  if (length(allQD)==1) { #if one table, make it not a list
    allQD=allQD[[1]]
  }

  if (is.null(outfile)) {
    message("Distance table(s) not written to file(s).")
  }
  invisible(allQD)
}

##################################################################################

#' Write a distance table to a file in nexus format
#'
#' Write a distance table to a file in nexus format.
#'
#' @param distMatrix a square matrix giving a distance table, with rows and columns labeled by taxon names
#' @param outfilename the name of an output file
#'
#' @return NULL
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' DQT=quartetTableDominant(RQT)
#' Dist=quartetDist(DQT)
#' nexusDist(Dist,outfile = file.path(tempdir(), "NANUQdist"))
#'
#' @export
nexusDist = function(distMatrix,
                     outfilename) {
  taxanames = colnames(distMatrix)
  ntaxa = length(taxanames)
  taxanames = matrix(taxanames, ntaxa, 1)
  sink(outfilename)
  cat("#nexus\n BEGIN Taxa;\n DIMENSIONS ntax=",
      ntaxa,
      ";\n TAXLABELS\n")
  cat(paste(taxanames, "\n"))
  cat(
    ";\n END;\n \n BEGIN Distances; \n DIMENSIONS ntax=",
    ntaxa,
    "; \n FORMAT labels=left diagonal triangle=both;\n MATRIX\n"
  )
  cat(t(cbind(taxanames, distMatrix)))
  cat("; \n END;")
  sink()
 message("Distance table written to file: ", outfilename)
}

############################################################################

#' Compute network quartet distance between taxa
#'
#' Produce network quartet distance table for the NANUQ algorithm, from a table of quartets and p-values,
#' and specified levels of quartet hypothesis tests. The network quartet distance, which
#' is described more fully by \insertCite{ABR19;textual}{MSCquartets}, generalizes
#' the quartet distance of \insertCite{Rho19;textual}{MSCquartets}.
#'
#' @details In case of a triple of quartet counts with the two largest equal and the third slighltly smaller,
#' along with \code{alpha} and \code{beta} leading to a star quartet being rejected and a tree not being rejected,
#' this function chooses a resolved quartet topology uniformly at random from the two largest counts. This is the only
#' stochastic element of the code, and its impact is usually negligible.
#'
#' @references
#' \insertRef{ABR19}{MSCquartets}
#'
#' \insertRef{Rho19}{MSCquartets}
#'
#' @param pTable a table of quartets and p-values, as computed by NANUQ, or
#' \code{quartetTreeTestInd} and \code{quartetStarTestInd}
#' @param alpha a scalar significance level for judging p-values \code{p_T3} indicating hybridization on quartet;
#'      smaller value gives fewer hybridization decisions
#' @param beta a scalar significance level for judging p-values \code{p_star} indicating quartet star tree;
#'      smaller value gives fewer resolved tree decisions
#'
#' @return a distance table
#'
#' @examples
#' data(pTableYeastRokas)
#' dist=quartetNetworkDist(pTableYeastRokas, alpha=.05, beta=.95)
#'
#' @seealso \code{\link{NANUQ}}, \code{\link{NANUQdist}}
#'
#' @export
quartetNetworkDist = function(pTable,
                              alpha,
                              beta) {
  if (!(is.numeric(alpha) && is.numeric(beta))) {
    stop("Critical values alpha and beta must be numeric.")
  }
  M = dim(pTable)[1] # number of quartet concordance factors
  ntaxa = which(colnames(pTable)=="12|34")-1 #number of taxa

  taxanames = colnames(pTable)[1:ntaxa]  #names of taxa
  qnames = c("12|34", "13|24", "14|23")

  quartetDist = matrix(0, ntaxa, ntaxa)
  #alocates space for dissimilarity matrix
  rownames(quartetDist) = taxanames
  colnames(quartetDist) = taxanames

  for (m in 1:M) {
    # consider each set of 4 taxa

    taxanums = which(pTable[m, 1:ntaxa] == 1) #determine taxa
    qcounts = pTable[m, qnames] # get counts

    if (pTable[m, "p_star"] > beta) {
      # if quartet judged as star tree
      quartetDist[taxanums, taxanums] = quartetDist[taxanums, taxanums] + 1 #  it separates all taxa
    } else {
      qcountspr = qcounts + runif(3) * .0001 # introduce random tie breaking
      majorq = which.max(qcountspr) #determine quartet with largest count
      if (pTable[m, "p_T3"] > alpha) {
        # if high p-value indicates tree
        if (majorq == 1) {
          quartetDist[taxanums[1:2], taxanums[3:4]] = quartetDist[taxanums[1:2], taxanums[3:4]] +
            1
        } else {
          if (majorq == 2) {
            quartetDist[taxanums[c(1, 3)], taxanums[c(2, 4)]] = quartetDist[taxanums[c(1, 3)], taxanums[c(2, 4)]] +
              1
          } else {
            quartetDist[taxanums[c(1, 4)], taxanums[c(2, 3)]] = quartetDist[taxanums[c(1, 4)], taxanums[c(2, 3)]] +
              1
          }
        }
      } else {
        # if low p-value indicates cycle
        minorq = which.min(qcountspr)
        middleq = setdiff(1:3, c(majorq, minorq))
        if (majorq == 1 || middleq == 1) {
          quartetDist[taxanums[1:2], taxanums[3:4]] = quartetDist[taxanums[1:2], taxanums[3:4]] +
            .5
        }
        if (majorq == 2 || middleq == 2) {
          quartetDist[taxanums[c(1, 3)], taxanums[c(2, 4)]] = quartetDist[taxanums[c(1, 3)], taxanums[c(2, 4)]] + .5
        }
        if (majorq == 3 || middleq == 3) {
          quartetDist[taxanums[c(1, 4)], taxanums[c(2, 3)]] = quartetDist[taxanums[c(1, 4)], taxanums[c(2, 3)]] + .5
        }
      }
    }
  }
  quartetDist = quartetDist + t(quartetDist)
  quartetDist = 2 * quartetDist + 2 * ntaxa - 4 # adjust to get appropriate distances off-diagonal
  quartetDist = quartetDist - diag(diag(quartetDist)) # but reset diagonal to 0

  return(quartetDist)
}
