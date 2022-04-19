#' Produce table of counts of quartets displayed on trees, in parallel for large data sets
#'
#' Compiles table of quartet count concordance factors (qcCFs) for topological quartets displayed on a
#' collection of trees. Gives the same output as \code{quartetTable}, but operates in parallel.
#'
#' @details
#'
#' The number of available cores can be determined by \code{parallel::detectCores()}.  
#' With overhead, tabulating quartets for a large data set (many taxa and/or many gene trees) on a 4-core 
#' computer using \code{numCores=4} may require less than half the elapsed time of the sequential \code{quartetTable}.
#' 
#' The names in \code{taxonnames} may be any subset of those on the trees.
#' Branch lengths of non-negative size less than or equal to \code{epsilon}
#' are treated as zero, giving polytomies.
#'
#' In the returned table, columns are labeled by taxon names and quartet names ("12|34", etc.).
#' 1s and 0s in taxon columns indicate the taxa in a quartet. Quartet 12|34
#' means the first and second indicated taxa form a cherry, 13|24 means the first and third form a cherry, 14|23 means
#' the first and fourth form a cherry, and 1234 means the quartet is unresolved.
#'
#' An error occurs if any branch length is negative.
#' Warnings are given if some of \code{taxonnames} are missing on some trees, or
#' if some 4-taxon set is not on any tree.
#'
#' If \code{random}>0, then for efficiency \code{random} should be much smaller then
#' the number of possible 4 taxon subsets.
#' 
#' If the quartet counts are to be used for NANUQ, or any other routines requiring resolved quartet counts,
#' \code{\link{quartetTableResolved}} must be run following \code{quartetTableParallel}. See example below.
#'
#' @param trees multiphylo object containing un/rooted metric/topological trees
#' @param taxonnames vector of \code{n} names of taxa of interest; if \code{NULL} then taken from taxa on \code{trees[[1]]}
#' @param epsilon minimum for branch lengths to be treated as non-zero
#' @param random number of random subsets of 4 taxa to consider; if 0, use all \code{n} choose 4 subsets
#' @param numCores number of cores to use for parallel calls
#' @return
#'     an (\code{n} choose 4)x(\code{n}+4) matrix (or (\code{random})x(\code{n}+4) matrix) encoding
#'     4 taxon subsets of \code{taxonnames} and counts of each of the
#'     quartets 12|34, 13|24, 14|23, 1234 across the trees
#'
#'
#' @seealso \code{\link{quartetTable}}, \code{\link{quartetTableResolved}}, \code{\link{quartetTableDominant}}, \code{\link{taxonNames}}
#'
#' @examples \donttest{
#
#' gtrees=read.tree(file=system.file("extdata","dataHeliconiusMartin",package="MSCquartets"))
#' QT=quartetTableParallel(gtrees,numCores=2)
#' RQT=quartetTableResolved(QT)
#' pTable=NANUQ(RQT,alpha=1e-40, beta=1e-30, outfile = file.path(tempdir(), "NANUQdist"))}
#'
#' @importFrom ape cophenetic.phylo compute.brlen
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
#' @export
quartetTableParallel <- function(trees,
                                 taxonnames = NULL,
                                 epsilon = 0,
                                 random = 0,
                                 numCores) {
  if (random < 0)
    stop("Parameter 'random' must be non-negative.")
  random = round(random)
  start_time <- Sys.time()
  
  #function to find all quartets displayed on 1 tree, to be called in parallel
  quartetDisplayed <- function (quartetIndices,
                                tree,
                                epsilon = 0) {
    N = ncol(quartetIndices)
    M = nrow(quartetIndices)
    Q = matrix(0, M, 4)
    colnames(Q) = c("12|34", "13|24", "14|23", "1234")
    taxonnames = colnames(quartetIndices)
    warnMissing = 0
    
    if (is.null(tree$edge.length)) {
      tree = compute.brlen(tree, 1)
    }
    else {
      if (sum(tree$edge.length < 0) > 0)
        stop("Error: Negative branch length in tree")
    }
    zeros = which(tree$edge.length <= epsilon)
    tree$edge.length[] = 1
    tree$edge.length[zeros] = 0
    d = cophenetic.phylo(tree)
    for (m in 1:M) {
      tax = as.character(taxonnames[which(quartetIndices[m, 1:N] ==
                                            1)])
      if (all(tax %in% colnames(d))) {
        a = d[tax[1], tax[2]] + d[tax[3], tax[4]]
        b = d[tax[1], tax[3]] + d[tax[2], tax[4]]
        c = d[tax[1], tax[4]] + d[tax[2], tax[3]]
        z = sort(c(a, b, c))
        if (z[1] == z[2]) {
          Q[m, "1234"] = Q[m, "1234"] + 1
        }
        else {
          if (z[1] == a) {
            Q[m, "12|34"] = Q[m, "12|34"] + 1
          }
          else {
            if (z[1] == b) {
              Q[m, "13|24"] = Q[m, "13|24"] + 1
            }
            else {
              Q[m, "14|23"] = Q[m, "14|23"] + 1
            }
          }
        }
      }
      else
        warnMissing = 1
    }
    return(list(Q, warnMissing))
  }
  
  # function for combining parallel call results
  combQuartets <-function(quartetsAndFlag1,quartetsAndFlag2){
    return(mapply('+',quartetsAndFlag1,quartetsAndFlag2))
  }
  
  # main code for parallel calls
  
  if (is.null(taxonnames)) {
    taxonnames = trees[[1]]$tip.label
    message("Using taxa that appear on first gene tree.")
  }
  taxonnames = sort(taxonnames)
  nt = length(trees)
  N = length(taxonnames)
  
  if (random > 0)
    M = random
  else
    M = choose(N, 4)
  Q = matrix(0, M, N)
  colnames(Q) = c(taxonnames)
  warnMissing = 0
  
  message(
    "Counting occurrences of displayed quartets for ",M,
    " four-taxon subsets of ",N,
    " taxa across ",nt," gene trees."
  )
  
  if (random == 0) {
    #encode all subsets of 4 taxa
    m = 0
    for (i in 1:(N - 3)) {
      for (j in (i + 1):(N - 2)) {
        for (k in (j + 1):(N - 1)) {
          for (l in (k + 1):N) {
            m = m + 1
            Q[m, c(i, j, k, l)] = 1
          }
        }
      }
    }
  }
  else { # encode random subsets of 4 taxa
    i = 1
    while (i <= random) {
      q = sample(N, size = 4, replace = FALSE)
      row = integer(N)
      row[q] = 1
      j = 1
      while (j < i) {
        if (identical(Q[j, 1:N], row)) # check if this subset already chosen
          j = i + 1
        else
          j = j + 1
      }
      if (j == i) { # new subset, so record it
        Q[i, 1:N] = row 
        i = i + 1
      }
    }
  }
  
  registerDoParallel(numCores)  # set number of cores to be used
  
 cfM <- foreach (i = 1:length(trees), .combine = 'combQuartets') %dopar% {
    quartetDisplayed(Q, trees[[i]], epsilon = epsilon)
  }
  
 stopImplicitCluster()
 
 
  cf = cfM[[1]]
  warnMissing = cfM[[2]]
  quartetTable <- cbind(Q,cf)
  
  
  if (warnMissing > 0)
    warning("Some taxa missing from some trees.")
  if ((N > 4) && (sum(rowSums(cf) == 0) > 0))
    warning("Some 4-taxon set not present on any tree.")
  if (sum(cf[, "1234"]) > 0)
    warning("Some quartets unresolved.")
  
  current_time = Sys.time()
  elapsedTime = as.numeric(difftime(current_time, start_time,
                                    units = "mins"))
  if (elapsedTime > 15) {
    message("Time to process quartets on gene trees was ",
            elapsedTime," minutes.")
  }
  
  return(quartetTable)
}
