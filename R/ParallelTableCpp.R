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
#' gtrees=read.tree(file=system.file("extdata","dataHeliconiusMartin",package="MSCquartets"))
#' QT=quartetTableParallel(gtrees,numCores=2)
#' RQT=quartetTableResolved(QT)
#' pTable=NANUQ(RQT,alpha=1e-40, beta=1e-30, outfile = file.path(tempdir(), "NANUQdist"))}
#'
#' @importFrom ape cophenetic.phylo compute.brlen
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#'
#' @export
quartetTableParallel <- function(trees,
                                 taxonnames = NULL,
                                 epsilon = 0,
                                 numCores) {

  if ( !(numCores>=2) ) stop("Argument 'numCores' must be at least 2.")
  numcores=floor(numCores)

  # function for combining parallel call results
  combQuartets <-function(quartets1,quartets2){
    quartets1[,c("12|34","13|24","14|23","1234")]=
      quartets1[,c("12|34","13|24","14|23","1234")] +
      quartets2[,c("12|34","13|24","14|23","1234")]
    return(quartets1)
  }

  # main code for parallel calls
  start_time = Sys.time()

  if (is.null(taxonnames)) {
    taxonnames = trees[[1]]$tip.label
    message("Using taxa that appear on first gene tree.")
  }
  taxonnames = sort(taxonnames)
  nt = length(trees)
  N = length(taxonnames)

  M = choose(N, 4)
  Q = matrix(0, M, N)
  colnames(Q) = c(taxonnames)

  registerDoParallel(numCores)  # set number of cores to be used

  numTrees=length(trees)
  numTs=ceiling(numTrees/numCores) #number of trees for each core (except last)

  i=0 #to avoid a CRAN check Note

  qTable <- foreach::foreach (i = 1:numCores, .combine = combQuartets ) %dopar% {
   Ts=trees[(1+(i-1)*numTs):min((i*numTs),numTrees)]
   suppressWarnings(quartetTable(Ts, taxonnames=taxonnames, epsilon = epsilon))
  }

  rSums=rowSums(qTable)
  if (sum(!(rSums=rSums[1])) >0)
    warning("Some taxa missing from some trees.")
  if ((N > 4) && (sum(rSums == 0) > 0))
    warning("Some 4-taxon set(s) not present on any tree.")
  if (sum(qTable[, "1234"]) > 0)
    warning("Some quartets unresolved.")

  current_time = Sys.time()
  elapsedTime = as.numeric(difftime(current_time, start_time,
                                    units = "mins"))
  if (elapsedTime > 15) {
    message("Time to process quartets on gene trees was ",
            elapsedTime," minutes.")
  }

  return(qTable)
}
