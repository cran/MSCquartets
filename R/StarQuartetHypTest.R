#' Multiple independent hypothesis tests for gene quartet counts fitting 
#' a species quartet star tree under the MSC
#'
#' Perform hypothesis tests for a species quartet star tree vs. any alternative for all quartet counts in an input table, 
#' as if the quartets are independent. 
#' 
#' @details This function assumes all quartets are resolved. 
#' The test performed is described in \code{quartetStarTest}.
#'
#' @param rqt  Table of resolved quartet counts, as produced by \code{quartetTableResolved}, or \code{quartetTreeTestInd}
#'
#' @return the same table as the input \code{rqt} with column \code{"p_star"} appended, containing p-values for 
#' judging fit to MSC on a star tree
#'
#' @seealso \code{\link{quartetStarTest}}, \code{\link{quartetTreeTest}}, \code{\link{quartetTreeTestInd}},
#' \code{\link{quartetTableResolved}}, \code{\link{quartetTestPlot}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' pTable=quartetStarTestInd(RQT)
#' quartetTablePrint(pTable[1:6,])
#'
#' @export
quartetStarTestInd = function(rqt) {
  M = dim(rqt)[1]
  qcols = c("12|34", "13|24", "14|23")
  rqt = cbind(rqt, p_star = 0)
  
  message("Applying hypothesis test for star tree model to ",M," quartets.")
  for (i in 1:M) {
    obs = rqt[i, qcols] # get quartet counts
    rqt[i, "p_star"] = quartetStarTest(obs) # store p-value for H0=star tree
  }
  return(rqt)
}

#' Hypothesis test for quartet counts fitting a star tree under the MSC
#'
#' Perform hypothesis test for star tree for a vector of  quartet counts to fit expected frequencies of (1/3,1/3,1/3).
#' The test performed is a standard Chi-square with 2 degrees of freedom.
#'
#' @param obs  vector of 3 counts of resolved quartet frequencies
#' @return p-value
#'
#' @examples
#' obs=c(16,72,12)
#' quartetStarTest(obs)
#'
#' @export
quartetStarTest = function(obs) {
  z = chisq.test(obs) #test agaist uniform expectated values
  p = z$p.value
  return(p)
}
