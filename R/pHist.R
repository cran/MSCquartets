#' Plot histogram of log p-values in table
#'
#' Graphical exploration of extreme p-values from quartet hypothesis tests, to aid in choosing critical values
#' for hypothesis tests. Log base 10 of p-values exceeding some minimum are plotted, to explore gaps in 
#' the tail of the distribution. 
#' 
#' @details 
#' Since logarithms are plotted, p-values close to 0 will appear as negative numbers of large magnitude, putting the tail of the distribution
#' to the left in the histogram.
#' 
#' When exploring possible critical values for the hypothesis tests in the NANUQ algorithm, use \code{model= "T3"} to 
#' choose values for \code{alpha} which distinguish treelikeness from hybridization, and \code{model= "star"} to 
#' choose values for \code{beta} which distinguish polytomies from resolved trees. 
#' In general, \code{alpha} should be chosen to be small and \code{beta}
#' to be large so that most quartets are interpreted as resolved trees.
#'
#'@param pTable a quartet table with p-values, such as from \code{NANUQ}, \code{quartetTreeTestInd}, 
#'or \code{quartetStarTestInd}
#'@param model  one of \code{"T1"}, \code{"T3"}, or \code{"star"}, where \code{pTable} 
#'contains a column \code{p_model} of p-values
#'@param pmin include only p-values above \code{pmin} in the histogram
#'
#'@return No return value, called for side effects
#'
#'@examples
#' pTable=NANUQ(system.file("extdata","dataYeastRokas",package="MSCquartets"), 
#'        alpha=0.05, beta=.95, outfile = file.path(tempdir(), "NANUQdist"))
#' pvalHist(pTable,"T3")
#'
#'@seealso \code{\link{NANUQ}}, \code{\link{NANUQdist}}
#'
#'@export
pvalHist = function(pTable, 
                    model, 
                    pmin = 0) {
  if (!(model %in% c("T1","T3","star"))) {
    stop('Argument model must be one of "T1", "T3", or "star".')
  }
  pcol=paste0("p_",model) 
  if (!(pcol %in% colnames(pTable))) {
    stop(c('Argument pTable has no column ',pcol,'.'))
  }
  pvals = pTable[, pcol]
  vals = pvals[which(pvals >= pmin)]
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mar = c(5.1, 4.1, 4.1, 2.1))# set margin
  hist(log10(vals),
       main = bquote("Histogram of log"[10] * "(" * .(pcol) * "), " * .(pcol) >
                       .(pmin)))
}
