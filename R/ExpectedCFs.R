#' Produce table of expected quartet concordance factors for a species tree
#'
#' Compiles table of expected quartet concordance factors (CFs) for gene trees under the MSC model 
#' on a metric species tree.
#'
#' @details
#' The species tree may be rooted or unrooted, but must have edge lengths given in coalescent units. Note that while the MSC
#' requires a rooted metric species tree parameter, as shown in \insertCite{Allman2011}{MSCquartets} the quartet CFs are independent 
#' of the root.
#'
#' In the returned table, columns are labeled by taxon names and quartet names ("12|34", etc.).
#' 1s and 0s in taxon columns indicate the taxa in a quartet. Quartet 12|34 means the first and second 
#' indicated taxa form a cherry, 13|24 means the first and third form a cherry, and 14|23 means
#' the first and fourth form a cherry. Unresolved quartets always have expectation 0 under the MSC.
#'
#' If a simplex plot is produced, for the \code{T1} model all CFs will lie on the vertical model line, 
#' and many choices of 4 taxa will give the same CFs. For model \code{T3} the model lines the CFs are 
#' plotted on depend on taxon names and are essentially arbitrary.
#'
#' @param speciestree phylo or character object specifying un/rooted metric species tree
#' @param plot \code{TRUE} (default) to produce simplex plot of CFs, or \code{FALSE} to omit plot
#' @param model "T1" or "T3" specifying model for plot
#' @param cex scaling factor for size of plotted symbols
#'
#' @return
#'     an (\code{n} choose 4)x(\code{n}+3) matrix encoding
#'     4 taxon subsets of taxa and expectation of each of the
#'     quartets 12|34, 13|24, 14|23 across gene trees
#'
#' @seealso \code{\link{quartetTable}}, \code{\link{quartetTableResolved}}
#'
#' @examples
#' stree="((((t5:5000,t6:5000):5000,t4:10000):2500,t7:12500):7500,((t8:3000,t9:3000):5000,
#' ((t1:4000,t2:4000):2500,t3:6500):1500):12000);"
#' st=read.tree(text=stree)
#' st$edge.length=st$edge.length/10000
#' expectedCFs(st)
#'
#' @references
#' 
#' \insertRef{Allman2011}{MSCquartets}
#' 
#' @importFrom methods is
#' @importFrom ape unroot read.tree cophenetic.phylo
#'
#' 
#' @export
expectedCFs <- function(speciestree,
                        plot = TRUE,
                        model = "T1",
                        cex = 1) {
  if (methods::is(speciestree,"phylo")) {
    stree = ape::unroot(speciestree)
  } else {
    if (methods::is(speciestree,"character")) {
      stree = ape::unroot(read.tree(text = speciestree)) #unrooted species tree
    } else {
      stop('Species tree must be of class "phylo"" or class "character".')
    }
  }
  if (is.null(stree$edge.length)) {
    stop('Species tree must have branch lengths.')
  }
  if (plot == TRUE & isFALSE(model %in% c("T1", "T3"))) {
    stop('Argument "model" must be "T1" or "T3" for plotting.')
  }
  
  D = ape::cophenetic.phylo(stree) # create distance table, for determining displayed quartets
  D = D[order(rownames(D)), order(colnames(D))]# put taxa in sorted order
  
  
  taxonnames = row.names(D)
  N = length(taxonnames)
  
  # set up table
  M = choose(N, 4)
  Q = matrix(0, M, N + 3) #allocate space for table
  qnames = c("12|34", "13|24", "14|23")
  colnames(Q) = c(taxonnames, qnames) #create column names
  
  # fill out table
  m = 0
  for (i in 1:(N - 3)) {
    # for each 4-taxon set
    for (j in (i + 1):(N - 2)) {
      for (k in (j + 1):(N - 1)) {
        for (l in (k + 1):N) {
          m = m + 1
          Q[m, c(i, j, k, l)] = 1 #encode set
          a = D[i, j] + D[k, l]
          b = D[i, k] + D[j, l]
          c = D[i, l] + D[j, k]
          
          if ((a <= b) & (a <= c)) {
            x = ((b + c) - 2 * a) / 2
            y = (1 / 3) * exp(-x)
            Q[m, (N + 1):(N + 3)] = c(1 - 2 * y, y, y)
          } else {
            if ((b <= a) & (b <= c)) {
              x = ((a + c) - 2 * b) / 2
              y = (1 / 3) * exp(-x)
              Q[m, (N + 1):(N + 3)] = c(y, 1 - 2 * y, y)
            } else {
              x = ((a + b) - 2 * c) / 2
              y = (1 / 3) * exp(-x)
              Q[m, (N + 1):(N + 3)] = c(y, y, 1 - 2 * y)
            }
          }
        }
      }
    }
  }
  
  
  
  if (plot == TRUE) {
    if (model == "T3") {
      simplexPrepare(model = model,
                     maintitle = "Expected Quartet Concordance Factors",
                     titletext = "Model T3")
      for (m in 1:M) {
        simplexPoint(Q[m, (N + 1):(N + 3)],
                     cex = cex,
                     pch = 8,
                     col = "blue")
      }
    } else {
      simplexPrepare(model = model,
                     maintitle = "Expected Quartet Concordance Factors",
                     titletext = "Model T1")
      for (m in 1:M) {
        simplexPoint(
          sort(Q[m, (N + 1):(N + 3)], decreasing = TRUE),
          cex = cex,
          pch = 8,
          col = "blue"
        )
      }
    }
  }
  return(Q)
}
