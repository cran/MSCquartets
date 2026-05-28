#' Build a Species Tree from Gene Trees Using the Quartet Distance
#'
#' This function takes a collection of gene trees and uses a
#' species tree by first computing the quartet distance of 
#' \insertCite{Rho19;textual}{MSCquartets}
#' from the gene trees and then applying a distance method.
#'
#' @references 
#' \insertRef{Rho19}{MSCquartets}
#' 
#' @param genedata gene tree data, supplied in any of 3 forms:
#' \enumerate{
#' \item as a character string giving the name of a file containing Newick gene trees,
#' \item as a \code{multiPhylo} object containing the gene trees, or
#' \item as a table of quartet counts on the gene trees, as produced by
#' \code{quartetTable} or \code{quartetTableResolved}
#' }
#' @param method a function for constructing a tree from a distance matrix, such as \code{NJ} or
#' \code{fastme.bal}
#' @param taxanames used only if \code{genedata} is a filneame or multiphylo object;
#' restrict gene trees to these taxa; if \code{NULL}, use taxa on first tree
#' @param omit used only if \code{genedata} is a filneame,\code{multiPhylo} object, 
#' or a quartet table with column \code{'1234'};
#' if TRUE omit unresolved quartets, if \code{FALSE} treat them as 1/3 of each resolution
#'
#' @return A \code{phylo} object representing the inferred tree.
#'
#' @seealso
#' \code{quartetDist}, \code{QDS}, \code{quartetTableDominant}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tree<-quartetDistTree(gtrees)
#' write.tree(tree)
#' plot(tree)
#'
#' @export
quartetDistTree = function(genedata,
                           method = ape::fastme.bal,
                           taxanames = NULL,
                           omit = FALSE) 
  {
  if ("matrix" %in% class(genedata)) {
    QT = genedata
    if (!(is.null(taxanames))) {
      message(
        "Since genedata supplied as quartet table, ignoring argument 'taxanames'."
      )
    }
  } else  {
    if ("multiPhylo" %in% class(genedata))  {
      gtrees = genedata
    } else {
      if ("character" %in% class(genedata)) {
        gtrees <- read.tree(genedata) #read gene trees
        message(paste("Read", length(gtrees), "gene trees from file."))
      } else {
        stop("Data must be supplied as an object of type multiPhylo, character, or matrix.")
      }
    }
    if (is.null(taxanames)) {
      # if no taxa names specified,
      taxanames = gtrees[[1]]$tip.label   # ... get them from first tree
    }
    taxanames = sort(taxanames)
    if (length(taxanames) <= 25) {
      namelist = paste0(taxanames, collapse = ", ")
    } else {
      namelist = paste0(paste0(taxanames[1:25], collapse = ", "),
                        ",...(see output table for full list)")
    }
    message("Analyzing ", length(taxanames), " taxa: ", namelist)
    QT <- quartetTable(gtrees, taxanames)
    
  }
  QT <- quartetTableResolved(QT, omit = omit)
  QT <- quartetTableDominant(QT)
  dist <- quartetDist(QT)
  tree <- method(dist)
  return(tree)
}


  
  
  

