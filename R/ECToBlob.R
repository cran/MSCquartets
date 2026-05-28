#' Combine p-values of tests
#'
#' Given a vector of p-values from not-necessarily-independent tests, compute a combined p-value
#' using the Cauchy combination method Liu and Xie (2020), the Bonferroni method (also called MinP, the minimum times the number of tests),
#' or the combined Cauchy/Bon variants BBC and CBC as proposed by Chen (2022).
#'
#' All \code{NA}a in vector are ignored.
#'
#' @references
#' \insertRef{LiuXie2020}{MSCquartets}
#'
#' \insertRef{Chen2022}{MSCquartets}
#' 
#' @param pvec a vector of p-values to be combined
#' @param method the combination method to be used, one of \code{"Cauchy"}, \code{"Bon"} (default), \code{"BBC"}, \code{"CBC"}, or \code{"all"}
#'
#' @returns If the argument \code{method} is \code{"Cauchy"}, \code{"Bon"}, \code{"BBC"}, or \code{"CBC"} a single p-value obtained from that method is returned.
#' If \code{method} is "All", a vector of all 4 are returned, with entries names "p_Cauchy", "p_Bon", "p_BBC", and "p_CBC"
#'
#' @examples
#' combineP(c(.01,.04,.11))
#' combineP(c(.01,.04,.11),"all")
#'
#'@export
combineP = function(pvec, method = "Bon")
{
  if (!(method %in% c("Cauchy", "Bon", "BBC", "CBC", "all")))  {
    stop("Invalid argument: 'method' ")
  }
  
  pvec=pvec[which(!is.na(pvec))] #remove NAs
  
  if (length(pvec) == 0)
    pvec = 1

  #Cauchy
  if (method != "Bon") {
    is.small = (pvec < 1e-15)
    tvec = rep(NULL, length(pvec))
    tvec[!is.small] = tan((0.5 - pvec[!is.small]) * pi)
    tvec[is.small] = 1 / (pvec[is.small] * pi)
    tstat = mean(tvec)
    pCauchy = pcauchy(tstat, lower.tail = F)
  }

  #Bon
  if (method != "Cauchy")
    pBon = min(1, length(pvec) * min(pvec))


  #BBC
  if (method %in% c("BBC", "all"))
    pBBC = min(1, 2 * min(pCauchy, pBon))

  #CBC
  if (method %in% c("CBC", "all")) {
    pvec2 = c(pCauchy, pBon)
    is.small = (pvec2 < 1e-15)
    tvec = rep(NULL, 2)
    tvec[!is.small] = tan((0.5 - pvec2[!is.small]) * pi)
    tvec[is.small] = 1 / (pvec2[is.small] * pi)
    tstat = mean(tvec)
    pCBC = pcauchy(tstat, lower.tail = F)
  }

  if (method == "all")
    pall = c(
      p_Cauchy = pCauchy,
      p_Bon = pBon,
      p_BBC = pBBC,
      p_CBC = pCBC
    )

  p = get(paste0("p", method))
  return(p)
}



##########################################

#' Hypothesis testing of quartet counts for ECToBlob
#'
#' Append p-values for the star and T1 tests of displayed quartets on a given
#' tree to resolved quartet count table. This is used in \code{ECToBlob}
#' to obtain an estimated Tree of Blobs for a data set from a supplied tree that
#' is assumed to be a resolution of the Tree of Blobs. Three columns are appended
#' to the quartet count table: "p_T1-D" holds the T1 test p-values, "qindexD" holds
#' the index (1,2,or 3) of the quartet topology displayed on the tree, and "p_star"
#' holds the star tree test p-value.
#'
#' If on the tree a quartet is unresolved, p_T1 and qindex column entry will be NA
#'
#' \code{tree} must be supplied as phylo object, but may have multifurcations.
#'
#' If either p_star or p_T1-D and qindex-D values are already in table, they are
#' not recomputed.
#'
#' @param rqt either a table of resolved quartet counts from a collection of gene trees,
#' or that table with additional columns of quartet p-values.
#' @param tree an unrooted tree determining the topologies of the displayed quartet trees to test
#' @param lambda parameter for power-divergence statistic (e.g., 0 for likelihood ratio statistic, 1 for Chi-squared statistic)
#' @param method "MLest", "conservative", or "bootstrap"
#' @param smallsample	"precomputed" or "bootstrap", method of obtaining p-value when sample is small (<30)
#' @param smallcounts	"precomputed" or "bootstrap", method of obtaining p-value when some (but not all) counts are small
#' @param bootstraps	number of samples for bootstrapping
#'
#' @returns the same table as \code{rqt} with three columns appended:
#' \code{qindex-D} containing the indicies (1,2,or 3)
#' of the quartet topology displayed on \code{tree}, \code{p_T1-D}
#' containing the p-value for the T1 test of the displayed quartet topologies.
#' \code{p_star} containing p-values for the star tree test.
#'
#' @examples
#' data(ASTRALtreeLeopardusLescroart)
#' Atree=read.tree(text=ASTRALtreeLeopardusLescroart)
#' data(tableLeopardusLescroart)
#' pTable=quartetECtestInd(tableLeopardusLescroart,Atree)
#'
#' @export
quartetECtestInd <- function (rqt,
                              tree,
                              lambda = 0,
                              method = "MLest",
                              smallsample = "precomputed",
                              smallcounts = "precomputed",
                              bootstraps = 10 ^ 4)
{
  if (is.null(tree) ||
      is.numeric(lambda) == FALSE || is.vector(lambda) == FALSE ||
      length(lambda) != 1 ||
      smallsample %in% c("bootstrap", "precomputed") == FALSE ||
      smallcounts %in% c("bootstrap", "precomputed") == FALSE ||
      is.numeric(bootstraps) ==  FALSE ||
      is.vector(bootstraps) == FALSE ||
      length(bootstraps) != 1 ||
      bootstraps <= 0 || bootstraps %% 1 != 0 ||
      method %in% c("MLest", "conservative", "bootstrap") == FALSE)
  {
    stop("Invalid argument(s)")
  }

  colnames = colnames(rqt)
  taxanames = tree$tip.label

  n = length(taxanames)
  if (!(all.equal(sort(colnames[1:n]), sort(taxanames))))
    stop("Taxa in table do not match those on tree.")

  M = dim(rqt)[1]
  tree = di2multi(unroot(tree)) #allow multifurcating tree
  nedges = dim(tree$edge)[1]
  tree$edge.length = rep(1, nedges) #to compute displayed quartets by 4pt, make all edges length 1
  D = cophenetic.phylo(tree)
  D = D[order(rownames(D)), order(colnames(D))] #order distance table alphabetically

  pTable = rqt
  if (all(c("p_T1-D", "qindex-D") %in%  colnames(pTable)))
    message("Not recomputing T1 p-values for displayed trees.")
  else
  {
    pTable = cbind(pTable, ct = 0, qi = 0)
    colnames(pTable)[which(colnames(pTable) == "ct")] = "p_T1-D"
    colnames(pTable)[which(colnames(pTable) == "qi")] = "qindex-D"

    message("Applying hypothesis tests for T1 model to ", M, " quartets.")
    for (m in 1:M) {
      #run through sets of 4 taxa; determine from 4pt condition what topology is displayed on tree
      qnames = which(rqt[m, 1:n] == 1)
      a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
      b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
      c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
      if (min(c(a, b, c)) == max(c(a, b, c))) {
        #unresolved quartet
        pTable[m, "qindex-D"] = NA
        pTable[m, "p_T1-D"] = NA
      }
      else {
        treeQuartet = which.min(c(a, b, c)) #for test, put count corresponding to tree topology first
        qcounts = rqt[m, c("12|34", "13|24", "14|23")]
        temp = qcounts[1]
        qcounts[1] = qcounts[treeQuartet]
        qcounts[treeQuartet] = temp
        pTable[m, "qindex-D"] = treeQuartet
      

      pvec = quartetTreeTest(
        unname(qcounts),
        model = "T1",
        lambda = lambda,
        method = method,
        smallsample = smallsample,
        smallcounts = smallcounts,
        bootstraps = bootstraps
      )
      pTable[m, "p_T1-D"] = pvec$p.value
      }
    }
  }
  if ("p_star" %in% colnames(pTable))
    # compute p_star if necessary
    message("Not recomputing p_star values.")
  else
    pTable = quartetStarTestInd(pTable)

  return(pTable)
}

####################################################################

#' Edge Contraction to infer a Tree of Blobs
#'
#' Infer the tree of blobs using edge contractions
#'
#' Given quartet gene tree count data from a network, and a starting tree that is 
#' a resolution of its tree of blobs, contract edges for which there is poor 
#' support due to, first, star-like quartet signals and then, successively the 
#' most extreme non-tree-like signals.
#'
#' p-values for star and T1 tests for individual quartets displayed on the tree are
#' combined by one of 4 test correction methods suitable for non-independent tests. Which quartet
#' p-values are combined depends on one of 3 modes indicated by argument \code{qType}.
#'
#' The input tree should be computed by ASTRAL, or another heuristic known to
#' produce a resolution of the tree of blobs in ideal circumstances.
#'
#' Plot levels allow for several graphical depictions of the edge contract process,
#' from the starting tree to a final star tree.
#'
#' For testing with simulated data, a reference tree (e.g., the true Tree of Blobs) 
#' may be provided and the RF distance between it and each tree will be shown.
#'
#' @references
#' \insertRef{ECTo2026}{MSCquartets}
#'
#' \insertRef{MAR19}{MSCquartets}
#'
#' \insertRef{LiuXie2020}{MSCquartets}
#'
#' \insertRef{Chen2022}{MSCquartets}
#' 
#' @param genedata one of:
#' \enumerate{
#' \item a character string giving the name of a file containing a list of newick 
#' gene trees (see argument \code{omit}),
#' \item a multiphylo object containing a collection of gene trees (see argument \code{omit}),
#' \item a table of resolved quartet counts  from a collection of
#' gene trees, or that table with additional columns of quartet p-values. If p_star
#' or both p_T1-D and qindex-D are present (such as from a previous call to this function) 
#' they are not recomputed.}
#' @param tree a tree as a phylo object, or a character string giving the name 
#' of a file with a newick tree, assumed to be a resolution of the true
#' tree of blobs
#' @param qType a choice of which quartet p-values are combined to determine one for a
#' tree edge, "bi" for all initial bipartition quartets, "quad" for all initial quadripartition
#' quartets, or "mul" for all multipartition quartets at current contracting step
#' @param testCorrection a method for combining possibly dependent p-values to 
#' get a single p-value, one of "Bon", "Cauchy", "BBC", "CBC"
#' @param alpha test level for choosing inferred tree of blobs
#' @param beta test level for collapsing edges showing no resolution
#' @param colorCutoff level for changing color of edge p-values for T1 test
#' @param plot 0 for no plot;
#'             1 to plot all trees as edges are contracted, along with plot of tree p-values;
#'             2 to add -log_10 of p-values < \code{colorCutoff} for T1 test, or > beta for star test, to edges in plots;
#'             3 to add red dots at multifurcations;
#'             4 to suppress red dots but show all edge p-values
#' @param refTree a reference tree, as a phylo object, to compute RF distance to 
#' all Tk (for simulation studies); will be unrooted   
#' @param omit \code{FALSE} to treat unresolved quartets on gene trees as 1/3 of each resolution;
#' \code{TRUE} to discard unresolved quartet data; ignored if genedata given as quartet table
#' @param lambda parameter for power-divergence statistic (e.g., 0 for likelihood
#' ratio statistic, 1 for Chi-squared statistic)
#' @param method "MLest", "conservative", or "bootstrap"
#' @param smallsample	"precomputed" or "bootstrap", method of obtaining p-value
#' when sample is small (<30)
#' @param smallcounts	"precomputed" or "bootstrap", method of obtaining p-value
#' when some (but not all) counts are small
#' @param bootstraps	number of samples for bootstrapping
#'
#' @returns An invisible list with 8 elements:
#' \describe{
#' \item{$treeList}{A list corresponding to edge contraction steps.}
#' \item{$indexEarly}{The index in \code{$treeList} of the first tree such that
#' the combined tree T1 is above \code{alpha}.}
#' \item{$indexLate}{The index in \code{$treeList} of the first tree such that
#' the combined tree T1 of it and all subsequent ones is above \code{alpha}.}
#' \item{$alpha}{The value of \code{alpha} determining the previous 2 entries.}
#' \item{$beta}{The value of \code{beta} for contracting edges using the 
#' combined star test.}
#' \item{$qType}{The quartet type mode used.}
#' \item{$testCorrection}{The p-value combination method used.}
#' \item{$pTable}{A table of quartet p-values for relevant
#' quartet tests.}
#' }
#' 
#' In \code{$treeList}, the
#' first is the input tree,
#' the second the result of contracting edges judged unresolved
#' by the combined star test,
#' and subsequent ones the result of collapsing edges with the current
#' smallest combined T1 test p-value, with the final one a star tree.
#' Each element in this list itself a
#' list of \describe{
#' \item{$tree}{A tree, with all edge lengths 1 and a node labels
#' of -log of the parent edge's combined p (star for the input tree, T1 for others).}
#' \item{$logpTree}{-log base 10 of the combined T1 p-value of the tree using \code{testCorrection}.}
#' \item{$logpT1Edge}{-log base 10 of the combined T1 p-value of the edge(s)
#'  that were contracted from previous tree (\code{NA}
#' for input tree and second tree).}
#' \item{$numEdgeCon}{The number of edges that have been contracted from the input tree}
#' \item{$RF}{The RF distance between this tree and \code{refTree}, if \code{refTree} is not NULL.}
#' }
#'
#' @examples
#' data(ASTRALtreeLeopardusLescroart)
#' Atree=read.tree(text=ASTRALtreeLeopardusLescroart)
#' data(tableLeopardusLescroart)
#' ECToBlob(tableLeopardusLescroart,Atree,"bi","Bon")
#'
#'
#' @export
#'
ECToBlob = function(genedata,
                    tree,
                    qType = "mul",
                    testCorrection = "Bon",
                    alpha = NA,
                    beta = .8,
                    colorCutoff = .05,
                    plot = 3,
                    refTree = NULL,
                    omit=FALSE,
                    lambda = 0,
                    method = "MLest",
                    smallsample = "precomputed",
                    smallcounts = "precomputed",
                    bootstraps = 10 ^ 4
                    )
{
  if (qType == "bi") {
    out = ECToBlob_bi(
      genedata = genedata,
      tree = tree,
      testCorrection = testCorrection,
      alpha = alpha,
      beta = beta,
      colorCutoff = colorCutoff,
      plot = plot,
      refTree = refTree,
      omit=omit,
      lambda = lambda,
      method = method,
      smallsample = smallsample,
      smallcounts = smallcounts,
      bootstraps = bootstraps
    )
    return(out)
  }
  else {
    if (qType == "quad") {
      out = ECToBlob_quad(
        genedata = genedata,
        tree = tree,
        testCorrection = testCorrection,
        alpha = alpha,
        beta = beta,
        colorCutoff = colorCutoff,
        plot = plot,
        refTree = refTree,
        omit=omit,
        lambda = lambda,
        method = method,
        smallsample = smallsample,
        smallcounts = smallcounts,
        bootstraps = bootstraps
      )
      return(out)
    }
    else {
      if (qType == "mul") {
        out = ECToBlob_mul(
          genedata = genedata,
          tree = tree,
          testCorrection = testCorrection,
          alpha = alpha,
          beta = beta,
          colorCutoff = colorCutoff,
          plot = plot,
          refTree = refTree,
          omit=omit,
          lambda = lambda,
          method = method,
          smallsample = smallsample,
          smallcounts = smallcounts,
          bootstraps = bootstraps
        )
        return(out)
      }
      else
        stop("Invalid argument: qType")
    }
  }
}

####################################################################
#' Edge Contraction to infer a Tree of Blobs - bipartition mode
#'
#' Infer the tree of blobs using edge contraction and bipartitions
#'
#' Given quartet gene tree count data from a network and a starting tree that is
#' a resolution of its tree of blobs, contract edges for which there is poor support
#' due to, first, star-like quartet signals and then successively the most extreme
#' non-tree-like signals.
#'
#' p-values for star and T1 tests for individual quartets displayed on the tree are
#' combined by one of 4 test correction methods suitable for non-independent tests. For star test,
#' only p-values for quartets defining each edge are combined. For T1 test,
#' p-values for quartets defining paths containing the edge (i.e., from the edge bipartitions)
#' are combined.
#'
#' The input tree should be computed by ASTRAL, or another heuristic known to
#' produce a resolution of the tree of blobs in ideal circumstances.
#'
#' Plot levels allow for several graphical depictions of the edge contract process,
#' from the starting tree to a final star tree.
#'
#' For testing with simulated data, a reference tree (e.g., the true Tree of Blobs) 
#' may be provided and the RF distance between it and each tree will be shown.
#'
#' @references
#' \insertRef{ECTo2026}{MSCquartets}
#'
#' \insertRef{MAR19}{MSCquartets}
#'
#' \insertRef{LiuXie2020}{MSCquartets}
#'
#' \insertRef{Chen2022}{MSCquartets}
#' 
#' @param genedata one of:
#' \enumerate{
#' \item a character string giving the name of a file containing a list of newick 
#' gene trees (see argument \code{omit}),
#' \item a multiphylo object containing a collection of gene trees (see argument \code{omit}),
#' \item a table of resolved quartet counts  from a collection of
#' gene trees, or that table with additional columns of quartet p-values. If p_star
#' or both p_T1-D and qindex-D are present (such as from a previous call to this function) 
#' they are not recomputed.}
#' @param tree a tree as a phylo object, or a character string giving the name of 
#' a file with a newick tree, assumed to be a resolution of the true
#' tree of blobs
#' @param testCorrection a method for combining possibly dependent p-values to 
#' get a single p-value, one of "Cauchy", "Bon", "BBC", "CBC"
#' @param alpha test level for choosing inferred tree of blobs
#' @param beta test level for collapsing edges showing no resolution
#' @param colorCutoff level for changing color of edge p-values for T1 test
#' @param plot 0 for no plot;
#'             1 to plot all trees as edges are contracted, along with plot of tree p-values;
#'             2 to add -log_10 of p-values < \code{colorCutoff} for T1 test, or > beta for star test, to edges in plots;
#'             3 to add red dots at multifurcations;
#'             4 to suppress red dots but show all edge p-values.
#' @param refTree a reference tree, as a phylo object, to compute RF distance to 
#' all Tk (for simulation studies); will be unrooted
#' @param omit \code{FALSE} to treat unresolved quartets on gene trees as 1/3 of each resolution;
#' \code{TRUE} to discard unresolved quartet data; ignored if genedata given as quartet table
#' @param method "MLest", "conservative", or "bootstrap"
#' @param lambda parameter for power-divergence statistic (e.g., 0 for likelihood
#' ratio statistic, 1 for Chi-squared statistic)
#'
#' @param smallsample	"precomputed" or "bootstrap", method of obtaining p-value
#' when sample is small (<30)
#' @param smallcounts	"precomputed" or "bootstrap", method of obtaining p-value
#' when some (but not all) counts are small
#' @param bootstraps	number of samples for bootstrapping
#'
#' @returns An invisible list with 8 elements:
#' \describe{
#' \item{$treeList}{A list corresponding to edge contraction steps.}
#' \item{$indexEarly}{The index in \code{$treeList} of the first tree such that
#' the combined tree T1 is above \code{alpha}.}
#' \item{$indexLate}{The index in \code{$treeList} of the first tree such that
#' the combined tree T1 of it and all subsequent ones is above \code{alpha}.}
#' \item{$alpha}{The value of \code{alpha} determining the previous 2 entries.}
#' \item{$beta}{The value of \code{beta} for contracting edges using the 
#' combined star test.}
#' \item{$qType}{The quartet type mode used.}
#' \item{$testCorrection}{The p-value combination method used.}
#' \item{$pTable}{A table of quartet p-values for relevant
#' quartet tests.}
#' }
#' 
#' In \code{$treeList}, the
#' first is the input tree,
#' the second the result of contracting edges judged unresolved
#' by the combined star test,
#' and subsequent ones the result of collapsing edges with the current
#' smallest combined T1 test p-value, with the final one a star tree.
#' Each element in this list itself a
#' list of \describe{
#' \item{$tree}{A tree, with all edge lengths 1 and a node labels
#' of -log of the parent edge's combined p (star for the input tree, T1 for others).}
#' \item{$logpTree}{-log base 10 of the combined T1 p-value of the tree using \code{testCorrection}.}
#' \item{$logpT1Edge}{-log base 10 of the combined T1 p-value of the edge(s)
#'  that were contracted from previous tree (\code{NA}
#' for input tree and second tree).}
#' \item{$numEdgeCon}{The number of edges that have been contracted from the input tree}
#' \item{$RF}{The RF distance between this tree and \code{refTree}, if \code{refTree} is not NULL.}
#' }
#'
#' @examples
#' data(ASTRALtreeLeopardusLescroart)
#' Atree=read.tree(text=ASTRALtreeLeopardusLescroart)
#' data(tableLeopardusLescroart)
#' ECToBlob_bi(tableLeopardusLescroart,Atree,"Bon")
#'
#' @export
#'
ECToBlob_bi = function(genedata,
                       tree,
                       testCorrection = "Bon",
                       alpha = NA,
                       beta = .8,
                       colorCutoff = .05,
                       plot = 3,
                       refTree = NULL,
                       omit=FALSE,
                       lambda = 0,
                       method = "MLest",
                       smallsample = "precomputed",
                       smallcounts = "precomputed",
                       bootstraps = 10 ^ 4
                       )
{
  if (!(testCorrection %in% c("Cauchy", "Bon", "BBC", "CBC")))
    stop("Invalid argument: testCorrection")
  if (!(plot %in% c(0,1,2,3,4)))
    stop("Invalid argument: plot")
  
  #get data on quartets
  if (!(any(c("matrix","character","multiPhylo") %in% class(genedata)))) 
            stop("Data must be supplied as an object of type multiPhylo, character, or matrix.")
  if (!("matrix" %in% class(genedata))) {
  if ("character" %in% class(genedata)) {
      genetrees <- read.tree(file=genedata) #read gene trees
      message(paste("Read", length(genetrees), "gene trees from file."))
  } else { genetrees = genedata} 
    taxanames = genetrees[[1]]$tip.label   # ... get taxa names from first tree
    taxanames = sort(taxanames)
    if (length(taxanames) <= 25) {
      namelist = paste0(taxanames, collapse = ", ")
    } else {
      namelist = paste0(paste0(taxanames[1:25], collapse = ", "),
                        ",...(see output table for full list)")
    }
    message("Analyzing ", length(taxanames), " taxa: ", namelist)
    genedata = quartetTable(genetrees)   # tally quartets on gene trees
    genedata = quartetTableResolved(genedata, omit)   # treat unresolved quartets
  }

  # prepare tree
  if ("character" %in% class(tree)) { tree <- read.tree(file=tree)}
  tree = unroot(tree)
  tree = di2multi(tree) #makes sure tree is unrooted and no zero length branches
  ntaxa = length(tree$tip.label) #number of taxa
  tree$edge.length = rep(1, dim(tree$edge)[1]) #make tree metric, with all edges length 1
  internalEdges = which(tree$edge[, 2] > length(tree$tip.label)) #numbers of internal edges
  if (length(internalEdges) == 0)
    stop("Tree has no internal edges.")
  treeTaxa=sort(tree$tip.label)
  dataTaxa=sort(colnames(genedata)[1:ntaxa])
  
  if (!all.equal(treeTaxa,dataTaxa))
    stop("Taxa on tree do not match those in genedata")
  numT0edges = dim(tree$edge)[[1]]

  out = c()# create empty output
  out$qType = "bi"
  out$testCorrection = testCorrection
  
  # apply T1-D and star test for all sets of 4 taxa
  pTable = quartetECtestInd(
    genedata,
    tree,
    lambda = lambda,
    method = method,
    smallsample = smallsample,
    smallcounts = smallcounts,
    bootstraps = bootstraps
  )
  out$pTable = pTable

  D = cophenetic(tree) #wipe out p_T1-D for quartets not resolved
  D = D[order(rownames(D)), order(colnames(D))]
  for (m in 1:dim(pTable)[1]) {
    qnames = which(pTable[m, 1:ntaxa] == 1)
    a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
    b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
    c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
    if (min(c(a, b, c)) == max(c(a, b, c)))
      pTable[m, "p_T1-D"] = NA
  }

  #prepare for indexing table entries
  {
    taxanames = colnames(pTable)[1:ntaxa] #names of taxa, as ordered in table
    np1 = ntaxa + 1 #number of taxa + 1
    #store binomial coefficients for efficiency in computing indices in table to access sets of 4 taxa
    C = matrix(nrow = ntaxa, ncol = 4) #store binomial coefficients for efficiency in computing indices in table to access sets of 4 taxa
    for (i in 0:(ntaxa - 1)) {
      for (j in 1:4) {
        C[i + 1, j] = choose(i, j)
      }
    }
    Cn4 = choose(ntaxa, 4) # need one extra value
  }

  #Compute combined p-values for star test for all edges
  {
    edgep = matrix(0, length(internalEdges), 2) #create array for combined p-values
    colnames(edgep) = c("edge", "pcomb")

    for (j in 1:length(internalEdges)) {
      #Combine p_star for all quartets defining this edge on original tree.
      pNode = tree$edge[internalEdges[j], 1] #find taxon groups around parent and child node of edge
      cNode = tree$edge[internalEdges[j], 2]
      pGroups = nodeGroups(tree, pNode)
      cGroups = nodeGroups(tree, cNode)

      #remove group above child (depends on nodeGroups returning non-descendants last in list)
      cGroups = cGroups[-length(cGroups)]
      cUnion = unlist(cGroups)

      for (i in 1:length(pGroups)) {
        #remove group below parent
        if (all(cUnion %in% pGroups[[i]]))
          remove = i
      }
      pGroups = pGroups[-remove]

      # compute combined p values for star
      pvecstar = c() #storage for p_star values

      #run through all choices of 2 groups at child and 2 at parent
      for (group1 in 2:length(pGroups)) {
        for (group2 in 1:(group1 - 1)) {
          for (group3 in 2:length(cGroups)) {
            for (group4 in 1:(group3 - 1)) {
              #run through all choices of 1 taxon per group
              for (at in pGroups[[group1]]) {
                for (bt in pGroups[[group2]]) {
                  for (ct in cGroups[[group3]]) {
                    for (dt in cGroups[[group4]]) {
                      #find a,b,c,d entry in table and get p_=star
                      positions = sort(match(tree$tip.label[c(at, bt, ct, dt)], taxanames)) # find table columns for these taxa
                      index = Cn4 - C[np1 - positions[1], 4] - C[np1 - positions[2], 3] - C[np1 - positions[3], 2] - C[np1 - positions[4], 1] #row index for this quartet
                      pstar = unname(pTable[index, "p_star"])
                      pvecstar = c(pvecstar, pstar)
                    }
                  }
                }
              }
            }
          }
        }
      }
      edgep[j, ] = c(internalEdges[j], combineP(pvecstar, testCorrection))
    }
    edgep = rbind(edgep) #make sure result is a matrix if only 1 row
  }

  #Put edge star p-value -logs in as node labels on T_0
  tree$node.label = rep("", tree$Nnode)
  tree$node.label[tree$edge[edgep[, "edge"], 2] - ntaxa] = -log(edgep[, "pcomb"], 10)

  #Put initial tree T_0 in list, with special treatment if reference Tree supplied
  pTree = combineP(na.omit(pTable[, "p_T1-D"]), testCorrection) # compute combined values for tree T_0
  pTreeList = c(pTree) # save in list
  indexTreeList = c(0)
  if (is.null(refTree))
    trees = list(list(
      "tree" = ape::write.tree(tree),
      "logpTree" = -log(pTree, 10),
      "logpT1Edge" = NA,
      "numEdgeCon" = 0
    ))
  else{
    refTree=unroot(refTree)
    RF = as.numeric(ape::dist.topo(refTree, tree))
    trees = list(
      list(
        "tree" = ape::write.tree(tree),
        "logpTree" = -log(pTree, 10),
        "logpT1Edge" = NA,
        "numEdgeCon" = 0,
        "RF" = RF
      )
    )# put it in list
    RFlist = c(RF)# and save RF for plotting
  }

  if (plot > 0) {
    #Plot initial tree
    line1 = ""
    if (!is.null(refTree))
      line1 = paste0("RF to reference Tree = ", RF, "; \n")

    plot(
      tree,
      type = "unrooted",
      cex.sub = .7,
      main = paste0("T0 (Input tree): ", testCorrection, ", bi"),
      sub = paste0(
        line1,
        "-log_10(p) for combined T1 over tree = ",
        formatC(-log(pTree, 10), format = 'e', digits = 3),
        "; \n Edge labels are -log_10(p) for combined star tests for that edge."
      )
    )
    if (plot %in% c(2, 3)) {
      badedgeindex = which(edgep[, "pcomb"] > beta)
      if (length(badedgeindex) > 0) {
        badedges = edgep[badedgeindex, "edge"]
        bgColor = rep("tomato1", length(badedges))
        ape::edgelabels(
          formatC(
            -log(edgep[badedgeindex, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = badedges,
          bg = bgColor,
          cex = .5
        )
      }
    }
    if (plot == 3) {
      polyList = c() #find multifurcations, starting at "top" in unrooted tree
      if (length(Descendants(tree, ntaxa + 1, 'children')) > 3)
        polyList = c(polyList, ntaxa + 1)
      if (tree$Nnode >= 2) {
        #loop over other internal nodes
        for (i in (ntaxa + 2):(ntaxa + tree$Nnode)) {
          if (length(Descendants(tree, i, 'children')) > 2)
            polyList = c(polyList, i)
        }
      }
      #add colored dots for these nodes
      lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
      XX <- lastPP$xx[polyList]
      YY <- lastPP$yy[polyList]
      points(XX,
             YY,
             pch = 16,
             col = "red",
             cex = 1)
    }
    if (plot == 4) {
      edgeind = which(edgep[, "pcomb"] != Inf)
      if (length(edgeind) > 0) {
        bgColor = rep("lightblue", dim(edgep)[[1]])
        bgColor[which(edgep[edgeind, "pcomb"] > beta)] = "tomato1"
        ape::edgelabels(
          formatC(
            -log(edgep[edgeind, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = edgep[edgeind, "edge"],
          bg = bgColor,
          cex = .5
        )
      }
    }
  }

  #Contract edges based on star test
  message('Contracting edges with star test at level beta = ',beta)
  edgepRow = which(edgep[, "pcomb"] > beta) #determine entries for edges to contract for star test
  edges2contract = edgep[edgepRow, "edge"] # and the edge numbers
  num2contract = length(edges2contract)
  tree$edge.length[edges2contract] = 0 #contract those edges
  
 
  D = cophenetic(tree) #wipe out p_T1-D for quartets no longer resolved
  D = D[order(rownames(D)), order(colnames(D))]
  for (m in 1:dim(pTable)[1]) {
    qnames = which(pTable[m, 1:ntaxa] == 1)
    a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
    b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
    c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
    if (min(c(a, b, c)) == max(c(a, b, c)))
      pTable[m, "p_T1-D"] = NA
  }
  
  # Computation of combined edge p for T1 test, once and for all, using T0
  internalEdges = which(tree$edge[, 2] > ntaxa) #numbers of internal edges
  edgep = matrix(0, length(internalEdges), 2) #create array for combined p-values
  colnames(edgep) = c("edge", "pcomb")
  for (ie in 1:length(internalEdges)) {
    if (ie %in% edges2contract) pvecT1=c(1) #ignore edges just contracted
    else { 
    pvecT1=c()
    cNode = tree$edge[internalEdges[ie], 2] # child node of edge
    cUnion = Descendants(tree, cNode)[[1]] # descendant taxa
    pUnion = setdiff(1:ntaxa, cUnion) # and other taxa
    for (i in 2:length(pUnion)) {
      for (j in 1:(i - 1)) {
        for (k in 2:(length(cUnion))) {
          for (l in 1:(k - 1)) {
            at = pUnion[[i]]
            bt = pUnion[[j]]
            ct = cUnion[[k]]
            dt = cUnion[[l]]
            positions = sort(match(tree$tip.label[c(at, bt, ct, dt)], taxanames)) # find table columns for these taxa
            index = Cn4 - C[np1 - positions[1], 4] - C[np1 - positions[2], 3] - C[np1 - positions[3], 2] - C[np1 - positions[4], 1] #row index for this quartet
            pvecT1 = c(pvecT1, unname(pTable[index, "p_T1-D"]))
          }
        }
      }
    }
    }
    edgep[ie, ] = c(internalEdges[ie], combineP(pvecT1, testCorrection))
  }
  edgep = rbind(edgep) #make sure result is a matrix if only 1 row

  #Put edge T1 p-value -logs in as node labels on Tk
  tree$node.label = rep("", tree$Nnode)
  tree$node.label[tree$edge[edgep[, "edge"], 2] - ntaxa] = -log(edgep[, "pcomb"], 10)
  treem = di2multi(tree) # tree for putting in output, with multifurcations
  numEdges = dim(treem$edge)[[1]]

  #Put star-contracted tree Tk in list, with special treatment if reference Tree supplied
  pTree = combineP(na.omit(pTable[, "p_T1-D"]), testCorrection) # compute combined values for tree T_0
  pTreeList = c(pTreeList, pTree) # save in list
  indexTreeList = c(indexTreeList, num2contract)
  if (is.null(refTree))
    trees = append(trees, list(
      list(
        "tree" = ape::write.tree(treem),
        "logpTree" = -log(pTree, 10),
        "logpT1Edge" = NA,
        "numEdgeCon" = num2contract
      )
    ))
  else{
    RF = as.numeric(ape::dist.topo(refTree, treem))
    trees = append(trees, list(
      list(
        "tree" = ape::write.tree(treem),
        "logpTree" = -log(pTree, 10),
        "logpT1Edge" = NA,
        "numEdgeCon" = num2contract,
        "RF" = RF
      )
    ))# put it in list
    RFlist = c(RFlist,RF)# and save RF for plotting
  }

  if (plot > 0) {
    #Plot this tree
    line1 = ""
    if (!is.null(refTree))
      line1 = paste0("RF to reference Tree = ", RF, "; \n")

    plot(
      tree,
      type = "unrooted",
      cex.sub = .7,
      main = paste0("T", num2contract, ": ", testCorrection, ", bi"),
      sub = paste0(
        line1,
        "-log_10(p) for combined T1 over tree = ",
        formatC(-log(pTree, 10), format = 'e', digits = 3),
        "; \n Edge labels are -log_10(p) for combined T1 test for that edge;\n",
        num2contract,
        " edges contracted by combined star test."
      )
    )
    if (plot %in% c(2, 3)) {
      badedgeindex = which(edgep[, "pcomb"] < colorCutoff)
      if (length(badedgeindex) > 0) {
        badedges = edgep[badedgeindex, "edge"]
        bgColor = rep("tomato1", length(badedges))
        ape::edgelabels(
          formatC(
            -log(edgep[badedgeindex, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = badedges,
          bg = bgColor,
          cex = .5
        )
      }
    }
    if (plot == 3) {
      polyList = c() #find multifurcations, starting at "top" in unrooted tree
      if (length(Descendants(treem, ntaxa + 1, 'children')) > 3)
        polyList = c(polyList, ntaxa + 1)
      if (tree$Nnode >= 2) {
        #loop over other internal nodes
        for (i in (ntaxa + 2):(ntaxa + tree$Nnode)) {
          if (length(Descendants(treem, i, 'children')) > 2)
            polyList = c(polyList, i)
        }
      }
      #add colored dots for these nodes
      lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
      XX <- lastPP$xx[polyList]
      YY <- lastPP$yy[polyList]
      points(XX,
             YY,
             pch = 16,
             col = "red",
             cex = 1)
    }
    if (plot == 4) {
      edgeind = which(edgep[, "pcomb"] != Inf)
      if (length(edgeind) > 0) {
        bgColor = rep("lightblue", dim(edgep)[[1]])
        bgColor[which(edgep[edgeind, "pcomb"] < colorCutoff)] = "tomato1"
        ape::edgelabels(
          formatC(
            -log(edgep[edgeind, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = edgep[edgeind, "edge"],
          bg = bgColor,
          cex = .5
        )
      }
    }
  }

   while (numEdges > ntaxa)
    #loop with successive contracting steps
  {
    minp = min(edgep[, "pcomb"],na.rm = TRUE) # p for edges to contract
    contract = which(edgep[, "pcomb"] == minp) # find all minimal edge p
    tree$edge.length[edgep[contract, 1]] = 0 # and set their lengths to 0
    edgep[contract, 2] = Inf # and never consider them again
    treem = di2multi(tree) #note this retains correct node label for edge p
    numEdges = dim(treem$edge)[[1]]
    
    D = cophenetic(treem) #wipe out p_T1-D for quartets no longer resolved
    D = D[order(rownames(D)), order(colnames(D))]
    for (m in 1:dim(pTable)[1]) {
      qnames = which(pTable[m, 1:ntaxa] == 1)
      a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
      b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
      c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
      if (min(c(a, b, c)) == max(c(a, b, c)))
        pTable[m, "p_T1-D"] = NA
    }

    countContracted = numT0edges - numEdges #total number of edges contracted
    pTree = combineP(na.omit(pTable[, "p_T1-D"]), testCorrection) # compute new combined values for tree
    pTreeList = c(pTreeList, pTree) # save new tree p
    indexTreeList = c(indexTreeList, countContracted)


    if (is.null(refTree))
      trees = append(trees, list(
        list(
          "tree" = ape::write.tree(treem),
          "logpTree" = -log(pTree, 10),
          "logpT1Edge" = -log(minp, 10),
          "numEdgeCon" = countContracted
        )
      ))

    else{
      RF = as.numeric(ape::dist.topo(refTree, treem))
      trees = append(trees, list(
        list(
          "tree" = ape::write.tree(treem),
          "logpTree" = -log(pTree, 10),
          "logpT1Edge" = -log(minp, 10),
          "numEdgeCon" = countContracted,
          "RF" = RF
        )
      ))
      RFlist = c(RFlist, RF)# and save RF for plotting
    }

    if (plot > 0) {
      line1 = ""
      if (!is.null(refTree))
        line1 = paste0("RF to reference Tree = ", RF, "; \n")
      plot(
        tree,
        type = "unrooted",
        cex.sub = .7,
        main = paste0("T", countContracted, ": ", testCorrection, ", bi"),
        sub = paste0(
          line1,
          "-log_10(p) for combined T1 over tree = ",
          formatC(-log(pTree, 10), format = 'e', digits = 3),
          "; \n Edge labels are -log_10(p) from combined edge T1;\n",
          length(contract),
          " more edge(s) contracted by combined T1 test."
        )
      )

      if (plot > 1)
      {
        if (plot %in% c(2, 3)) {
          badedgeindex = which(edgep[, "pcomb"] < colorCutoff)
          if (length(badedgeindex > 0))
          {
            badedges = edgep[badedgeindex, "edge"]
            bgColor = rep("tomato1", length(badedges))
            ape::edgelabels(
              formatC(
                -log(edgep[badedgeindex, "pcomb"], 10),
                format = 'e',
                digits = 0
              ),
              edge = badedges,
              bg = bgColor,
              cex = .5
            )
          }
        }
        if (plot == 3)
        {
          contractedEdges = which(tree$edge.length == 0)
          polyList = tree$edge[contractedEdges, 1]
          #colored dots for these nodes
          lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
          internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
          XX <- lastPP$xx[polyList]
          YY <- lastPP$yy[polyList]
          points(XX,
                 YY,
                 pch = 16,
                 col = "red",
                 cex = 1)
        }
        if (plot == 4) {
          edgeind = which(edgep[, "pcomb"] != Inf)
          if (length(edgeind) > 0) {
            bgColor = rep("lightblue", dim(edgep)[[1]])
            bgColor[which(edgep[edgeind, "pcomb"] < colorCutoff)] = "tomato1"
            ape::edgelabels(
              formatC(
                -log(edgep[edgeind, "pcomb"], 10),
                format = 'e',
                digits = 0
              ),
              edge = edgep[edgeind, "edge"],
              bg = bgColor,
              cex = .5
            )
          }
        }
      }
    }
    }
  
  out$treeList = trees

  if (!is.na(alpha)) {
    small = which(pTreeList < alpha)
    if (length(small) == 0)
      small = 0
    out$indexLate = max(small) + 1
    out$indexEarly = which(pTreeList > alpha)[1]
  } else {
    out$indexLate = NA
    out$indexEarly = NA
  }

  #final plot, of p-values for tree, and RF
  if (plot > 0) {
    y = -log(pTreeList, 10)
    yfinpos = which(is.finite(y))
    yinfpos = setdiff(1:length(y), yfinpos)
    yinf = max(y[yfinpos]) * 1.5
    plot(
      indexTreeList[yinfpos],
      rep(yinf, length(yinfpos)),
      ylim = c(0, yinf),
      xlim = c(-.5, max(indexTreeList) + .5),
      main = paste0("Tree p-values; ", testCorrection, ", bi"),
      xlab = "Number of contracted edges",
      ylab = "-log_10(p) (solid red = Inf)",
      col = "red",
      pch = 19
    )
    points(indexTreeList[yfinpos], y[yfinpos], col = "blue")


    if (!is.null(refTree))
      plot(
        indexTreeList,
        RFlist,
        main = paste0("RFdist to supplied ToB; ", testCorrection, ", bi"),
        xlab = "Tree",
        ylab = "RFdist",
        col = "blue"
      )
  }

  out$alpha=alpha
  out$beta=beta
  
  out = out[c("treeList",
              "indexEarly",
              "indexLate",
              "alpha",
              "beta",
              "qType",
              "testCorrection",
              "pTable")]

  invisible(out)
}

####################################################################
#' Edge Contraction to infer a Tree of Blobs - quadripartition mode
#'
#' Infer the tree of blobs using edge contraction and quadripartitions
#'
#' Given quartet gene tree count data from a network and a starting tree that is
#' a resolution of its tree of blobs, contract edges for which there is poor support
#' due to, first, star-like quartet signals and then successively the most extreme
#' non-tree-like signals.
#'
#' p-values for star and T1 tests for individual quartets displayed on the tree are
#' combined by one of 4 test correction methods suitable for non-independent tests. For both
#' star and T1 test, only p-values for quartets defining each edge are combined.
#' (For T1, this means defining the edge after any star edge contractions have been done.)
#'
#' The input tree should be computed by ASTRAL, or another heuristic known to
#' produce a resolution of the tree of blobs in ideal circumstances.
#'
#' Plot levels allow for several graphical depictions of the edge contract process,
#' from the starting tree to a final star tree.
#'
#' For testing with simulated data, a reference tree (e.g., the true Tree of Blobs) 
#' may be provided and the RF distance between it and each tree will be shown.
#'
#' @references
#' \insertRef{ECTo2026}{MSCquartets}
#'
#' \insertRef{MAR19}{MSCquartets}
#'
#' \insertRef{LiuXie2020}{MSCquartets}
#'
#' \insertRef{Chen2022}{MSCquartets}
#' 
#' @param genedata one of:
#' \enumerate{
#' \item a character string giving the name of a file containing a list of newick 
#' gene trees (see argument \code{omit}),
#' \item a multiphylo object containing a collection of gene trees (see argument \code{omit}),
#' \item a table of resolved quartet counts  from a collection of
#' gene trees, or that table with additional columns of quartet p-values. If p_star
#' or both p_T1-D and qindex-D are present (such as from a previous call to this function) 
#' they are not recomputed.}
#' @param tree a tree as a phylo object, or a character string giving the name of 
#' a file with a newick tree, assumed to be a resolution of the true
#' tree of blobs
#' @param testCorrection a method for combining possibly dependent p-values to 
#' get a single p-value, one of "Cauchy", "Bon", "BBC", "CBC"
#' @param alpha test level for choosing inferred tree of blobs
#' @param beta test level for collapsing edges showing no resolution
#' @param colorCutoff level for changing color of edge p-values for T1 test
#' @param plot 0 for no plot;
#'             1 to plot all trees as edges are contracted, along with plot of tree p-values;
#'             2 to add -log_10 of p-values < \code{colorCutoff} for T1 test, or > beta for star test, to edges in plots;
#'             3 to add red dots at multifurcations;
#'             4 to suppress red dots but show all edge p-values
#' @param refTree a reference tree, as a phylo object, to compute RF distance to 
#' all Tk (for simulation studies); will be unrooted
#' @param omit \code{FALSE} to treat unresolved quartets on gene trees as 1/3 of each resolution;
#' \code{TRUE} to discard unresolved quartet data; ignored if genedata given as quartet table
#' @param lambda parameter for power-divergence statistic (e.g., 0 for likelihood
#' ratio statistic, 1 for Chi-squared statistic)
#' @param method "MLest", "conservative", or "bootstrap"
#' @param smallsample	"precomputed" or "bootstrap", method of obtaining p-value
#' when sample is small (<30)
#' @param smallcounts	"precomputed" or "bootstrap", method of obtaining p-value
#' when some (but not all) counts are small
#' @param bootstraps	number of samples for bootstrapping
#'
#' @returns An invisible list with 8 elements:
#' \describe{
#' \item{$treeList}{A list corresponding to edge contraction steps.}
#' \item{$indexEarly}{The index in \code{$treeList} of the first tree such that
#' the combined tree T1 is above \code{alpha}.}
#' \item{$indexLate}{The index in \code{$treeList} of the first tree such that
#' the combined tree T1 of it and all subsequent ones is above \code{alpha}.}
#' \item{$alpha}{The value of \code{alpha} determining the previous 2 entries.}
#' \item{$beta}{The value of \code{beta} for contracting edges using the 
#' combined star test.}
#' \item{$qType}{The quartet type mode used.}
#' \item{$testCorrection}{The p-value combination method used.}
#' \item{$pTable}{A table of quartet p-values for relevant
#' quartet tests.}
#' }
#' 
#' In \code{$treeList}, the
#' first is the input tree,
#' the second the result of contracting edges judged unresolved
#' by the combined star test,
#' and subsequent ones the result of collapsing edges with the current
#' smallest combined T1 test p-value, with the final one a star tree.
#' Each element in this list itself a
#' list of \describe{
#' \item{$tree}{A tree, with all edge lengths 1 and a node labels
#' of -log of the parent edge's combined p (star for the input tree, T1 for others).}
#' \item{$logpTree}{-log base 10 of the combined T1 p-value of the tree using \code{testCorrection}.}
#' \item{$logpT1Edge}{-log base 10 of the combined T1 p-value of the edge(s)
#'  that were contracted from previous tree (\code{NA}
#' for input tree and second tree).}
#' \item{$numEdgeCon}{The number of edges that have been contracted from the input tree}
#' \item{$RF}{The RF distance between this tree and \code{refTree}, if \code{refTree} is not NULL.}
#' }
#'
#' @examples
#' data(ASTRALtreeLeopardusLescroart)
#' Atree=read.tree(text=ASTRALtreeLeopardusLescroart)
#' data(tableLeopardusLescroart)
#' ECToBlob_quad(tableLeopardusLescroart,Atree,"Bon")
#'
#' @export
#'
ECToBlob_quad = function(genedata,
                         tree,
                         testCorrection = "Bon",
                         alpha = NA,
                         beta = .8,
                         colorCutoff = .05,
                         plot = 3,
                         refTree = NULL, 
                         omit=FALSE,
                         lambda = 0,
                         method = "MLest",
                         smallsample = "precomputed",
                         smallcounts = "precomputed",
                         bootstraps = 10 ^ 4
                         )
{
  if (!(testCorrection %in% c("Cauchy", "Bon", "BBC", "CBC")))
    stop("Invalid argument: testCorrection")
  if (!(plot %in% c(0,1,2,3,4)))
    stop("Invalid argument: plot")
  
  #get data on quartets
  if (!(any(c("matrix","character","multiPhylo") %in% class(genedata)))) 
    stop("Data must be supplied as an object of type multiPhylo, character, or matrix.")
  if (!("matrix" %in% class(genedata))) {
    if ("character" %in% class(genedata)) {
      genetrees <- read.tree(file=genedata) #read gene trees
      message(paste("Read", length(genetrees), "gene trees from file."))
    } else { genetrees = genedata} 
    taxanames = genetrees[[1]]$tip.label   # ... get taxa names from first tree
    taxanames = sort(taxanames)
    if (length(taxanames) <= 25) {
      namelist = paste0(taxanames, collapse = ", ")
    } else {
      namelist = paste0(paste0(taxanames[1:25], collapse = ", "),
                        ",...(see output table for full list)")
    }
    message("Analyzing ", length(taxanames), " taxa: ", namelist)
    genedata = quartetTable(genetrees)   # tally quartets on gene trees
    genedata = quartetTableResolved(genedata, omit)   # treat unresolved quartets
  }
  
  # prepare tree
  if ("character" %in% class(tree)) { tree <- read.tree(file=tree) } #read input tree 
  tree = unroot(tree)
  tree = di2multi(tree) #makes sure tree is unrooted and no zero length branches
  ntaxa = length(tree$tip.label) #number of taxa
  tree$edge.length = rep(1, dim(tree$edge)[1]) #make tree metric, with all edges length 1
  internalEdges = which(tree$edge[, 2] > length(tree$tip.label)) #numbers of internal edges
  if (length(internalEdges) == 0)
    stop("Tree has no internal edges.")
  treeTaxa=sort(tree$tip.label)
  dataTaxa=sort(colnames(genedata)[1:ntaxa])
  if (!all.equal(treeTaxa,dataTaxa))
    stop("Taxa on tree do not match those in genedata")
  numT0edges = dim(tree$edge)[[1]]
  
  out = c()# create empty output
  out$qType = "quad"
  out$testCorrection = testCorrection

  # apply T1-D and star test for all sets of 4 taxa
  pTable = quartetECtestInd(
    genedata,
    tree,
    lambda = lambda,
    method = method,
    smallsample = smallsample,
    smallcounts = smallcounts,
    bootstraps = bootstraps
  )
  out$pTable = pTable

  D = cophenetic(tree) #wipe out p_T1-D for any quartets that are not resolved
  D = D[order(rownames(D)), order(colnames(D))]
  for (m in 1:dim(pTable)[1]) {
    qnames = which(pTable[m, 1:ntaxa] == 1)
    a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
    b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
    c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
    if (min(c(a, b, c)) == max(c(a, b, c)))
      pTable[m, "p_T1-D"] = NA
  }

  #prepare for indexing table entries
  {
    taxanames = colnames(pTable)[1:ntaxa] #names of taxa, as ordered in table
    np1 = ntaxa + 1 #number of taxa + 1
    #store binomial coefficients for efficiency in computing indices in table to access sets of 4 taxa
    C = matrix(nrow = ntaxa, ncol = 4) #store binomial coefficients for efficiency in computing indices in table to access sets of 4 taxa
    for (i in 0:(ntaxa - 1)) {
      for (j in 1:4) {
        C[i + 1, j] = choose(i, j)
      }
    }
    Cn4 = choose(ntaxa, 4) # need one extra value
  }

  #Compute combined p-values for star and T1 test for all edges
  {
    edgep = matrix(0, length(internalEdges), 3) #create array for combined p-values
    colnames(edgep) = c("edge", "pcombstar", "pcomb")

    for (j in 1:length(internalEdges)) {
      #Combine p_star for all quartets defining this edge on original tree.
      pNode = tree$edge[internalEdges[j], 1] #find taxon groups around parent and child node of edge
      cNode = tree$edge[internalEdges[j], 2]
      pGroups = nodeGroups(tree, pNode)
      cGroups = nodeGroups(tree, cNode)

      #remove group above child (depends on nodeGroups returning non-descendants last in list)
      cGroups = cGroups[-length(cGroups)]
      cUnion = unlist(cGroups)

      for (i in 1:length(pGroups)) {
        #remove group below parent
        if (all(cUnion %in% pGroups[[i]]))
          remove = i
      }
      pGroups = pGroups[-remove]

      # compute combined p values for star and T1
      pvecstar = c() #storage for p_star values
      pvecT1 = c() #storage for p_T1-D values

      #run through all choices of 2 groups at child and 2 at parent
      for (group1 in 2:length(pGroups)) {
        for (group2 in 1:(group1 - 1)) {
          for (group3 in 2:length(cGroups)) {
            for (group4 in 1:(group3 - 1)) {
              #run through all choices of 1 taxon per group
              for (at in pGroups[[group1]]) {
                for (bt in pGroups[[group2]]) {
                  for (ct in cGroups[[group3]]) {
                    for (dt in cGroups[[group4]]) {
                      #find a,b,c,d entry in table and get p_=star
                      positions = sort(match(tree$tip.label[c(at, bt, ct, dt)], taxanames)) # find table columns for these taxa
                      index = Cn4 - C[np1 - positions[1], 4] - C[np1 - positions[2], 3] - C[np1 - positions[3], 2] - C[np1 - positions[4], 1] #row index for this quartet
                      pstar = unname(pTable[index, "p_star"])
                      pvecstar = c(pvecstar, pstar)
                      pT1 = unname(pTable[index, "p_T1-D"])
                      pvecT1 = c(pvecT1, pT1)
                    }
                  }
                }
              }
            }
          }
        }
      }
      edgep[j, ] = c(
        internalEdges[j],
        combineP(pvecstar, testCorrection),
        combineP(pvecT1, testCorrection)
      )
    }
    edgep = rbind(edgep) #make sure result is a matrix if only 1 row
  }



  #Put edge star p-value -logs in as node labels on T_0
  tree$node.label = rep("", tree$Nnode)
  tree$node.label[tree$edge[edgep[, "edge"], 2] - ntaxa] = -log(edgep[, "pcombstar"], 10)

  #Put initial tree T_0 in list, with special treatment if reference Tree supplied
  pTree = combineP(na.omit(pTable[, "p_T1-D"]), testCorrection) # compute combined values for tree T_0
  pTreeList = c(pTree) # save in list
  indexTreeList = c(0)
  if (is.null(refTree))
    trees = list(list(
      "tree" = ape::write.tree(tree),
      "logpTree" = -log(pTree, 10),
      "logpT1Edge" = NA,
      "numEdgeCon" = 0
    ))
  else{
    refTree=unroot(refTree)
    RF = as.numeric(ape::dist.topo(refTree, tree))
    trees = list(
      list(
        "tree" = ape::write.tree(tree),
        "logpTree" = -log(pTree, 10),
        "logpT1Edge" = NA,
        "numEdgeCon" = 0,
        "RF" = RF
      )
    )# put it in list
    RFlist = c(RF)# and save RF for plotting
  }


  if (plot > 0) {
    #Plot initial tree
    line1 = ""
    if (!is.null(refTree))
      line1 = paste0("RF to reference Tree = ", RF, "; \n")

    plot(
      tree,
      type = "unrooted",
      cex.sub = .7,
      main = paste0("T0 (Input tree): ", testCorrection, ", quad"),
      sub = paste0(
        line1,
        "-log_10(p) for combined T1 over tree = ",
        formatC(-log(pTree, 10), format = 'e', digits = 3),
        "; \n Edge labels are -log_10(p) for combined star tests for that edge."
      )
    )
    if (plot %in% c(2, 3)) {
      badedgeindex = which(edgep[, "pcombstar"] > beta)
      if (length(badedgeindex) > 0) {
        badedges = edgep[badedgeindex, "edge"]
        bgColor = rep("tomato1", length(badedges))
        ape::edgelabels(
          formatC(
            -log(edgep[badedgeindex, "pcombstar"], 10),
            format = 'e',
            digits = 0
          ),
          edge = badedges,
          bg = bgColor,
          cex = .5
        )
      }
    }
    if (plot == 3) {
      polyList = c() #find multifurcations, starting at "top" in unrooted tree
      if (length(Descendants(tree, ntaxa + 1, 'children')) > 3)
        polyList = c(polyList, ntaxa + 1)
      if (tree$Nnode >= 2) {
        #loop over other internal nodes
        for (i in (ntaxa + 2):(ntaxa + tree$Nnode)) {
          if (length(Descendants(tree, i, 'children')) > 2)
            polyList = c(polyList, i)
        }
      }
      #add colored dots for these nodes
      lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
      XX <- lastPP$xx[polyList]
      YY <- lastPP$yy[polyList]
      points(XX,
             YY,
             pch = 16,
             col = "red",
             cex = 1)
    }
    if (plot == 4) {
      edgeind = which(edgep[, "pcombstar"] != Inf)
      if (length(edgeind) > 0) {
        bgColor = rep("lightblue", dim(edgep)[[1]])
        bgColor[which(edgep[edgeind, "pcombstar"] > beta)] = "tomato1"
        ape::edgelabels(
          formatC(
            -log(edgep[edgeind, "pcombstar"], 10),
            format = 'e',
            digits = 0
          ),
          edge = edgep[edgeind, "edge"],
          bg = bgColor,
          cex = .5
        )
      }
    }
  }

  #Contract edges based on star test
  message('Contracting edges with star test at level beta = ',beta)
  edgepRow = which(edgep[, "pcombstar"] > beta) #determine entries for edges to contract for star test
  edges2contract = edgep[edgepRow, "edge"] # and the edge numbers
  num2contract = length(edges2contract)
  tree$edge.length[edges2contract] = 0 #contract those edges

  D = cophenetic(tree) #wipe out p_T1-D for quartets no longer resolved
  D = D[order(rownames(D)), order(colnames(D))]
  for (m in 1:dim(pTable)[1]) {
    qnames = which(pTable[m, 1:ntaxa] == 1)
    a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
    b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
    c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
    if (min(c(a, b, c)) == max(c(a, b, c)))
      pTable[m, "p_T1-D"] = NA
  }


  #Put edge T1 p-value -logs in as node labels on Tk
  tree$node.label = rep("", tree$Nnode)
  tree$node.label[tree$edge[edgep[, "edge"], 2] - ntaxa] = -log(edgep[, "pcomb"], 10)


  #Put star-contracted tree Tk in list, with special treatment if reference Tree supplied
  pTree = combineP(na.omit(pTable[, "p_T1-D"]), testCorrection) # compute combined values for tree T_0
  pTreeList = c(pTreeList, pTree) # save in list
  indexTreeList = c(indexTreeList, num2contract)
  treem = di2multi(tree)
  if (is.null(refTree))
    trees = append(trees, list(
      list(
        "tree" = ape::write.tree(treem),
        "logpTree" = -log(pTree, 10),
        "logpT1Edge" = NA,
        "numEdgeCon" = num2contract
      )
    ))
  else{
    RF = as.numeric(ape::dist.topo(refTree, tree))
    trees = append(trees, list(
      list(
        "tree" = ape::write.tree(treem),
        "logpTree" = -log(pTree, 10),
        "logpT1Edge" = NA,
        "numEdgeCon" = num2contract,
        "RF" = RF
      )
    ))# put it in list
    RFlist = c(RFlist,RF)# and save RF for plotting
  }

  if (plot > 0) {
    #Plot this tree
    line1 = ""
    if (!is.null(refTree))
      line1 = paste0("RF to reference Tree = ", RF, "; \n")

    plot(
      treem,
      type = "unrooted",
      cex.sub = .7,
      main = paste0("T", num2contract, ": ", testCorrection, ", quad"),
      sub = paste0(
        line1,
        "-log_10(p) for combined T1 over tree = ",
        formatC(-log(pTree, 10), format = 'e', digits = 3),
        "; \n Edge labels are -log_10(p) for combined T1 test for that edge;\n",
        num2contract,
        " edges contracted from by combined star test."
      )
    )
    if (plot %in% c(2, 3)) {
      badedgeindex = which(edgep[, "pcomb"] < colorCutoff)
      if (length(badedgeindex) > 0) {
        badedges = edgep[badedgeindex, "edge"]
        bgColor = rep("tomato1", length(badedges))
        ape::edgelabels(
          formatC(
            -log(edgep[badedgeindex, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = badedges,
          bg = bgColor,
          cex = .5
        )
      }
    }
    if (plot == 3) {
      polyList = c() #find multifurcations, starting at "top" in unrooted tree
      if (length(Descendants(treem, ntaxa + 1, 'children')) > 3)
        polyList = c(polyList, ntaxa + 1)
      if (treem$Nnode >= 2) {
        #loop over other internal nodes
        for (i in (ntaxa + 2):(ntaxa + tree$Nnode)) {
          if (length(Descendants(treem, i, 'children')) > 2)
            polyList = c(polyList, i)
        }
      }
      #add colored dots for these nodes
      lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
      XX <- lastPP$xx[polyList]
      YY <- lastPP$yy[polyList]
      points(XX,
             YY,
             pch = 16,
             col = "red",
             cex = 1)
    }
    if (plot == 4) {
      edgeind = which(edgep[, "pcomb"] != Inf)
      if (length(edgeind) > 0) {
        bgColor = rep("lightblue", dim(edgep)[[1]])
        bgColor[which(edgep[edgeind, "pcomb"] < colorCutoff)] = "tomato1"
        ape::edgelabels(
          formatC(
            -log(edgep[edgeind, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = edgep[edgeind, "edge"],
          bg = bgColor,
          cex = .5
        )
      }
    }
  }

  numEdges = dim(treem$edge)[[1]]

  while (numEdges > ntaxa)
    #loop with successive contracting steps
  {
    minp = min(edgep[, "pcomb"]) # p for edges to contract
    contract = which(edgep[, "pcomb"] == minp) # find all minimal edge p
    tree$edge.length[edgep[contract, 1]] = 0 # and set their lengths to 0
    edgep[contract, "pcomb"] = Inf # and never consider them again
    treem = di2multi(tree) #note this retains correct node label for edge p
    numEdges = dim(treem$edge)[[1]] #how many edges? 
    
    D = cophenetic(tree) #wipe out p_T1-D for quartets no longer resolved
    D = D[order(rownames(D)), order(colnames(D))]
    for (m in 1:dim(pTable)[1]) {
      qnames = which(pTable[m, 1:ntaxa] == 1)
      a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
      b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
      c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
      if (min(c(a, b, c)) == max(c(a, b, c)))
        pTable[m, "p_T1-D"] = NA
    }
    
    countContracted = numT0edges - numEdges #total number of edges contracted
    pTree = combineP(na.omit(pTable[, "p_T1-D"]), testCorrection) # compute new combined values for tree
    pTreeList = c(pTreeList, pTree) # save new tree p
    indexTreeList = c(indexTreeList, countContracted)

    if (is.null(refTree))
      trees = append(trees, list(
        list(
          "tree" = ape::write.tree(treem),
          "logpTree" = -log(pTree, 10),
          "logpT1Edge" = -log(minp, 10),
          "numEdgeCon" = countContracted
        )
      ))

    else{
      RF = as.numeric(ape::dist.topo(refTree, treem))
      trees = append(trees, list(
        list(
          "tree" = ape::write.tree(treem),
          "logpTree" = -log(pTree, 10),
          "logpT1Edge" = -log(minp, 10),
          "numEdgeCon" = countContracted,
          "RF" = RF
        )
      ))
      RFlist = c(RFlist, RF)# and save RF for plotting
    }

    if (plot > 0) {
      line1 = ""
      if (!is.null(refTree))
        line1 = paste0("RF to reference Tree = ", RF, "; \n")
      plot(
        tree,
        type = "unrooted",
        cex.sub = .7,
        main = paste0("T", countContracted, ": ", testCorrection, ", quad"),
        sub = paste0(
          line1,
          "-log_10(p) for combined T1 over tree = ",
          formatC(-log(pTree, 10), format = 'e', digits = 3),
          "; \n Edge labels are -log_10(p) from combined edge T1;\n",
          length(contract),
          " more edge(s) contracted by combined T1 test."
        )
      )

      if (plot > 1)
      {
        if (plot %in% c(2, 3)) {
          badedgeindex = which(edgep[, "pcomb"] < colorCutoff)
          if (length(badedgeindex > 0))
          {
            badedges = edgep[badedgeindex, "edge"]
            bgColor = rep("tomato1", length(badedges))
            ape::edgelabels(
              formatC(
                -log(edgep[badedgeindex, "pcomb"], 10),
                format = 'e',
                digits = 0
              ),
              edge = badedges,
              bg = bgColor,
              cex = .5
            )
          }
        }
        if (plot == 3)
        {
          contractedEdges = which(tree$edge.length == 0)
          polyList = tree$edge[contractedEdges, 1]
          #colored dots for these nodes
          lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
          internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
          XX <- lastPP$xx[polyList]
          YY <- lastPP$yy[polyList]
          points(XX,
                 YY,
                 pch = 16,
                 col = "red",
                 cex = 1)
        }
        if (plot == 4) {
          edgeind = which(edgep[, "pcomb"] != Inf)
          if (length(edgeind) > 0) {
            bgColor = rep("lightblue", dim(edgep)[[1]])
            bgColor[which(edgep[edgeind, "pcomb"] < colorCutoff)] = "tomato1"
            ape::edgelabels(
              formatC(
                -log(edgep[edgeind, "pcomb"], 10),
                format = 'e',
                digits = 0
              ),
              edge = edgep[edgeind, "edge"],
              bg = bgColor,
              cex = .5
            )
          }
        }
      }
    }

  }
  out$treeList = trees

  if (!is.na(alpha)) {
    small = which(pTreeList < alpha)
    if (length(small) == 0)
      small = 0
    out$indexLate = max(small) + 1
    out$indexEarly = which(pTreeList > alpha)[1]
  } else {
    out$indexLate = NA
    out$indexEarly = NA
  }

  #final plots, of p-values for tree, and RF
  if (plot > 0) {
    y = -log(pTreeList, 10)
    yfinpos = which(is.finite(y))
    yinfpos = setdiff(1:length(y), yfinpos)
    yinf = max(y[yfinpos]) * 1.5
    plot(
      indexTreeList[yinfpos],
      rep(yinf, length(yinfpos)),
      ylim = c(0, yinf),
      xlim = c(-.5, max(indexTreeList) + .5),
      main = paste0("Tree p-values; ", testCorrection, ", quad"),
      xlab = "Number of contracted edges",
      ylab = "-log_10(p) (solid red = Inf)",
      col = "red",
      pch = 19
    )
    points(indexTreeList[yfinpos], y[yfinpos], col = "blue")

    if (!is.null(refTree))
      plot(
        indexTreeList,
        RFlist,
        main = paste0("RFdist to supplied ToB; ", testCorrection, ", quad"),
        xlab = "Tree",
        ylab = "RFdist",
        col = "blue"
      )
  }

  out$alpha=alpha
  out$beta=beta
  
  out = out[c("treeList",
              "indexEarly",
              "indexLate",
              "alpha",
              "beta",
              "qType",
              "testCorrection",
              "pTable")]
  
  invisible(out)
}


####################################################################
#' Edge Contraction to infer a Tree of Blobs - multipartition mode
#'
#' Infer the tree of blobs using edge contraction and multipartitions
#'
#' Given quartet gene tree count data from a network and a starting tree that is
#' a resolution of its tree of blobs, contract edges for which there is poor support
#' due to, first, star-like quartet signals and then successively the the most extreme
#' non-tree-like signals.
#'
#' p-values for star and T1 tests for individual quartets displayed on the tree are
#' combined by one of 4 test correction methods suitable for non-independent tests. For both
#' star test, only p-values for quartets defining each edge are combined. For
#' T1 test, quartets from the current multipartition defining an edge are used.
#'
#' The input tree should be computed by ASTRAL, or another heuristic known to
#' produce a resolution of the tree of blobs in ideal circumstances.
#'
#' Plot levels allow for several graphical depictions of the edge contract process,
#' from the starting tree to a final star tree.
#'
#' For testing with simulated data, a reference tree (e.g., the true Tree of Blobs) 
#' may be provided and the RF distance between it and each tree will be shown.
#'
#' @references
#' \insertRef{ECTo2026}{MSCquartets}
#'
#' \insertRef{MAR19}{MSCquartets}
#'
#' \insertRef{LiuXie2020}{MSCquartets}
#'
#' \insertRef{Chen2022}{MSCquartets}
#' 
#' @param genedata one of:
#' \enumerate{
#' \item a character string giving the name of a file containing a list of newick 
#' gene trees (see argument \code{omit}),
#' \item a multiphylo object containing a collection of gene trees (see argument \code{omit}),
#' \item a table of resolved quartet counts  from a collection of
#' gene trees, or that table with additional columns of quartet p-values. If p_star
#' or both p_T1-D and qindex-D are present (such as from a previous call to this function) 
#' they are not recomputed.}
#' @param tree a tree as a phylo object, or a character string giving the name of 
#' a file with a newick tree, assumed to be a resolution of the true
#' tree of blobs
#' @param testCorrection a method for combining possibly dependent p-values to 
#' get a single p-value, one of "Cauchy", "Bon", "BBC", "CBC"
#' @param alpha test level for choosing inferred tree of blobs
#' @param beta test level for collapsing edges showing no resolution
#' @param colorCutoff level for changing color of edge p-values for T1 test
#' @param plot 0 for no plot;
#'             1 to plot all trees as edges are contracted, along with plot of tree p-values;
#'             2 to add -log_10 of p-values < \code{colorCutoff} for T1 test, or > beta for star test, to edges in plots;
#'             3 to add red dots at multifurcations;
#'             4 to suppress red dots but show all edge p-values
#' @param refTree a reference tree, as a phylo object, to compute the RF distance to 
#' all Tk (for simulation studies); will be unrooted
#' @param omit \code{FALSE} to treat unresolved quartets on gene trees as 1/3 of each resolution;
#' \code{TRUE} to discard unresolved quartet data; ignored if genedata given as quartet table
#' @param lambda parameter for power-divergence statistic (e.g., 0 for likelihood
#' ratio statistic, 1 for Chi-squared statistic)
#' @param method "MLest", "conservative", or "bootstrap"
#' @param smallsample	"precomputed" or "bootstrap", method of obtaining p-value
#' when sample is small (<30)
#' @param smallcounts	"precomputed" or "bootstrap", method of obtaining p-value
#' when some (but not all) counts are small
#' @param bootstraps	number of samples for bootstrapping
#'
#' @returns An invisible list with 8 elements:
#' \describe{
#' \item{$treeList}{A list corresponding to edge contraction steps.}
#' \item{$indexEarly}{The index in \code{$treeList} of the first tree such that
#' the combined tree T1 is above \code{alpha}.}
#' \item{$indexLate}{The index in \code{$treeList} of the first tree such that
#' the combined tree T1 of it and all subsequent ones is above \code{alpha}.}
#' \item{$alpha}{The value of \code{alpha} determining the previous 2 entries.}
#' \item{$beta}{The value of \code{beta} for contracting edges using the 
#' combined star test.}
#' \item{$qType}{The quartet type mode used.}
#' \item{$testCorrection}{The p-value combination method used.}
#' \item{$pTable}{A table of quartet p-values for relevant
#' quartet tests.}
#' }
#' 
#' In \code{$treeList}, the
#' first is the input tree,
#' the second the result of contracting edges judged unresolved
#' by the combined star test,
#' and subsequent ones the result of collapsing edges with the current
#' smallest combined T1 test p-value, with the final one a star tree.
#' Each element in this list itself a
#' list of \describe{
#' \item{$tree}{A tree, with all edge lengths 1 and a node labels
#' of -log of the parent edge's combined p (star for the input tree, T1 for others).}
#' \item{$logpTree}{-log base 10 of the combined T1 p-value of the tree using \code{testCorrection}.}
#' \item{$logpT1Edge}{-log base 10 of the combined T1 p-value of the edge(s)
#'  that were contracted from previous tree (\code{NA}
#' for input tree and second tree).}
#' \item{$numEdgeCon}{The number of edges that have been contracted from the input tree}
#' \item{$RF}{The RF distance between this tree and \code{refTree}, if \code{refTree} is not NULL.}
#' }
#'
#' @examples
#' data(ASTRALtreeLeopardusLescroart)
#' Atree=read.tree(text=ASTRALtreeLeopardusLescroart)
#' data(tableLeopardusLescroart)
#' ECToBlob_mul(tableLeopardusLescroart,Atree,"Bon")
#'
#' @export
#'
ECToBlob_mul = function(genedata,
                        tree,
                        testCorrection = "Bon",
                        alpha = NA,
                        beta = .8,
                        colorCutoff = .05,
                        plot = 3,
                        refTree = NULL,
                        omit=FALSE,
                        lambda = 0,
                        method = "MLest",
                        smallsample = "precomputed",
                        smallcounts = "precomputed",
                        bootstraps = 10 ^ 4
                        )
{
  if (!(testCorrection %in% c("Cauchy", "Bon", "BBC", "CBC")))
    stop("Invalid argument: testCorrection")
  if (!(plot %in% c(0,1,2,3,4)))
    stop("Invalid argument: plot")
  
  #get data on quartets
  if (!(any(c("matrix","character","multiPhylo") %in% class(genedata)))) 
    stop("Data must be supplied as an object of type multiPhylo, character, or matrix.")
  if (!("matrix" %in% class(genedata))) {
    if ("character" %in% class(genedata)) {
      genetrees <- read.tree(file=genedata) #read gene trees
      message(paste("Read", length(genetrees), "gene trees from file."))
    } else { genetrees = genedata} 
    taxanames = genetrees[[1]]$tip.label   # ... get taxa names from first tree
    taxanames = sort(taxanames)
    if (length(taxanames) <= 25) {
      namelist = paste0(taxanames, collapse = ", ")
    } else {
      namelist = paste0(paste0(taxanames[1:25], collapse = ", "),
                        ",...(see output table for full list)")
    }
    message("Analyzing ", length(taxanames), " taxa: ", namelist)
    genedata = quartetTable(genetrees)   # tally quartets on gene trees
    genedata = quartetTableResolved(genedata, omit)   # treat unresolved quartets
  }

  # prepare tree
  if ("character" %in% class(tree)) { tree <- read.tree(file=tree) } #read input tree 
  tree = unroot(tree)
  tree = di2multi(tree) #makes sure tree is unrooted and no zero length branches
  ntaxa = length(tree$tip.label) #number of taxa
  tree$edge.length = rep(1, dim(tree$edge)[1]) #make tree metric, with all edges length 1
  internalEdges = which(tree$edge[, 2] > length(tree$tip.label)) #numbers of internal edges
  if (length(internalEdges) == 0)
    stop("Tree has no internal edges.")
  treeTaxa=sort(tree$tip.label)
  dataTaxa=sort(colnames(genedata)[1:ntaxa])
  if (!all.equal(treeTaxa,dataTaxa))
    stop("Taxa on tree do not match those in genedata")
  numT0edges = dim(tree$edge)[[1]]
  

  out = c()# create empty output
  out$qType = "mul"
  out$testCorrection = testCorrection

  # apply T1-D and star test for all sets of 4 taxa
  pTable = quartetECtestInd(
    genedata,
    tree,
    lambda = lambda,
    method = method,
    smallsample = smallsample,
    smallcounts = smallcounts,
    bootstraps = bootstraps
  )
  out$pTable = pTable

  D = cophenetic(tree) #wipe out p_T1-D for quartets that are resolved
  D = D[order(rownames(D)), order(colnames(D))]
  for (m in 1:dim(pTable)[1]) {
    qnames = which(pTable[m, 1:ntaxa] == 1)
    a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
    b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
    c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
    if (min(c(a, b, c)) == max(c(a, b, c)))
      pTable[m, "p_T1-D"] = NA
  }
  #prepare for indexing table entries
  {
    taxanames = colnames(pTable)[1:ntaxa] #names of taxa, as ordered in table
    np1 = ntaxa + 1 #number of taxa + 1
    #store binomial coefficients for efficiency in computing indices in table to access sets of 4 taxa
    C = matrix(nrow = ntaxa, ncol = 4) #store binomial coefficients for efficiency in computing indices in table to access sets of 4 taxa
    for (i in 0:(ntaxa - 1)) {
      for (j in 1:4) {
        C[i + 1, j] = choose(i, j)
      }
    }
    Cn4 = choose(ntaxa, 4) # need one extra value
  }

  #Compute combined p-values for star and T1 test for all edges
  {
    edgep = matrix(0, length(internalEdges), 2) #create array for combined p-values
    colnames(edgep) = c("edge", "pcomb")

    for (j in 1:length(internalEdges)) {
      #Combine p_star for all quartets defining this edge on original tree.
      pNode = tree$edge[internalEdges[j], 1] #find taxon groups around parent and child node of edge
      cNode = tree$edge[internalEdges[j], 2]
      pGroups = nodeGroups(tree, pNode)
      cGroups = nodeGroups(tree, cNode)

      #remove group above child (depends on nodeGroups returning non-descendants last in list)
      cGroups = cGroups[-length(cGroups)]
      cUnion = unlist(cGroups)

      for (i in 1:length(pGroups)) {
        #remove group below parent
        if (all(cUnion %in% pGroups[[i]]))
          remove = i
      }
      pGroups = pGroups[-remove]

      # compute combined p values for star
      pvecstar = c() #storage for p_star values

      #run through all choices of 2 groups at child and 2 at parent
      for (group1 in 2:length(pGroups)) {
        for (group2 in 1:(group1 - 1)) {
          for (group3 in 2:length(cGroups)) {
            for (group4 in 1:(group3 - 1)) {
              #run through all choices of 1 taxon per group
              for (at in pGroups[[group1]]) {
                for (bt in pGroups[[group2]]) {
                  for (ct in cGroups[[group3]]) {
                    for (dt in cGroups[[group4]]) {
                      #find a,b,c,d entry in table and get p_=star
                      positions = sort(match(tree$tip.label[c(at, bt, ct, dt)], taxanames)) # find table columns for these taxa
                      index = Cn4 - C[np1 - positions[1], 4] - C[np1 - positions[2], 3] - C[np1 - positions[3], 2] - C[np1 - positions[4], 1] #row index for this quartet
                      pstar = unname(pTable[index, "p_star"])
                      pvecstar = c(pvecstar, pstar)
                    }
                  }
                }
              }
            }
          }
        }
      }
      edgep[j, ] = c(internalEdges[j], combineP(pvecstar, testCorrection))
    }
    edgep = rbind(edgep) #make sure result is a matrix if only 1 row
  }

  #Put edge star p-value -logs in as node labels on T_0
  tree$node.label = rep("", tree$Nnode)
  tree$node.label[tree$edge[edgep[, "edge"], 2] - ntaxa] = -log(edgep[, "pcomb"], 10)

  #Put initial tree T_0 in list, with special treatment if reference Tree supplied
  pTree = combineP(na.omit(pTable[, "p_T1-D"]), testCorrection) # compute combined values for tree T_0
  pTreeList = c(pTree) # save in list
  indexTreeList = c(0)
  if (is.null(refTree))
    trees = list(list(
      "tree" = ape::write.tree(tree),
      "logpTree" = -log(pTree, 10),
      "logpT1Edge" = NA,
      "numEdgeCon" = 0
    ))
  else{
    refTree=unroot(refTree)
    RF = as.numeric(ape::dist.topo(refTree, tree))
    trees = list(
      list(
        "tree" = ape::write.tree(tree),
        "logpTree" = -log(pTree, 10),
        "logpT1Edge" = NA,
        "numEdgeCon" = 0,
        "RF" = RF
      )
    )# put it in list
    RFlist = c(RF)# and save RF for plotting
  }

  if (plot > 0) {
    #Plot initial tree
    line1 = ""
    if (!is.null(refTree))
      line1 = paste0("RF to reference Tree = ", RF, "; \n")

    plot(
      tree,
      type = "unrooted",
      cex.sub = .7,
      main = paste0("T0 (Input tree): ", testCorrection, ", mul"),
      sub = paste0(
        line1,
        "-log_10(p) for combined T1 over tree = ",
        formatC(-log(pTree, 10), format = 'e', digits = 3),
        "; \n Edge labels are -log_10(p) for combined star tests for that edge."
      )
    )
    if (plot %in% c(2, 3)) {
      badedgeindex = which(edgep[, "pcomb"] > beta)
      if (length(badedgeindex) > 0) {
        badedges = edgep[badedgeindex, "edge"]
        bgColor = rep("tomato1", length(badedges))
        ape::edgelabels(
          formatC(
            -log(edgep[badedgeindex, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = badedges,
          bg = bgColor,
          cex = .5
        )
      }
    }
    if (plot == 3) {
      polyList = c() #find multifurcations, starting at "top" in unrooted tree
      if (length(Descendants(tree, ntaxa + 1, 'children')) > 3)
        polyList = c(polyList, ntaxa + 1)
      if (tree$Nnode >= 2) {
        #loop over other internal nodes
        for (i in (ntaxa + 2):(ntaxa + tree$Nnode)) {
          if (length(Descendants(tree, i, 'children')) > 2)
            polyList = c(polyList, i)
        }
      }
      #add colored dots for these nodes
      lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
      XX <- lastPP$xx[polyList]
      YY <- lastPP$yy[polyList]
      points(XX,
             YY,
             pch = 16,
             col = "red",
             cex = 1)
    }
    if (plot == 4) {
      edgeind = which(edgep[, "pcomb"] != Inf)
      if (length(edgeind) > 0) {
        bgColor = rep("lightblue", dim(edgep)[[1]])
        bgColor[which(edgep[edgeind, "pcomb"] > beta)] = "tomato1"
        ape::edgelabels(
          formatC(
            -log(edgep[edgeind, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = edgep[edgeind, "edge"],
          bg = bgColor,
          cex = .5
        )
      }
    }
  }

  #Contract edges based on star test
  message('Contracting edges with star test at level beta = ',beta)
  edgepRow = which(edgep[, "pcomb"] > beta) #determine entries for edges to contract for star test
  edges2contract = edgep[edgepRow, "edge"] # and the edge numbers
  num2contract = length(edges2contract)
  tree$edge.length[edges2contract] = 0 #contract those edges
  tree = di2multi(tree) #replace tree


  D = cophenetic(tree) #wipe out p_T1-D for quartets no longer resolved
  D = D[order(rownames(D)), order(colnames(D))]
  for (m in 1:dim(pTable)[1]) {
    qnames = which(pTable[m, 1:ntaxa] == 1)
    a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
    b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
    c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
    if (min(c(a, b, c)) == max(c(a, b, c)))
      pTable[m, "p_T1-D"] = NA
  }

  # Computation of initial combined edge p for T1 test, using tree after star test contraction
  internalEdges = which(tree$edge[, 2] > ntaxa) #numbers of internal edges
  edgep = matrix(0, length(internalEdges), 2) #create array for combined p-values
  colnames(edgep) = c("edge", "pcomb")

  for (j in 1:length(internalEdges)) {
    #Combine p_T1-D for all quartets defining this edge on original tree.
    pNode = tree$edge[internalEdges[j], 1] #find taxon groups around parent and child node of edge
    cNode = tree$edge[internalEdges[j], 2]
    pGroups = nodeGroups(tree, pNode)
    cGroups = nodeGroups(tree, cNode)

    #remove group above child (depends on nodeGroups returning non-descendants last in list)
    cGroups = cGroups[-length(cGroups)]
    cUnion = unlist(cGroups)

    for (i in 1:length(pGroups)) {
      #remove group below parent
      if (all(cUnion %in% pGroups[[i]]))
        remove = i
    }
    pGroups = pGroups[-remove]

    # compute combined p values for T1
    pvecT1 = c() #storage for p_T1-D values

    #run through all choices of 2 groups at child and 2 at parent
    for (group1 in 2:length(pGroups)) {
      for (group2 in 1:(group1 - 1)) {
        for (group3 in 2:length(cGroups)) {
          for (group4 in 1:(group3 - 1)) {
            #run through all choices of 1 taxon per group
            for (at in pGroups[[group1]]) {
              for (bt in pGroups[[group2]]) {
                for (ct in cGroups[[group3]]) {
                  for (dt in cGroups[[group4]]) {
                    #find a,b,c,d entry in table and get p_=star
                    positions = sort(match(tree$tip.label[c(at, bt, ct, dt)], taxanames)) # find table columns for these taxa
                    index = Cn4 - C[np1 - positions[1], 4] - C[np1 - positions[2], 3] - C[np1 - positions[3], 2] - C[np1 - positions[4], 1] #row index for this quartet
                    pT1 = unname(pTable[index, "p_T1-D"])
                    pvecT1 = c(pvecT1, pT1)
                  }
                }
              }
            }
          }
        }
      }
    }
    edgep[j, ] = c(internalEdges[j], combineP(pvecT1, testCorrection))
  }
  edgep = rbind(edgep) #make sure result is a matrix if only 1 row


  #Put edge T1 p-value -logs in as node labels on Tk
  tree$node.label = rep("", tree$Nnode)
  tree$node.label[tree$edge[edgep[, "edge"], 2] - ntaxa] = -log(edgep[, "pcomb"], 10)


  #Put star-contracted tree Tk in list, with special treatment if reference Tree supplied
  pTree = combineP(na.omit(pTable[, "p_T1-D"]), testCorrection) # compute combined values for tree T_0
  pTreeList = c(pTreeList, pTree) # save in list
  indexTreeList = c(indexTreeList, num2contract)
  if (is.null(refTree))
    trees = append(trees, list(
      list(
        "tree" = ape::write.tree(tree),
        "logpTree" = -log(pTree, 10),
        "logpT1Edge" = NA,
        "numEdgeCon" = num2contract
      )
    ))
  else{
    RF = as.numeric(ape::dist.topo(refTree, tree))
    trees = append(trees, list(
      list(
        "tree" = ape::write.tree(tree),
        "logpTree" = -log(pTree, 10),
        "logpT1Edge" = NA,
        "numEdgeCon" = num2contract,
        "RF" = RF
      )
    ))# put it in list
    RFlist = c(RFlist,RF)# and save RF for plotting
  }

  if (plot > 0) {
    #Plot this tree
    line1 = ""
    if (!is.null(refTree))
      line1 = paste0("RF to reference Tree = ", RF, "; \n")

    plot(
      tree,
      type = "unrooted",
      cex.sub = .7,
      main = paste0("T", num2contract, ": ", testCorrection, ", mul"),
      sub = paste0(
        line1,
        "-log_10(p) for combined T1 over tree = ",
        formatC(-log(pTree, 10), format = 'e', digits = 3),
        "; \n Edge labels are -log_10(p) for combined T1 test for that edge;\n",
        num2contract,
        " edges contracted from by combined star test."
      )
    )
    if (plot %in% c(2, 3)) {
      badedgeindex = which(edgep[, "pcomb"] < colorCutoff)
      if (length(badedgeindex) > 0) {
        badedges = edgep[badedgeindex, "edge"]
        bgColor = rep("tomato1", length(badedges))
        ape::edgelabels(
          formatC(
            -log(edgep[badedgeindex, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = badedges,
          bg = bgColor,
          cex = .5
        )
      }
    }
    if (plot == 3) {
      polyList = c() #find multifurcations, starting at "top" in unrooted tree
      if (length(Descendants(tree, ntaxa + 1, 'children')) > 3)
        polyList = c(polyList, ntaxa + 1)
      if (tree$Nnode >= 2) {
        #loop over other internal nodes
        for (i in (ntaxa + 2):(ntaxa + tree$Nnode)) {
          if (length(Descendants(tree, i, 'children')) > 2)
            polyList = c(polyList, i)
        }
      }
      #add colored dots for these nodes
      lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
      XX <- lastPP$xx[polyList]
      YY <- lastPP$yy[polyList]
      points(XX,
             YY,
             pch = 16,
             col = "red",
             cex = 1)
    }
    if (plot == 4) {
      edgeind = which(edgep[, "pcomb"] != Inf)
      if (length(edgeind) > 0) {
        bgColor = rep("lightblue", dim(edgep)[[1]])
        bgColor[which(edgep[edgeind, "pcomb"] < colorCutoff)] = "tomato1"
        ape::edgelabels(
          formatC(
            -log(edgep[edgeind, "pcomb"], 10),
            format = 'e',
            digits = 0
          ),
          edge = edgep[edgeind, "edge"],
          bg = bgColor,
          cex = .5
        )
      }
    }
  }

  numEdges = length(internalEdges)

  while (numEdges > 0)
    #loop with successive contracting steps
  {
    minp = min(edgep[, "pcomb"]) # p for edges to contract
    contract = which(edgep[, "pcomb"] == minp) # find all minimal edge p
    tree$edge.length[edgep[contract, 1]] = 0 # and set their lengths to 0
    edgep[contract, 2] = Inf # and never consider them again
    tree = di2multi(tree)

    # Computation of combined edge p for T1 test, using tree after last contraction
    internalEdges = which(tree$edge[, 2] > ntaxa) #numbers of internal edges
    if (length(internalEdges) > 0) {
      edgep = matrix(0, length(internalEdges), 2) #create array for combined p-values
      colnames(edgep) = c("edge", "pcomb")

      for (j in 1:length(internalEdges)) {
        pNode = tree$edge[internalEdges[j], 1] #find taxon groups around parent and child node of edge
        cNode = tree$edge[internalEdges[j], 2]
        pGroups = nodeGroups(tree, pNode)
        cGroups = nodeGroups(tree, cNode)

        #remove group above child (depends on nodeGroups returning non-descendants last in list)
        cGroups = cGroups[-length(cGroups)]
        cUnion = unlist(cGroups)

        for (i in 1:length(pGroups)) {
          #remove group below parent
          if (all(cUnion %in% pGroups[[i]]))
            remove = i
        }
        pGroups = pGroups[-remove]

        # compute combined p values for T1
        pvecT1 = c() #storage for p_T1-D values

        #run through all choices of 2 groups at child and 2 at parent
        for (group1 in 2:length(pGroups)) {
          for (group2 in 1:(group1 - 1)) {
            for (group3 in 2:length(cGroups)) {
              for (group4 in 1:(group3 - 1)) {
                #run through all choices of 1 taxon per group
                for (at in pGroups[[group1]]) {
                  for (bt in pGroups[[group2]]) {
                    for (ct in cGroups[[group3]]) {
                      for (dt in cGroups[[group4]]) {
                        #find a,b,c,d entry in table and get p_=star
                        positions = sort(match(tree$tip.label[c(at, bt, ct, dt)], taxanames)) # find table columns for these taxa
                        index = Cn4 - C[np1 - positions[1], 4] - C[np1 - positions[2], 3] - C[np1 - positions[3], 2] - C[np1 - positions[4], 1] #row index for this quartet
                        pT1 = unname(pTable[index, "p_T1-D"])
                        pvecT1 = c(pvecT1, pT1)
                      }
                    }
                  }
                }
              }
            }
          }
        }
        edgep[j, ] = c(internalEdges[j], combineP(pvecT1, testCorrection))
      }
      edgep = rbind(edgep) #make sure result is a matrix if only 1 row

      #Put edge T1 p-value -logs in as node labels on Tk
      tree$node.label = rep("", tree$Nnode)
      tree$node.label[tree$edge[edgep[, "edge"], 2] - ntaxa] = -log(edgep[, "pcomb"], 10)
    }

    D = cophenetic(tree) #wipe out p_T1-D for quartets no longer resolved
    D = D[order(rownames(D)), order(colnames(D))]
    for (m in 1:dim(pTable)[1]) {
      qnames = which(pTable[m, 1:ntaxa] == 1)
      a = D[qnames[1], qnames[2]] + D[qnames[3], qnames[4]]
      b = D[qnames[1], qnames[3]] + D[qnames[2], qnames[4]]
      c = D[qnames[1], qnames[4]] + D[qnames[2], qnames[3]]
      if (min(c(a, b, c)) == max(c(a, b, c)))
        pTable[m, "p_T1-D"] = NA
    }

    countContracted = numT0edges - dim(tree$edge)[[1]] #total number of edges contracted
    pTree = combineP(na.omit(pTable[, "p_T1-D"]), testCorrection) # compute new combined values for tree
    pTreeList = c(pTreeList, pTree) # save new tree p
    indexTreeList = c(indexTreeList, countContracted)

    if (is.null(refTree))
      trees = append(trees, list(
        list(
          "tree" = ape::write.tree(tree),
          "logpTree" = -log(pTree, 10),
          "logpT1Edge" = -log(minp, 10),
          "numEdgeCon" = countContracted
        )
      ))

    else{
      RF = as.numeric(ape::dist.topo(refTree, tree))
      trees = append(trees, list(
        list(
          "tree" = ape::write.tree(tree),
          "logpTree" = -log(pTree, 10),
          "logpT1Edge" = -log(minp, 10),
          "numEdgeCon" = countContracted,
          "RF" = RF
        )
      ))
      RFlist = c(RFlist, RF)# and save RF for plotting
    }

    numEdges = dim(tree$edge)[[1]] - ntaxa  #how many edges left to collapse?

    if (plot > 0) {
      line1 = ""
      if (!is.null(refTree))
        line1 = paste0("RF to reference Tree = ", RF, "; \n")
      plot(
        tree,
        type = "unrooted",
        cex.sub = .7,
        main = paste0("T", countContracted, ": ", testCorrection, ", mul"),
        sub = paste0(
          line1,
          "-log_10(p) for combined T1 over tree = ",
          formatC(-log(pTree, 10), format = 'e', digits = 3),
          "; \n Edge labels are -log_10(p) from combined edge T1;\n",
          length(contract),
          " more edge(s) contracted by combined T1 test."
        )
      )

      if (plot > 1)
      {
        if (plot %in% c(2, 3)) {
          badedgeindex = which(edgep[, "pcomb"] < colorCutoff)
          if (length(badedgeindex > 0))
          {
            badedges = edgep[badedgeindex, "edge"]
            bgColor = rep("tomato1", length(badedges))
            ape::edgelabels(
              formatC(
                -log(edgep[badedgeindex, "pcomb"], 10),
                format = 'e',
                digits = 0
              ),
              edge = badedges,
              bg = bgColor,
              cex = .5
            )
          }
        }
        if (plot == 3)
        {
          polyList = c() # find multifurcations
          if (length(Descendants(tree, ntaxa + 1, 'children')) > 3) {
            polyList = c(polyList, ntaxa + 1)
          }
          if (tree$Nnode >= 2) {
            for (i in (ntaxa + 2):(ntaxa + tree$Nnode)) {
              #loop over internal nodes
              if (length(Descendants(tree, i, 'children')) > 2) {
                polyList = c(polyList, i)
              }
            }
          }
          
          #colored dots for these nodes
          lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
          internal_nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
          XX <- lastPP$xx[polyList]
          YY <- lastPP$yy[polyList]
          points(XX,
                 YY,
                 pch = 16,
                 col = "red",
                 cex = 1)
        }
        if (plot == 4) {
          edgeind = which(edgep[, "pcomb"] != Inf)
          if (length(edgeind) > 0) {
            bgColor = rep("lightblue", dim(edgep)[[1]])
            bgColor[which(edgep[edgeind, "pcomb"] < colorCutoff)] = "tomato1"
            ape::edgelabels(
              formatC(
                -log(edgep[edgeind, "pcomb"], 10),
                format = 'e',
                digits = 0
              ),
              edge = edgep[edgeind, "edge"],
              bg = bgColor,
              cex = .5
            )
          }
        }
      }
    }

  }
  out$treeList = trees

  if (!is.na(alpha)) {
    small = which(pTreeList < alpha)
    if (length(small) == 0)
      small = 0
    out$indexLate = max(small) + 1
    out$indexEarly = which(pTreeList > alpha)[1]
  } else {
    out$indexLate = NA
    out$indexEarly = NA
  }

  #final plot, of p-values for tree, and RF
  if (plot > 0) {
    y = -log(pTreeList, 10)
    yfinpos = which(is.finite(y))
    yinfpos = setdiff(1:length(y), yfinpos)
    yinf = max(y[yfinpos]) * 1.5
    plot(
      indexTreeList[yinfpos],
      rep(yinf, length(yinfpos)),
      ylim = c(0, yinf),
      xlim = c(-.5, max(indexTreeList) + .5),
      main = paste0("Tree p-values; ", testCorrection, ", mul"),
      xlab = "Number of contracted edges",
      ylab = "-log_10(p) (solid red = Inf)",
      col = "red",
      pch = 19
    )
    points(indexTreeList[yfinpos], y[yfinpos], col = "blue")

    if (!is.null(refTree))
      plot(
        indexTreeList,
        RFlist,
        main = paste0("RFdist to supplied ToB; ", testCorrection, ", mul"),
        xlab = "Tree",
        ylab = "RFdist",
        col = "blue"
      )
   
  }

  out$alpha=alpha
  out$beta=beta
  
  out = out[c("treeList",
              "indexEarly",
              "indexLate",
              "alpha",
              "beta",
              "qType",
              "testCorrection",
              "pTable")]
  
  invisible(out)
}
