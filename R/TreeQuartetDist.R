#' Compute quartet distance between taxa
#'
#' Compute the Quartet Distance of \insertCite{Rho19;textual}{MSCquartets} from a table specifying a collection of quartets on
#' \code{n} taxa.
#'
#' @references 
#' \insertRef{Rho19}{MSCquartets}
#' 
#' @param dqt an (\code{n} choose 4) x \code{n} (or \code{n+1}) matrix of form output by \code{quartetTableDominant};
#'  (Note: If present, the \code{n+1}th column of \code{dqt} is ignored.)
#'
#' @return a pairwise distance matrix on \code{n} taxa
#'
#' @seealso 
#' \code{\link{quartetTableDominant}},
#' \code{\link{QDS}},
#' \code{\link{QDC}},
#' \code{\link{quartetWeightedDist}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' DQT=quartetTableDominant(RQT)
#' Dist=quartetDist(DQT)
#' tree=NJ(Dist)
#' write.tree(tree)
#' plot(tree)
#'
#' @export
quartetDist = function(dqt) {
  Qcolnames = colnames(dqt)# get column names
  n = length(Qcolnames) - 1# number of taxa
  taxa = Qcolnames[1:n]# names of taxa
  m = dim(dqt)[1]# number of quartets
  
  D = matrix(0, n, n)# create distance matrix
  colnames(D) = taxa
  rownames(D) = taxa
  
  for (i in 1:m) {
    ones = which(dqt[i, ] == 1)
    if (length(ones) == 2) {
      # skip null quartets
      nones = which(dqt[i, ] == -1)
      D[ones, nones] = D[ones, nones] + 1 # increment 4 entries of D
      D[nones, ones] = D[nones, ones] + 1 # and another 4
    }
  }
  
  D = 2 * D + 2 * n - 4 # adjust to appropriate Q distance
  for (i in 1:n) {
    D[i, i] = 0 # and unadjust diagonal back to 0
  }
  return(D)
}

###################################################################

#' Compute Quartet Distance Supertree
#'
#' Apply the Quartet Distance Supertree method of \insertCite{Rho19;textual}{MSCquartets} to a table specifying a
#' collection of quartets on \code{n} taxa.
#'
#' @references 
#' \insertRef{Rho19}{MSCquartets}
#' 
#' @details This function is a wrapper which runs \code{quartetDist} and then builds a tree.
#'
#' @param dqt an (\code{n} choose 4) x \code{n} (or \code{n+1}) matrix of form output by \code{quartetTableDominant};
#'  (Note: If present, the \code{n+1}th column of \code{dqt} is ignored)
#' 
#' @param method tree building function (e.g., fastme.bal, nj)
#' 
#' @return
#'     an unrooted metric tree of type phylo. Edge lengths are not in interpretable units
#' @seealso 
#' \code{\link{quartetTableDominant}},
#' \code{\link{quartetDist}},
#' \code{\link{QDC}},
#' \code{\link{WQDS}}, 
#' \code{\link{WQDC}}, 
#' \code{\link{WQDCrecursive}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' DQT=quartetTableDominant(RQT)
#' tree=QDS(DQT)
#' write.tree(tree)
#' plot(tree)
#' 
#' @importFrom ape fastme.bal
#'
#' @export
QDS = function(dqt,
               method = fastme.bal) {
  D = quartetDist(dqt)
  QDS = unroot(method(D))# build tree
  return(QDS)
}

#################################################

#' Compute the Weighted Quartet Distance between taxa
#'
#' Compute the Weighted Quartet Distance between taxa of \insertCite{YR19;textual}{MSCquartets} from a table specifying a collection of quartets on
#' \code{n} taxa and the quartets' internal branch lengths.
#' 
#' @references 
#' \insertRef{YR19}{MSCquartets}
#' 
#' @param dqt an (\code{n} choose 4) x (\code{n+1}) matrix of the form output by \code{quartetTableDominant}
#'
#' @return
#'     a pairwise distance matrix on \code{n} taxa
#'
#' @seealso \code{\link{quartetTableDominant}},
#'           \code{\link{WQDSAdjustLengths}},
#'           \code{\link{WQDS}},
#'           \code{\link{WQDC}}, 
#'           \code{\link{WQDCrecursive}},
#'           \code{\link{quartetWeightedDist}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' DQT=quartetTableDominant(RQT,bigweights="finite")
#' D=quartetWeightedDist(DQT)
#' tree=NJ(D)
#' stree=WQDSAdjustLengths(tree)
#' write.tree(stree)
#'
#' @export
quartetWeightedDist = function(dqt) {
  Qcolnames = colnames(dqt)# get column names
  n = length(Qcolnames) - 1# number of taxa
  taxa = Qcolnames[1:n]# names of taxa
  m = dim(dqt)[1]# number of quartets
  w = dqt[, n + 1]# weights
  
  D = matrix(2, n, n)# create distance matrix
  colnames(D) = taxa
  rownames(D) = taxa
  
  for (i in 1:m) {
    ones = which(dqt[i,] == 1)
    if (length(ones) == 2) {
      # skip null quartets
      nones = which(dqt[i,] == -1)
      D[ones, nones] = D[ones, nones] + w[i] # increment 4 entries of D
      D[nones, ones] = D[nones, ones] + w[i] # and another 4
    }
  }
  
  for (i in 1:n) {
    D[i, i] = 0 # unadjust diagonal back to 0
  }
  return(D)
}


#################################################

#' Compute the Weighted Quartet Distance Supertree 
#'
#' Apply the Weighted Quartet Distance Supertree method of \insertCite{YR19;textual}{MSCquartets} to
#' a collection of quartets on \code{n} taxa together with internal 
#' quartet branch lengths, specified by a table.
#' 
#' @references 
#' \insertRef{YR19}{MSCquartets}
#'
#' @details This function is a wrapper which runs \code{quartetWeightedDist}, builds a tree, and then adjusts edge lengths
#' with \code{WQDSAdjustLengths}.
#'
#' @param dqt an (\code{n} choose 4) x \code{n+1}) matrix of form output by \code{quartetTableDominant}
#' @param method a distance-based tree building function (e.g., fastme.bal, NJ, etc.)
#'
#' @return an unrooted metric tree, of type phylo
#'          
#' @seealso \code{\link{quartetTableDominant}}, 
#'          \code{\link{quartetWeightedDist}},
#'          \code{\link{WQDSAdjustLengths}},
#'          \code{\link{WQDC}}, 
#'          \code{\link{WQDCrecursive}},
#'          \code{\link{QDS}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' DQT=quartetTableDominant(RQT,bigweights= "finite")
#' tree=WQDS(DQT)
#' write.tree(tree)
#' plot(tree)
#'
#' @importFrom ape fastme.bal
#' @export
WQDS = function(dqt,
                method = fastme.bal) {
  D = quartetWeightedDist(dqt)# compute distance
  WQDS = unroot(method(D))# build tree
  WQDS=WQDSAdjustLengths(WQDS)
  return(WQDS)
}

#################################################

#' Adjust edge lengths on tree built from Weighted Quartet distance 
#' to estimate metric tree
#'
#' Modify edge lengths of a tree built from a distance table produced by \code{quartetWeightedDist},
#' to remove scaling factors related to the size of the split associated to the edge.
#'
#' @details As explained by \insertCite{YR19;textual}{MSCquartets}, a metric tree produced from
#' the weighted quartet distance has edge lengths
#' inflated by a factor dependent on the associated split size. Removing these
#' factors 
#' yields a consistent estimate of the metric species tree displaying the weighted
#' quartets, if such a tree exists.
#'
#' This function should not be used on trees output from \code{WQDS}, \code{WQDC},
#' or \code{WQDCrecursive}, as 
#' their edges are already adjusted. It can be used on trees built from the distance
#' computed by \code{quartetWeightedDist}.
#' 
#' @references 
#' \insertRef{YR19}{MSCquartets}
#'  
#' @param tree an unrooted metric tree, of type phylo
#' 
#' @return an unrooted metric tree, of type phylo
#' 
#' @seealso 
#' \code{\link{WQDS}},
#' \code{\link{WQDC}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' DQT=quartetTableDominant(RQT,bigweights="finite")
#' D=quartetWeightedDist(DQT)
#' tree=NJ(D)
#' write.tree(tree)
#' plot(tree)
#' stree=WQDSAdjustLengths(tree)
#' write.tree(stree)
#' plot(stree)
#'
#' @importFrom ape is.rooted unroot root
#' @importFrom phangorn Descendants
#' @export
WQDSAdjustLengths = function(tree) {
  if (is.rooted(tree) == TRUE) {
    message("Tree is rooted; unrooting it.")
    tree = unroot(tree)
  }
  ntree = root(tree, outgroup = tree$tip.label[1])
  ntaxa = length(ntree$tip.label)
  nedges = dim(ntree$edge)[1]
  for (i in 1:nedges) {
    node = ntree$edge[i, 2] # get child of edge
    ndesc = length(Descendants(ntree, node, "tips")[[1]])# count descendents, using phangorn function
    if (ndesc == 1) {
      ntree$edge.length[i] = 1
    }
    else if (ndesc == ntaxa - 1) {
      ntree$edge.length[i] = 0
    }
    else {
      ntree$edge.length[i] = ntree$edge.length[i] / ((ndesc - 1) * (ntaxa - ndesc - 1)) # adjust edge length
    }
  }
  
  WQDSAdjustLengths = unroot(ntree)
  return(WQDSAdjustLengths)
}


#################################################


#' Compute Quartet Distance Consensus tree from gene tree data
#'
#' Compute the Quartet Distance Consensus \insertCite{Rho19}{MSCquartets} estimate of an unrooted 
#' topological species tree from gene tree data.
#'
#' @details This function is a wrapper which performs the steps of reading in a collection
#' of gene trees, tallying quartets, computing the quartet distance between taxa, building
#' a tree which consistently estimates the unrooted species tree topology under the MSC, and then possibly estimating edge
#' lengths using the "simpleML" method.
#'
#' @references 
#' \insertRef{Rho19}{MSCquartets}
#' 
#' @param genetreedata  gene tree data that may be supplied in any of 3 forms: 
#' \enumerate{
#' \item a character string giving the name of a file containing gene trees in Newick,
#' \item a multiPhylo object containing gene trees, or
#' \item a resolved quartet table, such as produced by \code{quartetTableResolved}
#' } 
#' @param taxanames if \code{genetreedata} is a file or a multiPhylo object, a vector of a subset  
#' of the taxa names on the gene trees 
#' to be analyzed, if \code{NULL} all taxa on the first gene tree are used; if \code{genetreedata} 
#' is a quartet table, this argument is ignored and all taxa in the table are used
#' @param method a distance-based tree building function, such as \code{fastme.bal} or \code{nj}
#' @param omit \code{TRUE} ignores unresolved quartets; \code{FALSE} treats them as 1/3 of each resolution; 
#' ignored if \code{genetreedata} is supplied as a quartet table
#' @param metric if \code{FALSE} return topological tree; if \code{TRUE} return metric tree with
#' internal edge lengths estimated by \code{estimateEdgeLengths} with \code{lambda=0}, and terminal branches of length 1
#' 
#' 
#' @return an unrooted tree of type phylo
#' 
#' @seealso \code{\link{quartetTable}},
#'          \code{\link{quartetTableResolved}},
#'          \code{\link{quartetTableDominant}},
#'          \code{\link{quartetDist}},
#'          \code{\link{QDS}},
#'          \code{\link{WQDC}}, 
#'          \code{\link{WQDCrecursive}}
#'          \code{\link{estimateEdgeLengths}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' stree=QDC(gtrees,tnames[1:6])
#' write.tree(stree)
#' plot(stree)
#' streeMetric=QDC(gtrees, tnames[1:6],metric=TRUE)
#' write.tree(streeMetric)
#' plot(streeMetric)
#' 
#' @export
QDC = function(genetreedata,
               taxanames = NULL,
               method=fastme.bal,
               omit = FALSE, 
               metric=FALSE) {
  if ("matrix" %in% class(genetreedata)) {
    RQT = genetreedata
    if (!is.null(taxanames)) {
      message("Ignoring argument taxanames since genetreedata supplied as quartet table.")
    }
  } else {
    if ("multiPhylo" %in% class(genetreedata))  {
      genetrees = genetreedata
    } else {
      if ("character" %in% class(genetreedata)) {
        genetrees <- read.tree(genetreedata) #read gene trees
        message(paste("Read", length(genetrees), "gene trees from file."))
      } else {
        stop("Data must be supplied as an object of class multiPhylo, character, or matrix.")
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
    
    # tally quartets on gene trees
    RQT = quartetTableResolved(quartetTable(genetrees, taxanames), omit = omit)         # treat unresolved quartets as 1/3 of each resolution
  }
  
  Q = quartetTableDominant(RQT)
  tree = QDS(Q,method=method)
  if (metric==FALSE) {
    tree$edge.length = NULL
  } else {
    tree = estimateEdgeLengths(tree, RQT, lambda = 0)
  }
  return(tree)
}


#############################################################

#' Compute Weighted Quartet Distance Consensus tree from gene tree data
#'
#' Compute the Weighted Quartet Distance Consensus \insertCite{YR19}{MSCquartets} estimate of a 
#' species tree from gene tree data. This is a consistent estimator of the unrooted 
#' species tree topology and all internal branch lengths.
#'
#' @details This function is a wrapper which performs the steps of reading in a collection
#' of gene trees, tallying quartets, estimating quartet internal branch lengths, computing the weighted
#' quartet distance between taxa, building
#' a tree, and adjusting edge lengths, to give a consistent estimate of the metric species tree in coalescent units
#' under the MSC.
#' 
#' If the gene tree data indicates some quartets experienced little to no incomplete lineage 
#' sorting, this algorithm tends to be less topologically accurate than \code{QDC} 
#' (which infers no metric information) or \code{WQDCrecursive} (which gives better topologies,
#' and reasonably accurate lengths for short edges, though long edge lengths may still be unreliable).
#' 
#' @references 
#' \insertRef{YR19}{MSCquartets}
#' 
#' @param genetreedata  gene tree data that may be supplied in any of 3 forms:
#' \enumerate{
#' \item a character string giving the name of a file containing gene trees in Newick
#' \item a multiPhylo object containing gene trees
#' \item a resolved quartet table, as produced by \code{quartetTableResolved}
#' }
#' @param taxanames if \code{genetreedata} is a file or a multiPhylo object, a vector of a subset  
#' of the taxa names on the gene trees 
#' to be analyzed, if \code{NULL} all taxa on the first gene tree are used; if \code{genetreedata} 
#' is a quartet table, this argument is ignored and all taxa in the table are used
#' @param method a distance-based tree building function, such as \code{fastme.bal} or \code{nj}
#' @param omit \code{TRUE} leaves out unresolved quartets, \code{FALSE} treats them as 1/3 of each resolution; ignored if 
#' \code{genetreedata} is given as  a resolved quartet table
#' @param terminal non-negative branch length to supply for terminal branches 
#' whose length cannot be inferred by \code{WQDC}
#' @return an unrooted metric tree of type phylo
#' @seealso \code{\link{quartetTable}}, 
#'          \code{\link{quartetTableResolved}}, 
#'          \code{\link{quartetTableDominant}}, 
#'          \code{\link{quartetWeightedDist}}, 
#'          \code{\link{WQDCrecursive}}, 
#'          \code{\link{WQDS}}, 
#'          \code{\link{QDC}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' stree=WQDC(gtrees,tnames[1:6])
#' write.tree(stree)
#' plot(stree)
#'
#' @export
WQDC = function(genetreedata,
                taxanames = NULL,
                method=fastme.bal,
                omit = FALSE,
                terminal = 1) {
  
  if ("matrix" %in% class(genetreedata)) {
    RQT = genetreedata
    if (!is.null(taxanames)) {
      message("Ignoring argument taxanames since genetreedata supplied as quartet table.")
    }
  } else {
    if ("multiPhylo" %in% class(genetreedata))  {
      genetrees = genetreedata
    } else {
      if ("character" %in% class(genetreedata)) {
        genetrees <- read.tree(genetreedata) #read gene trees
        message(paste("Read", length(genetrees), "gene trees from file."))
      } else {
        stop("Data must be supplied as an object of class multiPhylo, character, or matrix.")
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
    
    # tally quartets on gene trees
    RQT = quartetTableResolved(quartetTable(genetrees, taxanames), omit = omit)         # treat unresolved quartets as 1/3 of each resolution
  }
  
  Q=quartetTableDominant(RQT, bigweights ="finite")
  
  WQDC = WQDS(Q,method=method)
  term_edges = which(WQDC$edge[, 2] <= length(taxanames)) #modify terminal edge lengths
  WQDC$edge.length[term_edges] = terminal
  return(WQDC)
}
