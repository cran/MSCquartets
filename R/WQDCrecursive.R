#' Compute the Recursive Weighted Quartet Distance Consensus tree from gene tree data
#'
#' Infer a metric species tree from counts
#' of quartets displayed on a collection of gene trees, as described by 
#' \insertCite{YR19;textual}{MSCquartets}. Edge lengths are in coalescent units.
#'
#'
#' @details The algorithm counts quartets displayed on the gene trees, builds a tree using \code{WQDS},
#' determines the split corresponding to the longest edge in that tree,
#' and then recursively builds trees
#' on the taxa in each split set together with a `composite taxon' formed by all
#' taxa in the other split set.
#' This approach is slower than non-recursive \code{WQDC}, but increases topological accuracy. Shorter branch 
#' lengths tend to be more accurately estimated.
#'
#' This function must be called with its argument a resolved quartet
#' table of size (n choose 4)x(n+3). Its recursive nature
#' requires building smaller resolved quartet tables on split sets with an additional
#' composite taxon.
#'
#' @references
#' \insertRef{YR19}{MSCquartets}
#'
#' @param rqt a resolved quartet table as produced by \code{quartetTableResolved}
#' @param method a distance-based tree building function, such as \code{fastme.bal} or \code{nj}
#' @param stopAt a non-negative branch length in coalescent units; recursive calls stop when the longest 
#' branch in a recursively examined subtree is smaller than this value
#' @param terminal non-negative branch length to supply for terminal branches,
#' whose lengths cannot be inferred by \code{WQDCrecursive}
#'
#' @return an unrooted metric tree, of type phylo
#'
#' @seealso \code{\link{quartetTableResolved}},\code{\link{quartetTable}},
#' \code{\link{QDC}}, \code{\link{QDS}}, \code{\link{quartetTableCollapse}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' stree=WQDCrecursive(RQT)
#' write.tree(stree)
#' plot(stree)
#'
#' @importFrom ape Ntip drop.tip bind.tree
#'
#' @export
WQDCrecursive = function(rqt,
                         method = fastme.bal,
                         stopAt = 2,
                         terminal = 1) {
  nX = dim(rqt)[2] - 3 # number of taxa
  if (nX < 4) {
    stop("Fewer than 4 taxa.")
  }
  
  taxanames = colnames(rqt)[1:nX]# get taxa names
  DQT = quartetTableDominant(rqt, bigweights = "finite")
  roughtree <-
    root(WQDS(DQT,method=method), outgroup = taxanames[1], resolve.root = TRUE)# construct initial estimate of tree, with arbitrary root
  
  internaledgeindices = which(roughtree$edge[, 2] > Ntip(roughtree)) # determine internal edges
  length = max(roughtree$edge.length[internaledgeindices]) #get longest internal edge length
  if (length <= stopAt) {
    return(unroot(roughtree))
  }
  
  longestedges = which(roughtree$edge.length[internaledgeindices] == length)# and longest edges
  if (length(longestedges) > 1) {
    longestedge = sample(longestedges, 1)# random tie breaking
  } else {
    longestedge = longestedges
  }
  
  longestinternaledgeindex = internaledgeindices[longestedge]# ...to pick edge
  node = roughtree$edge[longestinternaledgeindex, 2] # get child of edge
  taxaAind = phangorn::Descendants(roughtree, node, "tips")[[1]] # determine split sets for that edge
  taxaA = roughtree$tip.label[taxaAind]
  nA = length(taxaAind)
  taxaB = setdiff(taxanames, taxaA)
  nB = nX - nA
  
  if (nA > 2) {
    # at least 3 taxa in A
    RQTA = quartetTableCollapse(rqt, taxaA, taxaB) #create smaller tables for set A \cup {"setB"}
    Ap1tree = WQDCrecursive(RQTA) #recursively build tree
    Bname = paste0(sort(taxaB), collapse = "")
    Atree = drop.tip(root(Ap1tree, Bname, resolve.root = TRUE), Bname)
  } else {
    Atree = drop.tip(root(roughtree, taxaB[1], resolve.root = TRUE), taxaB) #keep only cherry for A
  }
  
  
  if (nB > 2) {
    # at least 3 taxa in B
    RQTB = quartetTableCollapse(rqt, taxaB, taxaA) #create smaller tables for set B \cup {"setA"}
    Bp1tree = WQDCrecursive(RQTB) #recursively build tree
    Aname = paste0(sort(taxaA), collapse = "")
    Btree = drop.tip(root(Bp1tree, Aname, resolve.root = TRUE), Aname)
  } else {
    # only 2 taxa in B
    Btree = drop.tip(root(roughtree, taxaA[1], resolve.root = TRUE), taxaA) #keep only cherry for A
  }
  
  Atree$root.edge = length
  Xtree = unroot(bind.tree(Atree, Btree, position = length))
  term_edges = which(Xtree$edge[, 2] <= length(taxanames)) #modify terminal edge lengths
  Xtree$edge.length[term_edges] = terminal
  return(Xtree)
}


#' Reduce quartet table by combining some taxa 
#'
#' Form a smaller resolved quartet table by lumping some taxa into a composite taxon.
#'
#' @details This function is needed for the recursive calls in \code{WQDSrec}.
#' It should only be applied to a resolved quartet table which includes counts
#' for all possible quartets on the taxa (though counts can be zero).
#' 
#' The sets \code{taxaA} and \code{taxaB} must be disjoint. Their union need not be all taxa in \code{rqt}.
#'
#' @param rqt a resolved quartet table, as from \code{quartetTableResolved}
#' @param taxaA a vector of taxon names in \code{rqt} to be included in the output table
#' @param taxaB a vector of taxon names in \code{rqt} to form new composite taxon in the output table
#'
#' @return a resolved quartet table (of class matrix) with \code{length(taxaA)+1} taxa; the
#'  composite taxon is named as the concatenation of the sorted
#'  names in \code{taxaB}
#'
#' @seealso \code{\link{WQDCrecursive}}
#'
#' @export
quartetTableCollapse = function(rqt, 
                                taxaA, 
                                taxaB) {
  taxaAB = union(taxaA, taxaB)
  if (!setequal(intersect(taxaAB, colnames(rqt)), taxaAB)) {
    stop("Not all taxa in given sets in table.")
  } else {
    if (length(intersect(taxaA, taxaB)) != 0)
      stop("Given taxon sets must be disjoont.")
  }
  
  M = dim(rqt)[1]
  ntaxa = dim(rqt)[2] - 3
  coltaxa = colnames(rqt)[1:ntaxa]
  nA = length(taxaA)
  nB = length(taxaB)
  taxaA = coltaxa[sort(match(taxaA, coltaxa))] #order taxa A as in rqt
  
  #set up table for subset A + "setB"
  nAp1 = nA + 1
  mRQTA = choose(nAp1, 4)
  RQTA = matrix(0, mRQTA, nAp1 + 3)
  qnames = c("12|34", "13|24", "14|23")
  Bname = paste0(sort(taxaB), collapse = "") #create name for composite taxon for set B
  colnames(RQTA) = c(taxaA, Bname, qnames) #create column names
  
  # encode 4-taxon sets in table
  m = 0
  for (i in 1:(nAp1 - 3)) {
    # for each 4-taxon set
    for (j in (i + 1):(nAp1 - 2)) {
      for (k in (j + 1):(nAp1 - 1)) {
        for (l in (k + 1):nAp1) {
          m = m + 1
          RQTA[m, c(i, j, k, l)] = 1 #encode set
        }
      }
    }
  }
  
  # fill counts in new table, by running through old table
  for (m in 1:M) {
    memA = taxaA[which(rqt[m, taxaA] == 1)]
    numA = length(memA)
    
    if (numA == 4) {
      mA = which(rowSums(RQTA[, memA, drop = FALSE]) == 4)
      RQTA[mA, qnames] = rqt[m, qnames]
    } else {
      if (numA == 3)  {
        mA = which(rowSums(RQTA[, c(memA, Bname), drop = FALSE]) == 4)
        taxonB = which(is.element(coltaxa[which(rqt[m, coltaxa] == 1)], taxaB))
        z = rqt[m, qnames]
        if (taxonB == 1) {
          RQTA[mA, qnames] = RQTA[mA, qnames] + z[c(3, 2, 1)]
        } else{
          if (taxonB == 2) {
            RQTA[mA, qnames] = RQTA[mA, qnames] + z[c(2, 3, 1)]
          } else{
            if (taxonB == 3) {
              RQTA[mA, qnames] = RQTA[mA, qnames] + z[c(1, 3, 2)]
            } else{
              RQTA[mA, qnames] = RQTA[mA, qnames] + z
            }
          }
        }
      }
    }
  }
  return(RQTA)
}

#' Estimate edge lengths on a species tree from gene tree quartet counts
#'
#' Estimate edge lengths, in coalescent units, on an unrooted species tree 
#' from a table of resolved quartet counts from a collection of gene trees.
#'
#' @details While the argument \code{tree} may be supplied as rooted or unrooted, metric or topological,
#' only its unrooted topology will be used for determining new metric estimates.
#' 
#' Counts of quartets for all those quartets which define a single edge
#' on the tree (i.e., whose internal edge is the
#' single edge on the unrooted input tree) are summed, and from this an
#' estimate of the branch length is computed. If \code{method= "simpleML"} this is the maximum likelihood estimate.
#' If \code{method="simpleBayes"} this is the Bayesian estimate of Theorem 2
#' of \insertCite{SayMir16;textual}{MSCquartets}, using parameter \code{lambda}.
#' Using \code{lambda=1/2}
#' gives a flat prior on [1/3,1] for the probability of the quartet displayed on the species tree.
#' 
#' These methods are referred to as `simple' since they use only the quartets defining a single edge of the species tree.
#' Quartets with central edges composed of several edges in the species tree are ignored.
#'
#' Note that branch length estimates may be 0 (if the count for the quartet
#' displayed on the input tree is not dominant),
#' positive, or \code{Inf} (if
#' the counts for quartet topologies not displayed on the input tree are all 0, and \code{method="simpleML"}).
#'
#' @references
#'  \insertRef{SayMir16}{MSCquartets}
#'
#' @param tree a phylo object, giving a resolved tree on which to estimate edge lengths
#' @param rqt a resolved quartet table, as from \code{quartetTableResolved},
#'  in which all taxa on \code{tree} appear
#' @param terminal an edge length to assign to terminal edges, whose lengths cannot be estimated
#' @param method \code{"simpleML"} or \code{"simpleBayes"}
#' @param lambda a positive parameter for the \code{"simpleBayes"} method
#'
#' @return  an unrooted metric tree with the same topology as \code{tree}, of type phylo
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' taxanames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,taxanames[1:6])
#' RQT=quartetTableResolved(QT)
#' DQT=quartetTableDominant(RQT)
#' tree=QDS(DQT)
#' write.tree(tree)
#' plot(tree)
#' metricMTree=estimateEdgeLengths(tree,RQT,method="simpleML")
#' write.tree(metricMTree)
#' plot(metricMTree)
#' metricBTree=estimateEdgeLengths(tree,RQT,method="simpleBayes")
#' write.tree(metricBTree)
#' plot(metricBTree)
#' 
#' @importFrom phangorn Descendants
#' @export
estimateEdgeLengths = function(tree,
                               rqt,
                               terminal = 1,
                               method = "simpleML",
                               lambda = 1/2) {
  if (!(method %in% c("simpleML","simpleBayes"))) {
    stop('Argument method must be "simpleML" or "simpleBayes".')
  }
  
  if (method=="simpleBayes") {
    if (lambda<=0) {
      stop("Argument lambda must be positive.")
    }
  } else { # for "simpleML"
      lambda=0 
    }
  
  newroot = tree$tip.label[1]
  tree = root(tree, outgroup = newroot, resolve.root = TRUE)
  
  nedges = dim(tree$edge)[1]
  tree$edge.length = rep(0, nedges)
  taxanames = tree$tip.label
  ntaxa = length(taxanames)
  
  tableTaxa = colnames(rqt)[1:ntaxa]
  qlabels = c("12|34", "13|24", "14|23")
  
  for (i in 1:nedges) {
    if (tree$edge[i, 2] <= ntaxa) {
      tree$edge.length[i] = terminal # pendant edge
    } else {
      edgeParent = tree$edge[i, 1]
      if (edgeParent == (ntaxa + 1)) {
        tree$edge.length[i] = 0 #root edge
      } else {
        edgeChild = tree$edge[i, 2]
        edgesBelowChild = which(tree$edge[, 1] == edgeChild)
        setA = phangorn::Descendants(tree, tree$edge[edgesBelowChild[1], 2], "tips")[[1]]
        setA = taxanames[setA]
        setB = phangorn::Descendants(tree, tree$edge[edgesBelowChild[2], 2], "tips")[[1]]
        setB = taxanames[setB]
        edgesBelowParent = which(tree$edge[, 1] == edgeParent)
        if (length(edgesBelowParent) != 2) {
          stop("Something's wrong...")
        } else {
          otherEdgeBelowParent = setdiff(edgesBelowParent, i)
          setC = phangorn::Descendants(tree, tree$edge[otherEdgeBelowParent, 2], "tips")[[1]]
          setC = taxanames[setC]
          setD = setdiff(taxanames, union(setA, union(setB, setC)))
          CF = rep(0, 3)
          
          for (a in setA) {
            for (b in setB) {
              for (c in setC) {
                for (d in setD) {
                  kk = 1
                  k = 0
                  while (k == 0) {
                    # find row for this quartet
                    taxind = rqt[kk, a] * rqt[kk, b] * rqt[kk, c] * rqt[kk, d]
                    if (taxind == 1) {
                      k = kk
                    }
                    kk = kk + 1
                  }
                  
                  
                  acol = which(tableTaxa == a)
                  bcol = which(tableTaxa == b)
                  ccol = which(tableTaxa == c)
                  dcol = which(tableTaxa == d)
                  r = rank(c(acol, bcol, ccol, dcol))
                 
                  if (identical(r, c(1, 2, 3, 4)) |
                      identical(r, c(2, 1, 4, 3)) |
                      identical(r, c(3, 4, 1, 2)) |
                      identical(r, c(4, 3, 2, 1))) {
                    qorder = c(1, 2, 3)
                  } else {
                    if (identical(r, c(1, 2, 4, 3)) |
                        identical(r, c(2, 1, 3, 4)) |
                        identical(r, c(3, 4, 2, 1)) |
                        identical(r, c(4, 3, 1, 2))) {
                      qorder = c(1, 3, 2)
                    } else {
                      if (identical(r, c(1, 3, 2, 4)) |
                          identical(r, c(2, 4, 1, 3)) |
                          identical(r, c(3, 1, 4, 2)) |
                          identical(r, c(4, 2, 3, 1))) {
                        qorder = c(2, 1, 3)
                      } else{
                        if (identical(r, c(1, 3, 4, 2)) |
                            identical(r, c(2, 4, 3, 1)) |
                            identical(r, c(3, 1, 2, 4)) |
                            identical(r, c(4, 2, 1, 3))) {
                          qorder = c(2, 3, 1)
                        } else {
                          if (identical(r, c(1, 4, 2, 3)) |
                              identical(r, c(2, 3, 1, 4)) |
                              identical(r, c(3, 2, 4, 1)) |
                              identical(r, c(4, 1, 3, 2))) {
                            qorder = c(3, 1, 2)
                          } else {
                            if (identical(r, c(1, 4, 3, 2)) |
                                identical(r, c(2, 3, 4, 1)) |
                                identical(r, c(3, 2, 1, 4)) |
                                identical(r, c(4, 1, 2, 3))) {
                              qorder = c(3, 2, 1)
                            }
                          }
                        }
                      }
                    }
                  }
                  CF = CF + rqt[k, qlabels[qorder]]
                }
              }
            }
          }
          
          avCF = CF / (length(setA) * length(setB) * length(setC) * length(setD))
          f = avCF[1] / (sum(avCF) + 2 * lambda)
          if (f >= 1) {
            edgeweight = Inf
          } else {
            if (f >= (1 / 3)) {
              edgeweight = -log((1.5 * (1 - f)))
            } else {
              edgeweight = 0
            }
          }
          tree$edge.length[i] = edgeweight
        }
      }
    }
  }
  tree = unroot(tree)
  return(tree)
}

