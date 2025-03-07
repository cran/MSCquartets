#' Suppress Network Binary Nodes
#'
#' Suppress any non-root binary nodes in a phylogenetic network.
#'
#' @details
#' This function is similar to ape's \code{collapse.singles} which only works
#' on phylo objects that are trees.
#'
#'@param net a phylogenetic network, of class "evonet"
#'@param addlength if TRUE (default), add lengths of two incident edges for new edge.
#'
#'@export
suppressNetBinaryNodes = function(net, addlength = TRUE) {
  ntaxa = length(net$tip.label)
  edge = net$edge
  reticEdge = net$reticulation
  edgeLength = net$edge.length
  nIntNode = net$Nnode #number of internal nodes
  nodeLabels = net$node.label

  nodeDegrees=tabulate(c(unlist(list(edge)),unlist(list(reticEdge))))
  #nodeDegrees = evonetdegree(net, details=TRUE)

  binNodes = which(nodeDegrees == 2)
  binNodes = binNodes[which(binNodes > (ntaxa + 1))]# degree 2 nodes other than root

  for (node in binNodes) {
    parentEdge = which(edge[, 2] == node)
    childEdge = which(edge[, 1] == node)
    childNode = edge[childEdge, 2]
    edge[parentEdge, 2] = childNode
    edge = edge[-childEdge, ]
    if (addlength == TRUE) {
      edgeLength[parentEdge] = edgeLength[parentEdge] + edgeLength[childEdge]
    }
    edgeLength = edgeLength[-childEdge]
  }
  nIntNode=nIntNode-length(binNodes)
  nodeLabels=nodeLabels[-(binNodes-ntaxa)]


  # renumber nodes
  newNodeNums=1:(ntaxa+nIntNode)
  oldNodeNums=sort(unique(unlist(as.list(edge))))
  for (i in 1:2) {
    for (j in 1:dim(edge)[1]) {edge[j,i]=newNodeNums[which(oldNodeNums==edge[j,i])]}
    for (j in 1:dim(reticEdge)[1]) {reticEdge[j,i]=newNodeNums[which(oldNodeNums==reticEdge[j,i])]}
  }

  #prepare output
  net$edge = edge
  net$edge.length = edgeLength
  net$Nnode = nIntNode
  net$node.label = nodeLabels
  net$reticulation = reticEdge

  return(net)
}

evonetdegree=function(net,details=TRUE){
  return(degree(net,details))
}

