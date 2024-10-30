#'Label internal nodes on tree
#'
#'Label all internal nodes of tree, as "Node NN" where NN is the node number, and
#' plot tree with labels
#'
#'@param tree a rooted or unrooted tree (phylo object)
#'@param plot TRUE for plot, FALSE for no plot
#'@param type plot type (e.g.,"unrooted") to be passed to ape plot command
#'
#'@return a phylo object with internal node labels
#'
#'@examples
#' data(pTableYeastRokas)
#' out=TINNIK(pTableYeastRokas,test="T3",alpha=.01, beta=.05)
#' labelIntNodes(out$ToB)
#'
#' @seealso \code{\link{resolveCycle}}, \code{\link{combineCycleResolutions}},
#'   \code{\link{resolveLevel1}}
#'
#'@export
#'
labelIntNodes = function(tree, plot = TRUE, type = "unrooted")
{
  #Label internal nodes and plot ToB
  treeLab = tree
  treeLab$node.label = paste0("Node ", 1:treeLab$Nnode + length(treeLab$tip.label))

  if (plot > 0)
  {
    plot(
      treeLab,
      show.node.label = T,
      type = "unrooted",
      main = "Internal Node Labels"
    ) #plot tree with node numbers
  }
  return(treeLab)
}



####################################################
#' Resolve a node on a tree of blobs as a cycle
#'
#' Given a Tree of Blobs and quartet Concordance Factor data, resolve a specific
#' polytomy to a cycle. Resolution is performed by finding a least-squares
#' best-fit of an empirical distance to an expected distance related to the cycle,
#' as described in \insertCite{ABRW24;textual}{MSCquartets}.
#'
#' Possible distances to use are the NANUQ distance and a modified NANUQ distance
#' that weights quartet trees differently, but from which the cycle structure is
#' still identifiable.
#'
#' For multifucations of degree less than a designated cutoff, all possible
#' circular orders and choices of hybrid nodes are considered in choosing the
#' best. Above that cutoff, a heuristic method built on the modified NANUQ
#' distance is used to obtain a small number of orders likely to be good fits,
#' with the least-squares fitting applied only to those.
#'
#'
#'@references
#' \insertRef{ABRW24}{MSCquartets}
#'
#'@param ToB an unrooted tree of blobs (phylo object), as determined by TINNIK
#'or another method, with multifurcations labelled by \code{labelIntNodes}
#'@param node number of an internal node to be resolved
#'@param pTable a table of qcCFs, with columns p_star and p_test
#'@param test either "T3" or "cut", indicating the test to use for determining what
#'qcCFs indicate hybridization
#'@param alpha test level for p_test
#'@param beta test level for p_star
#'@param distance cycle resolution distance to be used; one of "NANUQ" or "modNANUQ"
#'@param hdegree resolve a multifurcation of this degree or larger by a heuristic
#'method; must be at least 5
#'@param plot if FALSE (0), no plots; if TRUE (>0), make plots of resolved cycle(s) considered
#'best and histogram of measure of fit for all hybrid/orders considered
#'@param delta cutoff for relative difference in squared residuals and smallest,
#'(RSS-minRSS)/minRSS, for determining near ties as "best" fit resolutions
#'
#'@return a list of resolution information, given as a list of:
#'\itemize{
#'\item \code{$node} node number,
#'\item \code{$cycleRes} list [[1]]-[[k]] of  best resolutions,
#'\item \code{$RSSs} RSSs from all cycle resolutions considered in choosing best.
#'}
#'Each resolution is itself a 5-element list with entries:
#'\itemize{
#'\item \code{$cycleNet} Newick network with 1 cycle (with all edge lengths 1)
#'\item \code{$cycleRSS} RSS for cycle,
#'\item \code{$taxonGroups} taxon groups for cycle,
#'\item \code{$order} order of groups around cycle,
#'\item \code{$nonRootEdges} logical vector indicating edges of \code{ToB}
#' where root cannot be.
#' }
#' (Items \code{$taxonGroups,$order,$nonRootEdges} are
#' needed to combine resolutions to form networks with multiple cycles using
#' \code{\link{combineCycleResolutions}}, and otherwise may not be of interest to users).
#'
#'
#'@examples
#' data(pTableYeastRokas)
#' out=TINNIK(pTableYeastRokas, alpha=.01, beta=.05)
#' ToB=labelIntNodes(out$ToB)
#' resolveCycle(ToB, node=9, pTable=out$pTable, alpha=.01, beta=.05, distance="NANUQ")
#'
#'@seealso \code{\link{TINNIK}}, \code{\link{labelIntNodes}}, \code{\link{combineCycleResolutions}},
#' \code{\link{resolveLevel1}}
#'
#'@importFrom ape keep.tip evonet write.evonet root bind.tree
#'@importFrom phangorn mrca.phylo Descendants
#'@importFrom graphics hist
#'
#'@export
resolveCycle = function(ToB,
                        node,
                        pTable,
                        test = "T3",
                        alpha,
                        beta,
                        distance = "NANUQ",
                        hdegree = 10,
                        plot = TRUE,
                        delta = 10 ^ -6)
{
  if (is.rooted(ToB)) {
    stop("Argument 'ToB' must be unrooted.")
  }
  if (is.null(ToB$node.label)) {
      stop("Argument 'ToB' must have internal nodes labelled by 'labelIntNodes'.")
    }
  if (!(distance %in% c("NANUQ", "modNANUQ"))) {
    stop("Argument 'distance' is not valid; must be 'NANUQ' or 'modNANUQ'.")
  }

  #Initialize
  nodeResolutions = list(node = node)# storage for Newick and order info for this node's resolutions
  cycleRes = list()# storage for all best cycle resolutions of this node we find

  taxa = ToB$tip.label
  ntaxa = length(taxa)

  groups = nodeGroups(ToB, node) #Determine groups of taxa around node to be resolved
  nGroups = length(groups) #degree of multifurcation
  if (nGroups <= 3) {
    stop(paste0("Node degree is ", nGroups, "; no need to resolve."))
  }
  hdegree = max(c(hdegree, 5)) #Heuristic not used for 4-blobs
  heurflag = (nGroups >= hdegree) # set flag, TRUE means use heuristic method

  outgroup = groups[[nGroups]] # make a copy of tree rerooted so last group is the outgroup
  ToBr = ape::root(ToB, outgroup = outgroup)

  #convert groups to a vector encoding
  groupvec = rep(0, ntaxa) #make it easy to look up group by taxon name
  names(groupvec) = taxa
  for (j in 1:nGroups) {
    groupvec[groups[[j]]] = j # put group number in list by taxon name
  }

  if (heurflag) {
    D = blobDistance(pTable,
                     taxa,
                     groupvec,
                     test = 'T3',
                     alpha,
                     beta,
                     dist = "modNANUQ") #compute modNANUQ distance for these blob groups

    orders = ordersHeuristicmodNANUQ(D, delta) #  use Min modNANUQ distance approach to find orders
  } else {
    orders = circHybOrders(nGroups) # consider all orders
  }

  #Find measures of fit
  D = blobDistance(pTable,
                   taxa,
                   groupvec,
                   test = 'T3',
                   alpha,
                   beta,
                   dist = distance) #compute empirical distance for these blob groups

  if (distance == 'NANUQ') {
    #Find expected sunlet distance for the distance
    E = expNANUQCycleDist(dim(D)[1])
  } else {
    E = expmodNANUQCycleDist(dim(D)[1])
  }


  fitMeasure = fitCycleOrders(D, E, orders) # measure fit with orders

  minfit = min(fitMeasure) #find best fit measure
  bestorders = which((fitMeasure - minfit) <= minfit * delta)

  if (length(bestorders) > 1) {
    message(paste0("For node ", node, ", ", length(bestorders)," resolutions found."))
    bestorders = bestorders[order(fitMeasure[bestorders])] #sort so best are first
  }

  for (k in bestorders) {
    # consider all best fit orders
    nonRootEdges = rep(F, dim(ToB$edge)[1]) #indicator vector for edges we know can't contain root

    order = orders[k, ]
    RSS = fitMeasure[k]
    if (nGroups > 4) {
      # for large cycles, flag edges below hybrid
      hybGroup = which(groupvec == match(1, order))
      if (length(hybGroup) == 1) {
        hybDesc = hybGroup
      }
      else{
        #first make sure ToB root is outside of hybrid group subtree
        mrca = phangorn::mrca.phylo(ToB, hybGroup)
        hybDesc = c(mrca, phangorn::Descendants(ToB, mrca, "all"))
      }

      if (length(hybDesc) != ToB$Nnode + ntaxa)
        #if (arbitrary) root of ToB was not in hybrid group subtree
      {
        nonRootEdges[match(hybDesc, ToB$edge[, 2])] = TRUE # flag edges ending at these nodes as below a hybrid
      } else
        #root of ToB is in hybrid group subtree (due to unfortunate rooting of ToB)
      {
        hyfitMeasureonDesc = phangorn::Descendants(ToB, node, "all")
        nonRootEdges[which(!(ToB$edge[, 2] %in% hyfitMeasureonDesc))] = TRUE # flag all edges in hybrid group subtree as below a hybrid
      }
    }

    #generate tree that will represent cycle once reticulation edge added
    groupt = read.tree(text = paste0("(", -nGroups, ",-1);")) # generate phylo object, mostly to be written over, using negative numbers as taxa names
    groupt$edge = matrix(0, nGroups + 2, 2)
    groupt$edge[, 1] = c(3, 3:(3 + nGroups))
    groupt$edge[, 2] = c(1, 4:(3 + nGroups), 2)
    groupt$edge.length = rep(1, nGroups + 2)
    groupt$Nnode = nGroups + 1
    groupt$node.label = c("r", "TH1", -(nGroups - 1):-2, "#H1")

    for (i in 1:nGroups) {
      #attach subtrees for each group
      gtaxa = which(groupvec == i)
      subtree = ape::keep.tip(ToBr, gtaxa)
      subtree$root.edge = 1 #add a root edge for resolution when attached
      attachlabel = -order[i]# attach at this label in cycle tree just created
      if (attachlabel %in% c(-1, -nGroups)) {
        attachPoint = which(groupt$tip.label == attachlabel) #will attach at a tip
      } else {
        attachPoint = length(groupt$tip.label) + which(groupt$node.label == attachlabel) #will attach at an internal node
      }
      groupt = ape::bind.tree(groupt, subtree, where = attachPoint)# attach the subtree
    }

    # add reticulation edge, and clean up internal node labels, edge lengths
    rt = which(groupt$node.label == "rt")
    H1 = which(groupt$node.label == "#H1")
    TH1 = which(groupt$node.label == "TH1")
    groupt$node.label = rep("", length(groupt$node.label))
    groupt$node.label[c(H1, TH1)] = c("#H1", "TH1")
    nTips = length(groupt$tip.label)
    groupn = ape::evonet(groupt, from = nTips + TH1, to = nTips + H1)
    groupn$edge.length = rep(1, length(groupn$edge.length))

    if (plot > 0)
    {
      titletext = bquote("distance: " * .(distance) * "," * ~ alpha * "=" * .(alpha) * "," ~ beta * "=" * .(beta))

      subtitle = "Network should be semidirected; Rooting is arbitrary."
      if (nGroups == 4)
        subtitle = "Network should be semidirected; Rooting is arbitrary.\n Hybrid node on 4-cycle is not identified."
      plot(
        groupn,
        show.node.label = T,
        type = "unrooted",
        main = paste0(
          "Resolution ",
          which(bestorders == k),
          " of Node ",
          node,
          "; RSS=",
          signif(RSS, digits = 6)
        ),
        sub = subtitle
      )
      mtext(eval(bquote(.(titletext))),
            side = 3,
            line = 0,
            cex = 1) # add rest of title text
    }

    netnwk = ape::write.evonet(groupn)
    netnwk=gsub(":-?[0123456789]+",":1",netnwk) #replace ape's odd branch lengths by 1
    cycleRes = append(cycleRes, list(
      list(
        cycleNet = netnwk,
        cycleRSS = fitMeasure[k],
        taxonGroups = groupvec,
        order = order,
        nonRootEdges = nonRootEdges
      )
    ))
    #save this resolution and all info
  }

  nodeResolutions = list(node = node,
                         cycleRes = cycleRes,
                         RSSs = list(as.double(fitMeasure))) #save for output

  if (plot > 0) {
    #
    graphics::hist(fitMeasure,
                   main = paste0(
                     "Node ",
                     node,
                     "; RSSs for ",
                     length(fitMeasure),
                     " orders considered"
                   ))
    mtext(eval(bquote(.(titletext))),
          side = 3,
          line = 0,
          cex = 1) # add rest of title text
  }
  return(nodeResolutions)
}

####################################################
#' Resolve Tree of Blobs to Level-1 network
#'
#' Given a Tree of Blobs and qcCF information, resolve all multifurcations to
#' cycles. Resolution is performed by finding a least-squares best-fit of an
#' empirical distance to an expected distance related to the cycle, as described
#' in \insertCite{ABRW24;textual}{MSCquartets}.
#'
#' Possible distances to use are the NANUQ distance and a modified NANUQ distance
#' that weights quartet trees differently, but from which the cycle structure is
#' still identifiable.
#'
#' For multifucations of degree less than a designated cutoff, all possible
#' circular orders and choices of hybrid nodes are considered in choosing the
#' best. Above that cutoff, a heuristic method  is used to obtain
#' a small number of orders likely to be good fits, with the least-squares
#' fitting applied only to those.
#'
#'@references
#' \insertRef{ABRW24}{MSCquartets}
#'
#'@param ToB an unrooted tree of blobs (phylo object) as determined by TINNIK or
#'another method
#'@param pTable a table of qcCFs, with columns p_star and p_test
#'@param test either "T3" or "cut", indicating test to use for determining what
#'qcCFs indicate hybridization
#'@param alpha test level for p_test
#'@param beta test level value for p_star
#'@param distance cycle resolution distance to be used ("NANUQ" or "modNANUQ")
#'@param hdegree resolve a multifurcation of this degree or larger by a heuristic
#'method; must be at least 5
#'@param plot if 0, no plots; if 1, plot only possible root locations on ToB and
#'full resolution; if 2, include plots of each individual blob resolution, if 3
#'include histograms of measure of fit for all hybrid/orders considered in choosing
#'best
#'@param delta cutoff for relative difference in squared residuals and smallest,
#'(RSS-minRSS)/minRSS, for determining near ties as "best" fit resolutions
#'@param fullResMax maximum number of full resolutions (all multifurcations at once)
#' to form; if the product of the number of resolutions of individual multifurcations
#' exceeds this, no full resolutions are produced, although \code{\link{combineCycleResolutions}}
#' can be applied to produce them.
#'
#'
#'@return a list of resolutions and squared residuals:
#'\itemize{
#'\item [[1]] is a list of Newick
#'resolutions of entire network, with all edge lengths 1 (NULL if one cannot be produced or \code{fullResMax} is exceeded),
#'\item [[2]]-[[n]] are individual resolutions of each  multifurcation on \code{ToB},
#'each given as a list as output from \code{\link{resolveCycle}}.
#'}
#'
#'@examples
#' data(pTableYeastRokas)
#' out=TINNIK(pTableYeastRokas, alpha=.01, beta=.05)
#' ToB=labelIntNodes(out$ToB)
#' resolveLevel1(ToB, pTable=out$pTable, alpha=.01, beta=.05, distance="NANUQ")
#'
#'@seealso \code{\link{TINNIK}}, \code{\link{labelIntNodes}}, \code{\link{resolveCycle}},
#' \code{\link{combineCycleResolutions}}
#'
#'@importFrom ape keep.tip evonet write.evonet root bind.tree
#'@importFrom phangorn mrca.phylo Descendants
#'@importFrom graphics hist
#'
#'@export
resolveLevel1 = function(ToB,
                         pTable,
                         test = "T3",
                         alpha,
                         beta,
                         distance = "NANUQ",
                         hdegree = 10,
                         plot = 2,
                         delta = 10 ^ -6,
                         fullResMax = 10)
{
  if (is.rooted(ToB)) {
    stop("Argument 'ToB' must be unrooted.")
  }

  #Label internal nodes and plot ToB
  ToBLab = labelIntNodes(ToB, plot, type = "unrooted") #Label internal nodes, and plot

  taxa = ToB$tip.label
  ntaxa = length(taxa)
  resolutions = c() #list for all resolutions

  #Determine multifurcations on ToB to be resolved into cycles
  polyList = c()
  if (length(Descendants(ToB, ntaxa + 1, 'children')) > 3) {
    polyList = c(polyList, ntaxa + 1)
  }
  if (ToB$Nnode >= 2) {
    for (i in (ntaxa + 2):(ntaxa + ToB$Nnode)) {
      #loop over internal nodes
      if (length(Descendants(ToB, i, 'children')) > 2) {
        polyList = c(polyList, i)
      }
    }
  }
  npoly = length(polyList)

  if (npoly == 0) {
    message("No multifurcations on tree.")
  } else {
    if (npoly == 1) {
      message("1 multifurcation on tree, at node ",
              as.character(polyList))
    } else {
      message(
        npoly,
        " multifurcations on tree, at nodes ",
        paste(as.character(polyList), collapse = ", "),
        "."
      )
    }

    #Initialize
    nonRootEdges = rep(FALSE, dim(ToB$edge)[1]) #indicator vector for edges that can't contain root
    #resolveAll = TRUE # flag to produce one resolved network at end

    # Resolve each multifurcation
    numFullRes = 1 #compute number of Full resolutions
    for (node in polyList) {
      resCyc = resolveCycle(ToBLab,
                            node,
                            pTable,
                            test,
                            alpha,
                            beta,
                            distance,
                            hdegree,
                            plot - 1,
                            delta)
      numFullRes = numFullRes * length(resCyc$cycleRes)
      resolutions = append(resolutions, list(resCyc)) #add these resolutions to list for all polytmoies
    }
    # combine resolutions of polytomies
    if (numFullRes > fullResMax)
    {
      warning("Exceeded fullResMax; use 'combineCycleResolutions' to see full resolutions")
      nets = list("fullResMax exceeded; full resolutions not produced")
    }
    else {
      titletext = bquote("distance: " * .(distance) * "," * ~ alpha * "=" * .(alpha) * "," ~ beta * "=" * .(beta))
      nets = (combineCycleResolutions(ToBLab, resolutions, plot, titletext))
    }
    resolutions = append(list(nets), resolutions)
  }
  return(resolutions)
}


#############################################################
#' Combine several cycle resolutions on a tree of blobs to create a network
#'
#' Given a list of resolutions of different multifurcations on a tree of blobs
#' (each as produced by \code{\link{resolveCycle}}), combine these with the tree
#' of blobs to form a network.
#'
#' This function is useful for forming near-optimal networks when there are
#' several resolutions that have similar fit for some of the multifurcations.
#'
#' @param ToB an unrooted tree of blobs for the network, with multifurcating nodes
#' labelled by \code{\link{labelIntNodes}}
#' @param resolutions a list of resolutions (each of which may be a list)
#' for different nodes with elements in format described in output of \code{\link{resolveCycle}}
#' @param plot if FALSE (0), no plots; if TRUE (>0) plot networks
#' @param titletext a string of text for plot
#'
#' @return a list of Newick strings for the networks, with all edge lengths 1
#'
#' @examples
#' data(pTableYeastRokas)
#' out=TINNIK(pTableYeastRokas, alpha=.01, beta=.05)
#' ToB=labelIntNodes(out$ToB)
#' R9=resolveCycle(ToB, node=9, pTable=out$pTable, alpha=.01, beta=.05, distance="NANUQ")
#' R10=resolveCycle(ToB, node=10, pTable=out$pTable, alpha=.01, beta=.05, distance="NANUQ")
#' combineCycleResolutions(ToB, resolutions=list(R9,R10),plot=TRUE)
#'
#'@seealso \code{\link{TINNIK}}, \code{\link{labelIntNodes}}, \code{\link{resolveCycle}},
#'  \code{\link{resolveLevel1}}
#'
#' @export
combineCycleResolutions = function(ToB,
                                   resolutions,
                                   plot = 1,
                                   titletext = NULL)
{
  if (is.null(ToB$node.label)) {
    stop("Argument 'ToB' must have internal nodes labelled by 'labelIntNodes'.")
  }
  savedToB = ToB
  nEdgeToB = dim(ToB$edge)[1] #save number of edges on original ToB

  #Check for duplicate nodes among resolutions, 4-cycles
  cycle4 = FALSE
  nodes = c()
  nRes = length(resolutions) # number of nodes resolved
  resols = vector(mode = "list", length = nRes)
  for (i in 1:nRes)
    #check for repeated nodes and expand resolution lists
  {
    if (resolutions[[i]]$node %in% nodes) {
      stop("Input error: Node is repeated in given resolutions.")
    }
    else {
      nodes = c(nodes, resolutions[[i]]$node)
    }

    if (length(resolutions[[i]]$cycleRes[[1]]$order) == 4) {
      cycle4 = TRUE
    }

    # expand resolutions of node into list of individual resolutions
    curNodeRes = resolutions[[i]] # focus on resolutions of single node
    k = length(curNodeRes$cycleRes) # number of resolutions of this node
    curNodeResList = vector(mode = "list", length = k) # list for these
    for (j in 1:k) {
      curNodeResList[[j]] = list(
        node = resolutions[[i]]$node,
        cycleRes = resolutions[[i]]$cycleRes[[j]],
        RSSs = resolutions[[i]]$RSSs
      )
    }
    resols[[i]] = curNodeResList
  }
  if (nRes==1) { #if only 1 node resolved, expand.grid doesn't work}
    numFullRes=length(resols[[1]])
    allFullRes=array(resols[[1]],c(numFullRes,1))
  } else {
   allFullRes = do.call(expand.grid, resols) #create big table with rows giving all full resolutions
   numFullRes = nrow(allFullRes)
  }

  Nets = vector(mode = "list", length = numFullRes)
  for (ri in 1:numFullRes) {
    # loop through the full resolutions

    resols = allFullRes[ri, ] #pick up infor for one full resolultion
    if (length(resols)==1) { resols=list(resols)} #if only 1 node resolution, need to make list
    rest=resols

    nonRootEdges = rep(FALSE, nEdgeToB) #indicator vector for edges that can't contain root

    for (j in 1:length(resols))
      # check for rootability
    {
      nonRootEdges = nonRootEdges |
        resols[[j]][[1]]$cycleRes$nonRootEdges
    }

    if (all(nonRootEdges))
      #if not rootable
    {
      Nets[[ri]] = NULL
      message(paste0("No root location is consistent for resolution ", ri))
    } else
    {
      if (plot > 0)
        # plot allowable root locations
      {
        #plot possible root locations
        edgecolor = rep("black", length(nonRootEdges))
        edgecolor[which(nonRootEdges == T)] = "red"

        if (numFullRes > 1)
          maintitle = paste0("Possible root locations; resolution ", ri)
        else
          maintitle = "Possible root locations"

        plot(
          savedToB,
          edge.color = edgecolor,
          type = "unrooted",
          main = maintitle,
          sub = "May be rooted on black edges; within cycles also possible"
        )
        mtext(eval(bquote(.(titletext))),
              side = 3,
              line = 0,
              cex = 1) # add rest of title text
      }

      ToB = savedToB #restore original ToB
      pendantRootEdges = (!nonRootEdges) &
        (ToB$edge[, 2] < length(ToB$tip.label)) #find possible pendant edges for rooting
      possibleOutTaxa = ToB$edge[which(pendantRootEdges), 2] # pick taxon to be outgroup for rooting
      outTaxon = possibleOutTaxa[1]
      outTaxonLabel = ToB$tip.label[outTaxon]

      ToB = root(ToB, outgroup = outTaxon, resolve.root = T)
      ToB$edge.length = rep(1, length(ToB$edge.length))

      Net = ToB
      TaxaFromToB = Net$tip.label # need original order to interpret multifurcation data

      for (ii in 1:length(resols))
        # loop through multifurcations to resolve
      {
        res = resols[[ii]][[1]] #get info for this multifurcation resolution
        nodeLab = res$node
        groupvec = res$cycleRes$taxonGroups
        order = res$cycleRes$order
        n = length(order)

        if (n == 4) {
          # Since hybrid node not determined for 4-cycles,
          # may have to switch order if outTaxon descended from
          # hybrid group

          if (order[groupvec[outTaxon]] == 1) {
            #check if need to reorder to avoid outTaxon as hybrid descendant
            order = order[c(2:length(order), 1)] # if so, shift order
            if (order[2] > order[4])
            {
              order[2:4] = order[4:2]
            }
          }
        }



        #generate tree that will represent cycle once reticulation edge added
        groupt = read.tree(text = paste0("(", -n, ",-1);")) # generate phylo object, mostly to be written over (negatives for taxa names)
        groupt$edge = matrix(0, n + 2, 2)
        groupt$edge[, 1] = c(3, 3:(3 + n))
        groupt$edge[, 2] = c(1, 4:(3 + n), 2)
        groupt$edge.length = rep(1, n + 2)
        groupt$Nnode = n + 1
        groupt$node.label = c("r", paste0("TH", ii), -(n - 1):-2, paste0("#H", ii))

        Net = ape::root(Net,
                        node = length(Net$tip.label) + which(Net$node.label == paste0("Node ", nodeLab))) #temporarily root at multifurcation

        for (i in 1:n)
          #loop through groups
        {
          gtaxa = TaxaFromToB[which(groupvec == i)] #Get taxa in this group
          subtree = drop.tip(Net,
                             setdiff(Net$tip.label, gtaxa),
                             collapse.singles =
                               F)
          subtree$node.label[1] = ""
          subtree$root.edge = 1 #add a root edge for resolution when attached
          attachlabel = -order[i]# attach at this label in cycle tree just created

          if (attachlabel %in% c(-1, -n)) {
            attachPoint = which(groupt$tip.label == attachlabel) #will attach at a tip
          } else {
            attachPoint = length(groupt$tip.label) + which(groupt$node.label == attachlabel) #will attach at an internal node
          }
          groupt = ape::bind.tree(groupt, subtree, where = attachPoint)# attach the subtree
          groupt$node.label[which(groupt$node.label == attachlabel)] = "" #wipe out this label
        }
        Net = groupt
      }

      Net$node.label[which(!(substr(Net$node.label, 1, 2) %in% c("TH", "#H")))] =
        "" # clean up irrelevent node labels

      Net = ape::root(Net, outgroup = outTaxonLabel, resolve.root = TRUE) # root properly
      Net$node.label[which(Net$node.label == "Root")] = "" #remove label for root
      if (Net$node.label[1] != "")
        #if root is labeled (THn), need to move root. Ape requires a hack for this
      {
        newLeaf <- list(
          edge = matrix(c(2, 1), 1, 2),
          #creat edge to extra taxon
          tip.label = c("***"),
          edge.length = 1,
          Nnode = 1
        )
        class(newLeaf) <- "phylo"

        Net$root.edge = 1
        Net = bind.tree(newLeaf, Net) #Attach an extra taxon above root
        Net = drop.tip(Net, outTaxonLabel, collapse.singles = F) # drop old outtaxon
        Net$tip.label[which(Net$tip.label == "***")] = outTaxonLabel #and rename new one
        Net$node.label[1] = ""
      }

      Net$edge.length = rep(1, length(Net$edge.length))

      fromNode = match(paste0("TH", 1:length(resols)), Net$node.label) + length(Net$tip.label)
      toNode = match(paste0("#H", 1:length(resols)), Net$node.label) + length(Net$tip.label)

      Net = ape::evonet(Net, from = fromNode, to = toNode)

      if (plot > 0) {
        subtitle = "Network should be semidirected; Rooting is arbitrary"
        if (cycle4)
          subtitle = "Network should be semidirected; Rooting is arbitrary\nHybrid node on 4-cycle is arbitrary"

        if (numFullRes > 1)
          maintitle = paste0("Inferred Level-1 Network, resolution ", ri)
        else
          maintitle = "Inferred Level-1 Network"

       # plot(
       #    Net,
       #    show.node.label = T,
       #    main = maintitle,
       #    sub = subtitle
       #  )
       #  mtext(eval(bquote(.(titletext))),
       #        side = 3,
       #        line = 0,
       #        cex = 1) # add rest of title text


        plot(
          Net,
          type = "unrooted",
          show.node.label = T,
          main = maintitle,
          sub = subtitle
        )
        mtext(eval(bquote(.(titletext))),
              side = 3,
              line = 0,
              cex = 1) # add rest of title text
      }

      Nets[[ri]] = gsub(":-?[0123456789]+",":1",ape::write.evonet(Net))# write Newick for net, making sure all branch lengths are 1
    }
  }

  return(Nets)
}



###################################################
#'Compute empirical distance between taxon groups.
#'
#'From gene quartet counts, computes NANUQ or modNANUQ distances between groups
#'of taxa (which should be those around a multifurcation in a tree of blobs.
#'If these groups are not singletons, averaging is done over group elements.
#'
#'@param pTable table of giving empirical gene quartet counts for the taxa on tree, with columns p_star and p_test
#'@param taxa a list of taxon names, who positions are used in 'groups'
#'@param groupvec taxon groups encoded in vector
#'@param test to be used for detecting hybridizations in quartete ("T3" or "cut")
#'@param alpha test level for p_test
#'@param beta test level for p_star
#'@param dist the distance to compute, either "NANUQ" or "modNANUQ"
#'
#'@return the distance matrix, ordered by taxon group number
#'
#'@export
blobDistance <-
  function(pTable,
           taxa,
           groupvec,
           test = 'T3',
           alpha,
           beta,
           dist = "NANUQ") {
    if (!(is.numeric(alpha) &&
          is.numeric(beta))) {
      stop("Test levels alpha and beta must be numeric.")
    }
    if (!(dist %in% c("NANUQ", "modNANUQ"))) {
      stop('Argument "dist" must be either "NANUQ" or "modNANUQ".')
    }

    pcol = paste0("p_", test)
    output = c()

    M = dim(pTable)[1] # number of quartet concordance factors
    ntaxa = length(taxa) #number of taxa
    if (!(M == choose(ntaxa, 4))) {
      stop("Improper number of rows in pTable.")
    }

    qnames = c("12|34", "13|24", "14|23")

    nGroups = length(unique(groupvec))
    groupSizes = as.vector(table(groupvec))

    Dist = matrix(0, nGroups, nGroups) #initialize distance matrix

    for (m in 1:M) {
      # consider each set of 4 taxa
      taxanames = colnames(pTable)[which(pTable[m, 1:ntaxa] == 1)] #determine taxa
      groupNums = groupvec[taxanames]    #get groups these are in

      if (length(unique(groupNums)) == 4)  {
        #if all 4 in different groups...
        factor = 1 / prod(groupSizes[groupNums]) #compute factor to downweight dist entry by
        qcounts = pTable[m, qnames] # get counts

        if (pTable[m, "p_star"] > beta) {
          # if quartet judged as star tree
          Dist[groupNums, groupNums] = Dist[groupNums, groupNums] + 1 *
            factor #  it separates all taxa
        } else {
          qcountspr = qcounts + runif(3) * .0001 # introduce random tie breaking
          majorq = which.max(qcountspr) #determine quartet with largest count
          if (pTable[m, pcol] > alpha) {
            # if high p-value indicates tree
            if (majorq == 1) {
              Dist[groupNums[1:2], groupNums[3:4]] = Dist[groupNums[1:2], groupNums[3:4]] +
                1 * factor
              if (dist == 'modNANUQ') {
                Dist[groupNums[1], groupNums[2]] = Dist[groupNums[1], groupNums[2]] + .5 *
                  factor
                Dist[groupNums[3], groupNums[4]] = Dist[groupNums[3], groupNums[4]] +
                  .5 * factor
              }
            } else {
              if (majorq == 2) {
                Dist[groupNums[c(1, 3)], groupNums[c(2, 4)]] = Dist[groupNums[c(1, 3)], groupNums[c(2, 4)]] +
                  1 * factor
                if (dist == 'modNANUQ') {
                  Dist[groupNums[1], groupNums[3]] = Dist[groupNums[1], groupNums[3]] + .5 *
                    factor
                  Dist[groupNums[2], groupNums[4]] = Dist[groupNums[2], groupNums[4]] +
                    .5 * factor
                }
              } else {
                Dist[groupNums[c(1, 4)], groupNums[c(2, 3)]] = Dist[groupNums[c(1, 4)], groupNums[c(2, 3)]] +
                  1 * factor
                if (dist == 'modNANUQ') {
                  Dist[groupNums[1], groupNums[4]] = Dist[groupNums[1], groupNums[4]] + .5 *
                    factor
                  Dist[groupNums[2], groupNums[3]] = Dist[groupNums[2], groupNums[3]] +
                    .5 * factor
                }
              }
            }
          } else {
            # if low p-value indicates cycle
            minorq = which.min(qcountspr)
            middleq = setdiff(1:3, c(majorq, minorq))
            if (majorq == 1 || middleq == 1) {
              Dist[groupNums[1:2], groupNums[3:4]] = Dist[groupNums[1:2], groupNums[3:4]] +
                .5 * factor
            }
            if (majorq == 2 || middleq == 2) {
              Dist[groupNums[c(1, 3)], groupNums[c(2, 4)]] = Dist[groupNums[c(1, 3)], groupNums[c(2, 4)]] + .5 *
                factor
            }
            if (majorq == 3 || middleq == 3) {
              Dist[groupNums[c(1, 4)], groupNums[c(2, 3)]] = Dist[groupNums[c(1, 4)], groupNums[c(2, 3)]] + .5 *
                factor
            }
          }
        }
      }
    }
    Dist = Dist + t(Dist)
    Dist = 2 * Dist + 2 * nGroups - 4 # adjust to get appropriate distances off-diagonal
    Dist = Dist - diag(diag(Dist)) # but reset diagonal to 0

    return(Dist)
  }

########################################
#'Expected NANUQ cycle distance
#'
#'Compute expected NANUQ distance for a sunlet network with 4 or more taxa. This is used
#' to resolve multifurcations in a tree of blobs by NANUQ+ functions
#'
#'@param n number of edges in cycle
#'
#'@return an nxn distance matrix with rows/columns ordered from hybrid following circular order
#'
#'@examples
#' expNANUQCycleDist(5)
#'
#'@seealso \code{\link{resolveCycle}}
#'
#'@export
expNANUQCycleDist = function(n) {
  if (n < 4)
    stop("Cycle size argument must be >3.")
  d = matrix(0, n, n)
  c = n * (n - 5) / 2
  for (j in 2:n) {
    #compute distances from hybrid to all others
    d[1, j] = c + j * n - (j - 1) ^ 2
  }
  e = n ^ 2 - 3 * n + 2
  for (i in 2:(n - 1)) {
    #compute distances between non-hybrids
    for (j in (i + 1):n) {
      d[i, j] = e - (i - 2) ^ 2 - (n - j) ^ 2
    }
  }
  return(d + t(d))
}

####################################################
#' Compute fit of circular orders to distance with least squares
#'
#' Compute residual sum of squares (RSS) comparing empirical distance for a blob
#' to an expected one for a cycle with each given order/designated hybrid. This is used
#' in NANAUQ+ commands for resolving multifurcations in a tree of blobs to a cycle
#'
#'@param D an empirical distance table
#'@param E an expected distance table, to be reordered
#'@param orders a vector indicating an order, or matrix whose rows give orders, to fit
#'
#'@return vector of RSSs, one for each order
#'
#'@seealso \code{\link{resolveCycle}}, \code{\link{resolveLevel1}}
#'
#'@export
fitCycleOrders = function(D, E, orders) {
  if (is.vector(orders)) {
    # if only 1 order...
    orders = matrix(orders, nrow = 1) # make it a row of a matrix
  }
  m = dim(orders)[1]# number of orders
  res2 = rep(0, m) # storage for residuals


  for (i in 1:m) {
    res2[i] = sum((D - E[orders[i, ], orders[i, ]]) ^ 2)
  }
  return(res2)
}

#################################################
#' Generate permutations
#'
#' Generate all permutations of 1 to n, as rows of a matrix
#'
#'@param n size of permutations
#'
#'@return an n!xn matrix whose rows give permutations
#'
#'@examples
#' allPerms(4)
#'
#'@export
allPerms = function(n) {
  if (n == 2) {
    return(matrix(c(1, 2, 2, 1), 2, 2))
  }
  else {
    a = allPerms(n - 1)
    asize = factorial(n - 1)
    b = matrix(0, factorial(n), n)
    ncol = matrix(n, asize, 1)
    b[1:asize, 1:(n - 1)] = a # first block
    b[1:asize, n] = ncol
    for (i in 1:(n - 1)) {
      b[(1 + i * asize):((i + 1) * asize), 1:(n - i - 1)] = a[, 1:(n - i - 1)]
      b[(1 + i * asize):((i + 1) * asize), n - i] = ncol
      b[(1 + i * asize):((i + 1) * asize), (n - i + 1):n] = a[, (n - i):(n - 1)]
    }
  }
  return(b)
}

#################################################
#' Generate all circular orders with designated hybrid
#'
#' Generate a matrix whose rows give all circular orders with a designated hybrid.
#' The order is encoded as a vector with entries from 1 to n, where the position
#' corresponds to a node/group of taxa.  The location in the vector of the 1
#' indicates the hybrid, the positions of 2, n its neighbors, etc.
#'
#'@details To avoid duplication of circular orders, the entry 2 in each vector
#'always occurs before the entry n.
#'
#'Since in using first-order quartet-based methods to infer 4-cycles the hybrid
#'node is not identifiable, for n=4 only 3 orders are given, with 1 as hybrid for each
#'
#'@param n size of order, with n>3
#'
#'@return an (n!/2)xn  (or 3xn if n=4) matrix whose rows give all circular orders.
#'
#'@examples
#' circHybOrders(4)
#' circHybOrders(5)
#'
#'@export
circHybOrders = function(n) {
  if (n < 4) {
    stop("Circular orders produced only for cycle size of at least 4.")
  }
  if (n == 4) {
    a = matrix(c(1, 1, 1, 2, 2, 3, 3, 4, 2, 4, 3, 4), 3, 4)
  } else {
    a = allPerms(n) # generate twice what we need
    for (i in (dim(a)[1]:2)) {
      loc1 = which(a[i, ] == 2)
      locn = which(a[i, ] == n)
      if (loc1 > locn)  {
        a = a[-i, ] #delete second representation of circular order w hybrid
      }
    }
  }
  return(a)
}

###########################################################
#'Expected modNANUQ cycle distance
#'
#'Compute expected modNANUQ distance for a sunlet network with 4 or more taxa.
#'This is used in a hueristic method for resolving multifurcation in a tree
#'of blobs to a cycle by NANUQ+ commands.
#'
#'@param n number of edges in cycle (at least 4)
#'
#'@return an nxn distance matrix with rows/columns ordered from hybrid following circular order
#'
#'@examples
#' expmodNANUQCycleDist(5)
#'
#'#'@seealso \code{\link{resolveCycle}} \code{\link{ordersHeuristicmodNANUQ}}
#'
#'@export
expmodNANUQCycleDist = function(n) {
  if (n < 4)
    stop('modNANUQ cycle distance requires at least 4 taxa.')
  D = matrix(0, n, n)

  # get row 1 first, d(x_1,x_k)
  for (k in 2:n) {
    D[1, k] = 2 * ((k - 2) * (n - k) + 1 / 2 * (choose(n - k, 2) +
                                                  choose(k - 2, 2))) + 2 * n - 4
  }
  # get other rows next
  for (k in 2:(n - 1)) {
    for (p in (k + 1):n) {
      D[k, p] =  2 * ((p - k - 1) * (n - p + k - 2) + choose(p - k - 1, 2) + (p -
                                                                                k - 1) + 1 / 2 * (n - p + k - 2) + (n - p) * (k - 2) + 1 / 2 *
                        (choose(n - p, 2) + choose(k - 2, 2))
      ) + 2 * n - 4
    }
  }
  return(D + t(D))
}

##############################################
#' Choose cycle orders heuristically from  empirical modNANUQ distance
#'
#' Find candidates for best hybrid node and circular order fitting the modNANUQ
#' distance.
#'
#' Candidadte orders are obtained  by first picking the hybrid node (from the
#' minimum column sum of the distance matrix), then ordering nodes by distance
#' from the hybrid, and for each consecutive pair picking  nodes in the cycle
#' closest to the previous node. This constructs one or more orders since ties
#' may occur. For more details, see \insertCite{ABRW24;textual}{MSCquartets}.
#'
#' This function is used by NANUQ+ commands
#' to resolve multifurcations in a tree of blobs of high degree.
#'
#'@references
#' \insertRef{ABMR24}{MSCquartets}
#'
#'@param M an empirical modNANUQ distance table
#'@param delta cutoff for relative difference in distances for determining near ties for "best" orders
#'
#'@return a list of circular orders
#'
#' @seealso \code{\link{expmodNANUQCycleDist}} \code{\link{resolveCycle}} \code{\link{resolveLevel1}}
#'
#'@export
ordersHeuristicmodNANUQ = function(M, delta = 10 ^ -6) {
  m = dim(M)[1] #size of matrix = size of cycle
  n = colSums(M)

  minn = min(n)
  hybrids = which((n - minn) <= delta * minn) # find all possible hybrid node
  allords = matrix(0, length(hybrids), m) #allocate space for orders
  allords[, 1] = hybrids # store hybrid groups first

  for (i in 2:m)
  {
    #fill out rest of entries in orders
    numords = dim(allords)[1] # get number of orders being filled
    for (j in 1:numords)
      # for each order
    {
      v = allords[j, ] # copy order
      lastentry = v[i - 1] # get last filled entry in current order
      currentrow = M[lastentry, ] #get row of distances to see what's closest to last filled entry
      currentrow[lastentry] = Inf #make 0 entry large to ignore for min
      vals = unique(sort(currentrow)) # get sorted values
      closest = NULL
      while (length(closest) == 0) {
        #look for closest entries not yet used
        closest = which((currentrow - vals[1]) <= delta * vals[1]) # find closest groups
        closest = closest[which(!(closest %in% v))] # removed any already used
        vals = vals[2:length(vals)] #and remove this value in case we need to look again
      }
      allords[j, i] = closest[1] #fill in one
      if (length(closest) > 1) {
        # if more than 1...
        v = allords[j, ] # copy order
        for (k in 2:length(closest))
        {
          v[i] = closest[k] # add new entry
          allords = rbind(allords, v) # and tack onto order array
        }
      }
    }
  }
  rownames(allords) = NULL

  #put circular orders in canonical form,  and remove duplicates
  if (m == 4) {
    # if only 4 groups, make group 1 hybrid for canonical order
    for (k in 1:dim(allords)[1]) {
      while (allords[k, 1] > 1) {
        allords[k, ] = c(allords[k, 2:4], allords[k, 1])
      }
    }
  }
  for (j in 1:dim(allords)[1]) {
    #make sure 2nd entry smaller than last
    if (allords[j, 2] > allords[j, m]) {
      allords[j, ] = allords[j, c(1, m:2)]
    }
  }
  allords = unique(allords) #remove duplicates

  #convert orders to format where kth entry give position of kth group
  for (i in 1:dim(allords)[1])
  {
    v = allords[i, ]
    for (j in 1:m)
      allords[i, j] = which(v == j)
  }

  return(allords)
}
