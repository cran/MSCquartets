#' Produce table of counts of quartets displayed on trees
#'
#' Compiles table of quartet count concordance factors (qcCFs) for topological quartets displayed on a 
#' collection of trees.
#' 
#'
#' @details
#' The names in \code{taxonnames} may be any subset of those on the trees.
#' Branch lengths of non-negative size less than or equal to \code{epsilon}
#' are treated as zero, giving polytomies.
#'  
#'  In the returned table, columns are labeled by taxon names and quartet names ("12|34", etc.).
#'  1s and 0s in taxon columns indicate the taxa in a quartet. Quartet 12|34
#'  means the first and second indicated taxa form a cherry, 13|24 means the first and third form a cherry, 14|23 means
#'  the first and fourth form a cherry, and 1234 means the quartet is unresolved.
#' 
#' An error occurs if any branch length is negative. 
#' Warnings are given if some of \code{taxonnames} are missing on some trees, or
#' if some 4-taxon set is not on any tree.
#'
#' If \code{random}>0, then for efficiency \code{random} should be much smaller then
#' the number of possible 4 taxon subsets.
#'
#' @param trees multiphylo object containing un/rooted metric/topological trees
#' @param taxonnames vector of \code{n} names of taxa of interest; if \code{NULL} then taken from taxa on \code{trees[[1]]}
#' @param epsilon minimum for branch lengths to be treated as non-zero
#' @param random number of random subsets of 4 taxa to consider; if 0, use all \code{n} choose 4 subsets
#' @return
#'     an (\code{n} choose 4)x(\code{n}+4) matrix (or (\code{random})x(\code{n}+4) matrix) encoding
#'     4 taxon subsets of \code{taxonnames} and counts of each of the
#'     quartets 12|34, 13|24, 14|23, 1234 across the trees
#'     
#'     
#' @seealso \code{\link{quartetTableResolved}}, \code{\link{quartetTableDominant}}, \code{\link{taxonNames}}
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' DQT=quartetTableDominant(RQT)
#'
#' @importFrom ape cophenetic.phylo compute.brlen
#'
#' @export
quartetTable = function(trees,
                        taxonnames=NULL,
                        epsilon = 0,
                        random = 0) {
  if (random < 0)
    stop("Parameter 'random' must be non-negative.")
  random = round(random)
  
  if (is.null(taxonnames)) {
    taxonnames=trees[[1]]$tip.label #get taxa names from 1st tree
    message("Using taxa that appear on first gene tree.")
  }
  
  taxonnames = sort(taxonnames) #put taxa in fixed order
  nt = length(trees) # number of gene trees
  N = length(taxonnames) # number of taxa
  if (random > 0)
    M = random # number of 4-taxon sets to consider
  else
    M = choose(N, 4)
  Q = matrix(0, M, N + 4) #allocate space for table
  qnames = c("12|34", "13|24", "14|23", "1234")
  colnames(Q) = c(taxonnames, qnames) #create column names
  warnMissing = 0 # flag for some taxa not on all gene trees
  
  # encode 4-taxon sets in matrix
  if (random == 0) {
    #use all subsets of 4
    m = 0
    for (i in 1:(N - 3)) {
      # for each 4-taxon set
      for (j in (i + 1):(N - 2)) {
        for (k in (j + 1):(N - 1)) {
          for (l in (k + 1):N) {
            m = m + 1
            Q[m, c(i, j, k, l)] = 1 #encode set
          }
        }
      }
    }
  }
  else {
    #use random sets of 4 taxa
    i = 1
    while (i <= random)  {
      q = sample(N, size = 4, replace = FALSE)# 4 random taxa
      row = integer(N) # encode as 0/1 vector
      row[q] = 1
      
      j = 1 #check earlier rows for repeated quartet
      while (j < i) {
        if (identical(Q[j, 1:N], row))
          j = i + 1 #set flag for match
        else
          j = j + 1 # go to next row
      }
      if (j == i) {
        # accept this choice of 4 taxa
        Q[i, 1:N] = row
        i = i + 1
      }
    }
  }
  
  message(
      'Counting occurrences of displayed quartets for ',
      M,' four-taxon subsets of ',
      N,' taxa across ',
      nt,' gene trees.'
    )

  
  # timing for large gt datasets
  start_time <- Sys.time()
  
  
  for (numt in 1:nt) {
    # for each gene tree
    t = trees[[numt]] #get next tree
    if (is.null(t$edge.length)) {
      t = compute.brlen(t, 1) #if no branch lengths then make them all 1
    } else {
      if (sum(t$edge.length < 0) > 0)
        stop("Error: Negative branch length in tree") #if any negative lengths abort
    }
    
    zeros = which(t$edge.length <= epsilon) # find locations of any zero branch lengths
    t$edge.length[] = 1 # set all lengths to 1, to avoid numerical error in 4-point below
    t$edge.length[zeros] = 0  # reset zeros to 0
    
    d = cophenetic.phylo(t) # get distance matrix for the tree
    for (m in 1:M) {
      # for each 4-taxon set
      tax = as.character(taxonnames[which(Q[m, 1:N] == 1)]) # names of 4 taxa (or fewer)
      if (all(tax %in% colnames(d))) {
        # if all 4 taxa present on this tree
        a = d[tax[1], tax[2]] + d[tax[3], tax[4]] # use 4-point condition
        b = d[tax[1], tax[3]] + d[tax[2], tax[4]]  #    to determine quartet
        c = d[tax[1], tax[4]] + d[tax[2], tax[3]]  #
        z = sort(c(a, b, c))
        if (z[1] == z[2]) {
          # tie means a polytomy
          Q[m, "1234"] = Q[m, "1234"] + 1
        } else {
          if (z[1] == a) {
            #quartet is ab|cd
            Q[m, "12|34"] = Q[m, "12|34"] + 1
          } else {
            if (z[1] == b) {
              #quartet is ac|bd
              Q[m, "13|24"] = Q[m, "13|24"] + 1
            } else {
              Q[m, "14|23"] = Q[m, "14|23"] + 1
            } #quartet is ad|bc
          }
        }
      } else
        warnMissing = 1 # some taxa not on this tree
    }
  }
  
  if (warnMissing == 1)
    warning("Some taxa missing from some trees.")
  if ((N > 4) && (sum(rowSums(Q[, qnames, drop = FALSE]) == 0) > 0))
    warning("Some 4-taxon set not present on any tree.")
  if (sum(Q[, "1234"]) > 0)
    warning("Some quartets unresolved.")
  
  current_time = Sys.time()
  
  elapsedTime=as.numeric(difftime(current_time,start_time,units="mins")) #compute timing
  if (elapsedTime>15) {
    message("Time to process quartets on gene trees was ",elapsedTime," minutes.")
  }
  
  return(Q)
}

#################################################

#' Modify quartet table to show only resolved quartets
#'
#' Converts table of all quartet counts, including unresolved ones, 
#' by either dropping unresolved ones, or distributing them uniformly 
#' among the three resolved counts.
#'
#' @param qt table, as produced by \code{quartetTable} for \code{n} taxa, with \code{n+4} columns
#' @param omit \code{TRUE} deletes unresolved quartets column;
#'            \code{FALSE} deletes the column but redistributes unresolved counts as (1/3,1/3,1/3) to resolved counts
#' @return
#'      a table of quartet counts similar to \code{qt}, but with columns showing only resolved quartet counts
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' QT[1:6,]
#' RQT=quartetTableResolved(QT)
#' RQT[1:6,]
#'
#' @seealso \code{\link{quartetTable}}, \code{\link{quartetTableDominant}}
#' @export
quartetTableResolved = function(qt, omit = FALSE) {
  if ("1234" %in% colnames(qt)) { # if unresolved column present
    RT = qt[,-which(colnames(qt) == "1234"), drop = FALSE]  #copy table without unresolved column
    if (omit == FALSE) {
      qs = c("12|34", "13|24", "14|23") #columns to be adjusted
      RT[, qs] = RT[, qs, drop = FALSE] + qt[, "1234", drop = FALSE] %*% matrix(1 / 3, 1, 3) #redistribute counts
    }
  } else { #no unresolved column present
    RT = qt
  }
  return(RT)
}

#################################################


#' Produce table of dominant quartets, with estimates of internal edge lengths
#'
#' Converts table of counts of resolved quartets on \code{n} taxa to show only dominant one, with
#' maximum likelihood estimate of internal edge weight under the MSC.
#'
#' @details 
#' If \code{bigweights="finite"}, when for a set of 4 taxa the quartet counts are (m,0,0) then
#' the edge weight is computed as if the relative frequency of the dominant topology were m/(m+1).
#' 
#' @param rqt a table, as produced by \code{quartetTableResolved} of size (n choose 4)x(n+3);
#' @param bigweights \code{"infinite"} or \code{"finite"}, to indicate whether the weight (internal edge length) 
#' of a quartet for which only one
#' topology appears is given as \code{Inf} or a finite, but large, numerical value
#' 
#' @return
#'   an (n choose 4)x(n+1) array with dominant quartet topology encoded by 1,1,-1,-1 in
#' taxon columns, with signs indicating cherries; the (n+1)th column \code{"weight"} contains the maximum likelihood estimates,
#' under MSC on a 4-taxon tree, of the quartets' central edge lengths, in coalescent units
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#' QT=quartetTable(gtrees,tnames[1:6])
#' RQT=quartetTableResolved(QT)
#' RQT[1:6,]
#' DQT=quartetTableDominant(RQT)
#' DQT[1:6,]
#'
#' @seealso \code{\link{quartetTable}}, \code{\link{quartetTableResolved}}
#' @importFrom stats runif
#' @export
quartetTableDominant = function(rqt, bigweights="infinite") {
  M = dim(rqt)[1] # number of 4-taxon sets
  N = dim(rqt)[2] - 3 # number of taxa
  taxonnames = colnames(rqt)[1:N] # names of taxaR
  
  
  DQT = rqt[, 1:(N + 1)] # allocate space for most frequent quartet topologies and weight
  if (is.vector(DQT)){ # if only 1 row,
    dim(DQT)=c(1,N+1) # make it a matrix,
    colnames(DQT)=colnames(rqt)[1:(N+1)] #with column names
    }
  for (m in 1:M) {
    # for each 4-taxon set
    taxa = which(rqt[m, 1:N] == 1) # taxa numbers
    z = rqt[m, (N + 1):(N + 3)] # frequencies of quartets
    tz = sum(z) # total number of quartets
    if (tz > 0) {
      # if some quartets were found
      mq = which.max(z + runif(3)) # which is biggest? (random tie-breaking)
      if (mq == 1) {
        DQT[m, c(taxa[3], taxa[4])] = -1
      } else {
        if (mq == 2) {
          DQT[m, c(taxa[2], taxa[4])] = -1
        } else {
          DQT[m, c(taxa[2], taxa[3])] = -1
        }
      }
      f = z[mq] / tz #relative frequency of most frequent quartet
      if (f<1) { # if some ILS
        DQT[m, N + 1] = -log((1.5 * (1 - f))) # save weight
      } else { # if no ILS
      if (bigweights =="infinite"){
        DQT[m, N + 1] = Inf
      } else { #despite no ILS give finite weight
      f = tz/(tz+1) # pretend f is a bit less than 1
      DQT[m, N + 1] = -log((1.5 * (1 - f))) # save weight
    }
  }
  
    }
  }
  colnames(DQT) = c(taxonnames, "weight") #add column names
  return(DQT)
}

######################################################################################

#' Get all taxon names from a collection of trees
#'
#' Create a vector of all taxon names appearing on a collection of trees, with no repeats.
#'
#' @param trees a multiPhylo object containing a collection of trees
#'
#' @return a vector of unique names of taxa appearing on the trees
#'
#' @examples
#' gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#' tnames=taxonNames(gtrees)
#'
#' @export
taxonNames = function(trees) {
  taxa = c()
  for (i in 1:length(trees)) {
    taxa = union(taxa, trees[[1]]$tip.label)
  }
  return(taxa)
}

###################################################################################
#' Print a quartet table with nice formatting
#' 
#' Print a quartet table with the taxa in each quartet shown by name.
#'
#' @param qt a table such as returned by \code{quartetTable}, \code{quartetTableResolved}, 
#' or \code{quartetTableDominant}, possibly with extra columns added by other functions
#'
#' @return NULL
#'
#' @examples 
#'  gtrees=read.tree(file=system.file("extdata","dataGeneTreeSample",package="MSCquartets"))
#'  tnames=taxonNames(gtrees)
#'  QT=quartetTable(gtrees,tnames[1:6])
#'  QT[1:6,]
#'  quartetTablePrint(QT[1:6,])  
#'  RQT=quartetTableResolved(QT)
#'  RQT[1:6,]
#'  quartetTablePrint(RQT[1:6,])
#'  pTable=quartetTreeTestInd(RQT,"T3")
#'  pTable[1:6,] 
#'  quartetTablePrint(pTable[1:6,])
#'  DQT=quartetTableDominant(RQT)
#'  DQT[1:6,]
#'  quartetTablePrint(DQT[1:6,])
#'  
#'@export
quartetTablePrint = function(qt) {
  # get number of rows and columns
  nrows = dim(qt)[1]
  ncols = dim(qt)[2]
  # best guess as to what columns to display
  colsToDisplay=rep(0,nrows)
  for (i in 1:nrows){
    colsToDisplay[i]=which(qt[i,]!=0)[5]
  }
  minColNumber=min(colsToDisplay)
  colsToDisplay=minColNumber:ncols
  
  # display depends on table entries
  
  if (-1 %in% qt[1,1:(minColNumber-1)]) {
    # formatting for dominant quartet table
    prettyDisplay=matrix(0,nrows,1+length(colsToDisplay))
    for (i in 1:nrows){
      cherry1 = which(qt[i,1:(minColNumber-1)] == 1)
      cherry2 = which(qt[i,1:(minColNumber-1)] == -1)
      qtree = paste0(colnames(qt)[cherry1[1]]," ",colnames(qt)[cherry1[2]]," | ",colnames(qt)[cherry2[1]]," ",colnames(qt)[cherry2[2]]," ")
      prettyDisplay[i,]=c(qtree,qt[i,colsToDisplay])
    }
    colnames(prettyDisplay) <- c(' ',colnames(qt)[colsToDisplay])
  }
  else {
    # formatting for qt and rQT
    prettyDisplay=matrix(0,nrows,4+length(colsToDisplay))
    for (i in 1:nrows){
      prettyDisplay[i,]=c(colnames(qt)[which(qt[i, 1:(minColNumber-1)]!=0)],qt[i,colsToDisplay])
    }
    colnames(prettyDisplay) <- c(rep(' ',4),colnames(qt)[colsToDisplay])
  }
  rownames(prettyDisplay) <-rep('',nrows)
  print(as.table(prettyDisplay))
}
