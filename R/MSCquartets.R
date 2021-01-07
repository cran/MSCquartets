#' Multispecies Coalescent Model Quartet Package
#'
#' A package for analyzing quartets displayed on gene trees, under the multispecies coalescent (MSC) model.
#'
#' @details
#' This package contains routines to analyze a collection of gene trees through the displayed 
#' quartets on them.
#' 
#' A quartet count concordance factor (qcCF) for a set of 4 taxa is the triple of counts of the three
#' possible resolved quartet trees on those taxa across some set of gene trees. The major routines in this package can:
#' \enumerate{
#' \item Tabulate all qcCFs for a collection of gene trees.
#' \item Perform hypothesis tests of whether one or more qcCFs are consistent with the MSC model on
#' a species tree \insertCite{MAR19}{MSCquartets}.
#' \item Produce simplex plots showing all estimated CFs as well as results of hypothesis tests \insertCite{AMR2020}{MSCquartets}.
#' \item Infer a species tree using the qcCFs via the QDC and WQDC methods \insertCite{Rho19,YR19}{MSCquartets}. 
#' \item Infer a level-1 species network via the NANUQ method \insertCite{ABR19}{MSCquartets}.
#' }
#'As discussed in the cited works, the inference methods for species trees and networks are
#'statistically consistent under the MSC and Network MSC respectively. 
#'
#'This package, and the theory on which it is based, allows gene trees to have 
#'missing taxa (i.e., not all gene trees display all the taxa). It does require 
#'that each subset of 4 taxa is displayed on at least one of the gene trees.
#'
#'
#' Several gene tree data sets, simulated and empirical, are included.
#'   
#' In publications please cite the software announcement \insertCite{RBMA2020}{MSCquartets}, as well as the
#' appropriate paper(s) developing the theory behind the routines you used.
#' 
#' @references
#' 
#' \insertRef{RBMA2020}{MSCquartets}
#' 
#' \insertRef{MAR19}{MSCquartets}
#' 
#' \insertRef{AMR2020}{MSCquartets}
#'    
#' \insertRef{Rho19}{MSCquartets}
#' 
#' \insertRef{YR19}{MSCquartets}
#' 
#' \insertRef{ABR19}{MSCquartets}
#' 
#' @importFrom graphics hist legend mtext par plot points segments text
#' @importFrom stats chisq.test dmultinom pchisq qnorm
#' @importFrom Rdpack reprompt
#'
#' @docType package
#' @name MSCquartets-package
NULL





#' Yeast gene tree dataset
#'
#' A text file dataset for Yeast containing 106 gene trees on 8 taxa  (7 Saccharomyces and 1 Candida outgroup). This is a subset of the data
#' of \insertCite{Rokas03;textual}{MSCquartets}.
#'
#' @docType data
#' 
#' @name dataYeastRokas
#' 
#' @references \insertRef{Rokas03}{MSCquartets}
#' 
#' @details File is accessed as \code{system.file("extdata","dataYeastRokas",package="MSCquartets")}, for example
#' via the \code{ape} command:
#' 
#' \code{ gts=read.tree(file = system.file("extdata","dataYeastRokas",package="MSCquartets"))}
#' 
#' @format A text file with 106 topological Newick gene trees on the taxa: 
#'  Sbay, Scas, Scer, Sklu, Skud, Smik, Spar, and Calb (outgroup).
#'
#' @source \url{https://wiki.rice.edu/confluence/download/attachments/8898533/yeast.trees?version=1&modificationDate=1360603275797&api=v2}
NULL


#' Heliconius gene tree dataset
#'
#' A text file dataset for Heliconius butterflies containing 2909 gene trees on 7 taxa, with 4 individuals sampled for each of 3 of the taxa, for a total of 16
#' leaves per gene tree. This is a subset of the data of \insertCite{Martin2013;textual}{MSCquartets}.
#'
#' @references \insertRef{Martin2013}{MSCquartets}
#' 
#' @docType data
#' 
#' @name dataHeliconiusMartin
#' 
#' @details File is accessed as \code{system.file("extdata","dataHeliconiusMartin",package="MSCquartets")}, for 
#' example via the \code{ape} command:
#' 
#' \code{gts = read.tree(file=system.file("extdata","dataHeliconiusMartin",package="MSCquartets"))}
#' 
#' @format A text file with 2909 metric Newick gene trees each with 16 leaves labelled:\cr
#' chioneus.553, chioneus.560, chioneus.564, chioneus.565, \cr
#' ethilla.67, hecale.273, melpomeneFG.13435, melpomeneFG.9315, \cr 
#' melpomeneFG.9316, melpomeneFG.9317, pardalinus.371, rosina.2071, \cr 
#' rosina.531, rosina.533, rosina.546, sergestus.202
#'
#' @source \doi{10.5061/dryad.dk712}
#' 
NULL 


#' Simulated gene tree dataset from species tree
#'
#' A text file dataset containing 1000 gene trees on 9 taxa simulated under the MSC on a species tree
#'
#' @details This simulated dataset was produced by SimPhy \insertCite{simphy}{MSCquartets}, using the species tree
#' 
#' ((((t5:5000,t6:5000):5000,t4:10000):2500,t7:12500):7500,((t8:3000,t9:3000):5000,\cr
#' ((t1:4000,t2:4000):2500,t3:6500):1500):12000);
#' 
#' with a population size of 10,000 throughout the tree.
#' 
#' @docType data
#' 
#' @name dataGeneTreeSample
#' 
#' @details File is accessed as  \code{system.file("extdata","dataGeneTreeSample",package="MSCquartets")}, for example
#' via the ape command:
#' 
#' \code{gts=read.tree(file = system.file("extdata","dataGeneTreeSample",package="MSCquartets") )}
#'  
#' @format A text file with 1000 metric Newick gene trees on the taxa t1-t9
#' 
#' @references 
#' \insertRef{simphy}{MSCquartets}
#' 
NULL 
