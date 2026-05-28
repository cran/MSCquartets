## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")
 
knitr::opts_chunk$set(fig.align = "center",
                      fig.show = "hold",
                      out.width = "55%",
                      fig.width = 7,  
                      fig.height = 6)



## ----setup--------------------------------------------------------------------
library(MSCquartets)

## ----eval=FALSE---------------------------------------------------------------
# # Read gene trees and compute quartet counts
# gts <- read.tree(file = "genetreefile")
# tableLeopardusLescroart <- quartetTable(gts)

## -----------------------------------------------------------------------------
AstralLeopard_tree <-   read.tree(text = ASTRALtreeLeopardusLescroart)

## ----eval=FALSE---------------------------------------------------------------
# gts <- read.tree(file = "genetreefile")
# alternative_tree <- quartetDistTree(gtrees)

## -----------------------------------------------------------------------------
out <- ECToBlob(tableLeopardusLescroart, AstralLeopard_tree, alpha = 0.05)

pTable <- out$pTable

## -----------------------------------------------------------------------------
out_1 <- ECToBlob(pTable, AstralLeopard_tree, alpha=0.05, qType="mul",  testCorrection="Bon",    plot=0)
out_2 <- ECToBlob(pTable, AstralLeopard_tree, alpha=0.05, qType="quad", testCorrection="Cauchy", plot=0)
out_3 <- ECToBlob(pTable, AstralLeopard_tree, alpha=0.05, qType="bi",   testCorrection="BBC",    plot=0)
out_4 <- ECToBlob(pTable, AstralLeopard_tree, alpha=0.05, qType="mul",  testCorrection="CBC",    plot=0)

## -----------------------------------------------------------------------------
early_tree <- out_2$treeList[[out_2$indexEarly]]$tree
plot(read.tree(text = early_tree))

## ----eval=FALSE---------------------------------------------------------------
# late_tree <- out_2$treeList[[out_2$indexLate]]$tree
# plot(read.tree(text = late_tree))

