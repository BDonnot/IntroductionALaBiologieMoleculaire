install.packages("cghseg")
install.packages("jointSeg", repos="http://R-Forge.R-project.org")
library("cghseg")
library("jointSeg")

signal <- rnorm(1e3)
maxBreakpoints <- 10

classic = cghseg:::segmeanCO(signal, maxBreakpoints)
pruned = jointSeg:::pruneByDP(as.matrix(signal,nrow = 1), K=maxBreakpoints)
pruned$rse
