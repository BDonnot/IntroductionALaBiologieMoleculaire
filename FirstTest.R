install.packages("cghseg")
install.packages("jointseg", repos="http://R-Forge.R-project.org")
library("cghseg")
library("jointseg")

signal <- rnorm(1e3)
maxBreakpoints <- 10

cghseg:::segmeanCO(signal, maxBreakpoints) jointseg:::pruneByDP(signal, K=maxBreakpoints)
