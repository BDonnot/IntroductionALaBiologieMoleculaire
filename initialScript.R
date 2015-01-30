install.packages("cghseg")
install.packages("jointSeg", repos="http://R-Forge.R-project.org")
library("cghseg")
library("jointSeg")

signal <- rnorm(1e3)
maxBreakpoints <- 10

classic = cghseg:::segmeanCO(signal, maxBreakpoints)
pruned = jointSeg:::pruneByDP(as.matrix(signal,nrow = 1), K=maxBreakpoints)
pruned$rse

#TO DO :
##Comparer les deux algos, avec temps de calcul et precision
#des estimations -> donnees artificielles
#-> verifier que la solution est la meme
##faire un truc avec des vraies donnees

#Rapport :
##presenter la problematique
##presenter l'algo classique
##donner l'intuition de l'algo speed-up 'ed
##montrer les resultats
##conclure

install.packages("cghseg")
install.packages("jointseg", repos="http://R-Forge.R-project.org")
library("cghseg")
library("jointseg")


##
##source("http://bioconductor.org/biocLite.R")
##biocLite("DNAcopy")

signal <- rnorm(1e3)
maxBreakpoints <- 10

cghseg:::segmeanCO(signal, maxBreakpoints) 
jointseg:::pruneByDP(signal, K=maxBreakpoints)

##source("http://bioconductor.org/biocLite.R")
##biocLite("affy")


 library("jointseg")

## load known real copy number regions
affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=0.7)

## generate a synthetic CN profile

K <- 10
len <- 1e4
sim <- getCopyNumberDataByResampling(len, K, regData=affyDat)
datS <- sim$profile
signal=datS$c
Res=cghseg:::segmeanCO(signal, K)
sim$bkp
Res$t.est[10,] 
jointseg:::pruneByDP(signal,candCP =sim$bkp,K=K)
RES_RBS=jointSeg(signal,method="RBS", K=K)
RES_RBS$dpBkpList[10]
sim$bkp
Res$t.est[10,] 
par(mfrow=c(4,4))
plotSeg(datS, sim$bkp)
plotSeg(datS,Res$t.est[10,])


affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=0.5)

## generate a synthetic CN profile
K <- 10
len <- 1e4
sim <- getCopyNumberDataByResampling(len, K, regData=affyDat)
datS <- sim$profile

## run binary segmentation (+ dynamic programming) resRBS <-
resRBS <- PSSeg(data=datS, method="RBS", stat=c("c", "d"), K=2*K, profile=TRUE)
resRBS$prof

getTpFp(resRBS$bestBkp, sim$bkp, tol=5)
plotSeg(datS, breakpoints=list(sim$bkp, resRBS$bestBkp))



