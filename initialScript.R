install.packages("cghseg")
install.packages("jointSeg", repos="http://R-Forge.R-project.org")
library("cghseg")
library("jointseg")

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
sim$bkp

datS <- sim$profile
signal=datS$c
Res_1=cghseg:::segmeanCO(signal, K)
sim$bkp
Res_1$t.est[10,] 
Res_2=doDynamicProgramming(signal,K=K)
Res_2$dpseg$bkp[[10]]
sum(Res_2$dpseg$rse)

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

<<results=tex>>
    xtable(t(rbind(c(1:11),Res_2$dpseg$rse)))
@
sim$dk
signal

if (F){
install.packages("reshape2")
install.packages("data.table",repos="http://R-Forge.R-project.org")
}
library(data.table)
library("cghseg")
library("jointseg")

## load known real copy number regions
set.seed(42)
affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=0.7)
K <- 10
len <- 1e4
sim <- getCopyNumberDataByResampling(len, K, regData=affyDat)
segMean = cghseg:::segmeanCO(sim$profile$c, K+1)
jointSeg=doDynamicProgramming(sim$profile$c,K=K)
sepM = segMean$t.est[K+1,]
sepJ = jointSeg$bkp
sepTh = sim$bkp
mySim = data.table(c = sim$profile$c,
                   clustTh = rep(1:11,c(sepTh[1],diff(c(sepTh,len)))),
                   clustJ = rep(1:11,c(sepJ[1],diff(c(sepJ,len)))),
                   clustM = rep(1:11,c(sepM[1],diff(sepM))))

# mySim[,var(c)*.N,by = clustJ]
# mySim[,var(c)*.N,by = clustM]
# mySim[,var(c)*.N,by = clustTh]
sum(mySim[,var(c)*.N,by = clustJ],na.rm = T) #there is sometimes only 1 point in a cluster, hence the na.rm
sum(mySim[,var(c)*.N,by = clustM],na.rm = T)
sum(mySim[,var(c)*.N,by = clustTh],na.rm = T)



length(rep(1:11
           ,c(sim$bkp[1]
              ,diff(c(sim$bkp,len))))
       )


