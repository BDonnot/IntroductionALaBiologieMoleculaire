install.packages("cghseg")
source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
install.packages("jointseg",repos="http://R-Forge.R-project.org")
library("cghseg")
library("jointseg")
library("rbenchmark")
# library("Segmentor3IsBack")


set.seed(42)
allSizes = c(10,100,200,300,500,1000)
res = rep(0,length(allSizes))
i = 1
for (size in allSizes){
sizeEach = c(350,250,50,350)*size
sigma = .3
gamma = rep(c(2,1,2,3), sizeEach)

c = gamma + rnorm(length(gamma),0,sigma)
# matC= as.matrix(c,ncol = 1)
maxBreakpoints = 3
# classic = cghseg:::segmeanCO(c, maxBreakpoints+1)
pruned = doDynamicProgramming(c, K=maxBreakpoints)
res[i] = system.time(doDynamicProgramming(c, K=maxBreakpoints))["user.self"]
i = i+1
}
# bestPruned = Segmentor(c,keep=T,model = 2,Kmax = maxBreakpoints+1)
if(F){
  plot(c)
  colorsclassic = rep(c("blue","red","darkgreen","black")
                      ,c(classic$t.est[maxBreakpoints+1,1],diff(classic$t.est[maxBreakpoints+1,]))
  )
  colorPruned = rep(c("blue","red","darkgreen","black")
                    ,c(pruned$dpBkpList[[maxBreakpoints]][1]
                       ,diff(c(pruned$dpBkpList[[maxBreakpoints]],sum(sizeEach)))
                    )
  )
  colorPrunedOpt = rep(c("blue","red","darkgreen","black")
                    ,c(bestPruned@breaks[maxBreakpoints+1,1],diff(bestPruned@breaks[maxBreakpoints+1,]))
  )
  
  plot(c,col = rep(c("blue","red","darkgreen","black"),sizeEach), main ="True value", pch = '.')
  plot(c,col = colorsclassic,main = "Calssic DP")
  plot(c,col = colorPruned,main = "Pruned DP")
  plot(c,col = colorPrunedOpt,main = "Pruned DP opt")
}



classicDP = function() cghseg:::segmeanCO(c, maxBreakpoints+1)
prunedDP = function() jointSeg:::jointseg(matC, K=maxBreakpoints,flavor = "DP")
# prunedDPOpt = function() Segmentor(c,keep=T,model = 2,Kmax = maxBreakpoints+1)
benchmark(
  classicDP(),
#   prunedDP(),
  prunedDPOpt(),
  order = "relative",
  columns = c("test","elapsed","relative"),
  replications = 10)
