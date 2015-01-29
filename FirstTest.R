library("cghseg")
library("jointSeg")
library("rbenchmark")
library("Segmentor3IsBack")

set.seed(42)
sizeEach = c(350,250,50,350)*1e4
sigma = .3
gamma = rep(c(2,1,2,3), sizeEach)

c = gamma + rnorm(length(gamma),0,sigma)
# matC= as.matrix(c,nrow = 1)
maxBreakpoints = 3
classic = cghseg:::segmeanCO(c, maxBreakpoints+1)
pruned = jointSeg:::jointseg(matC, K=maxBreakpoints,flavor = "DP")
c(pruned$dpBkpList[[maxBreakpoints]][1]
  ,diff(c(pruned$dpBkpList[[maxBreakpoints]],sum(sizeEach)))) == sizeEach

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
