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
