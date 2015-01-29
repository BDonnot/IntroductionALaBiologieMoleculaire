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


addressPit="C:\\Users\\Peter-Jack\\Downloads\\GSE17359_RAW\\GSM433918.cel"
addressBenj=""

library(affy) 
affy.data = ReadAffy(filenames=addressPit) 
mypm <- pm(affy.data) 
mymm <- mm(affy.data) 
myaffyids <- probeNames(affy.data) 