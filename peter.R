source("http://bioconductor.org/biocLite.R")
biocLite()

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")
biocLite("affy")
biocLite("gcrma")
biocLite("hugene10stv1cdf")
biocLite("hugene10stv1probe")
n
biocLite("hugene10stprobeset.db")
n
biocLite("hugene10sttranscriptcluster.db")
n

#Load the necessary libraries
library(GEOquery)
library(affy)
library(gcrma)
library(hugene10stv1cdf)
library(hugene10stv1probe)
library(hugene10stprobeset.db)
library(hugene10sttranscriptcluster.db)
library(R.utils)



addressPit="C:/Users/Peter-Jack/Desktop/GSE17359_RAW"
setwd(addressPit)
untar("GSE17359_RAW.tar", exdir="data")
cels = list.files( pattern = "txt")
table=read.table(cels[1],header=T)
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")
setwd("C:/Users/Peter-Jack/Desktop/GSE17359_RAW/data")
cels = list.files("C:/Users/Peter-Jack/Desktop/GSE17359_RAW/data", pattern = "CEL")
cels[1]
ReadAffy(filenames=cels[1])
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="hugene10stv1") #From bioconductor

data.rma.norm=rma(raw.data)

#Get the important stuff out of the data - the expression estimates for each array
rma=exprs(data.rma.norm)

#Format values to 5 decimal places
rma=format(rma, digits=5)

#Map probe sets to gene symbols or other annotations
#To see all available mappings for this platform
ls("package:hugene10stprobeset.db") #Annotations at the exon probeset level
ls("package:hugene10sttranscriptcluster.db") #Annotations at the transcript-cluster level (more gene-centric view)

#Extract probe ids, entrez symbols, and entrez ids
probes=row.names(rma)
Symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hugene10sttranscriptclusterENTREZID, ifnotfound=NA))

#Combine gene annotations with raw data
rma=cbind(probes,Symbols,Entrez_IDs,rma)
attach(rma)
rma$GSM433916.CEL
head(rma[,4])

library("cghseg")
library("jointseg")

signal <- rnorm(1e3)
maxBreakpoints <- 10
install.packages("matrixStats")

cghseg:::segmeanCO(y, maxBreakpoints)
jointseg:::pruneByDP(signal, K=maxBreakpoints)