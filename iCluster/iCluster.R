##### iCluster.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# cross level clustering between CNV, RNA, and Proteome data

#setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/cross_correlation/")
source("~/bin/LIB_exp.R")

baseD = "/gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer/"
setwd(paste(baseD, "pan3can_analysis/iCluster", sep=""))

#requirements

library(iCluster)
library(caTools)
library(gdata)
library(gtools)
library(gplots)
library(lattice)

# example
data(gbm)
str(gbm)
#setting the penalty parameter lambda=0 returns non-sparse fit
fit=iCluster(datasets=gbm, k=3, lambda=list(0.44,0.33,0.28))
plotiCluster(fit=fit, label=rownames(gbm[[1]]))
compute.pod(fit)
data(coord)
chr=coord[,1]
plotHeatmap(fit=fit, data=gbm, feature.order=c(FALSE,TRUE,TRUE),
sparse=c(FALSE,TRUE,TRUE),plot.chr=c(TRUE,FALSE,FALSE), chr=chr)
