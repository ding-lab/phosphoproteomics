##### iClusterPlus_tuneModel.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# cross level clustering between CNV, RNA, and Proteome data
# tune iCluster model for best results at each k

setwd("/gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer/pan3can_analysis/iCluster")
source("~/bin/LIB_exp.R")
baseD = "/gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer/"
cancer_genes = read.table(file='/gscmnt/gc2524/dinglab/Proteomics/projects/reference_files/cgenes_and_druggable.list', header=FALSE, stringsAsFactors = F)

if (FALSE){ # for use on my macpro
  setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/iCluster")
  baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"
  source("/Users/khuang/bin/LIB_exp.R")
  cancer_genes = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/cgenes_and_druggable.list', header=FALSE, stringsAsFactors = F)
}

cgenes = as.vector(t(cancer_genes))

#system("mkdir clusterRdata")

# libraries
library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)

args=commandArgs(TRUE)

k=as.numeric(args[1])
#cancer=args[2]

cat("############################")
cat(date())
cat(paste("Tuning iCluster model in BRCA of cluster number: ",k))

# functions
unfactorize = function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

##### MAIN CODE #####
# load data stored from preprocessData script
load("clusterRdata/2015-11-13_BRCA.Rdata")

##### model tuning #####
if (sample_aligned){
  #tune_icluster = function(){}
  set.seed(123) # the greedy algorithm requires setting seed
  date()
  
  cv.fit = tune.iClusterPlus(cpus=12,dt1=tmut,dt2=tCNV,dt3=tRNA,dt4=tPRO, K=k,
                             type=c("binomial","gaussian","gaussian","gaussian"),
                             scale.lambda=c(0.05,0.20,1.00,0.80),n.lambda=307,maxiter=30)
  save(cv.fit, file=paste("clusterRdata/cv.fit.k_lambda005_020_100_080", "_",k,".Rdata",sep=""))
  
  date()
} else {
  warning("Samples across data files are not aligned!!")
}
