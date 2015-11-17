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

system("mkdir clusterRdata")

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
if (FALSE){ # don't load directly anymore
  # BRCA data: load
  BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep=""))
  BRCA_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_CNV_formatted_normalized.txt",sep=""))
  BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))
  BRCA_PRO = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
  #BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
  BRCA_clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))
  # sync samples
  s = intersect(colnames(BRCA_mut),colnames(BRCA_CNV))
  s = intersect(s,colnames(BRCA_RNA))
  s = intersect(s,colnames(BRCA_PRO))
  s = intersect(s,colnames(BRCA_clin))
  
  BRCA_mut.s = BRCA_mut[,s]
  BRCA_CNV.s = BRCA_CNV[,s]
  BRCA_RNA.s = BRCA_RNA[,s]
  BRCA_PRO.s = BRCA_PRO[,s]
  BRCA_clin.s = BRCA_clin[,s]
  
  #preProcess = function(mut =, cnv = , rna =, pro = )
  # pre-clustering feature selection
  BRCA_mut.s.ch = unfactorize(BRCA_mut.s[row.names(BRCA_mut.s) %in% cgenes,])
  BRCA_mut.s.n = BRCA_mut.s.ch
  for (i in 1:nrow(BRCA_mut.s.ch)){
    BRCA_mut.s.n[i,][BRCA_mut.s.n[i,] == "wt"] = 0
    BRCA_mut.s.n[i,][BRCA_mut.s.n[i,] == "intronic"] = 0
    BRCA_mut.s.n[i,][BRCA_mut.s.n[i,] == "silent"] = 0
    BRCA_mut.s.n[i,][BRCA_mut.s.n[i,] == "RNA"] = 0
    BRCA_mut.s.n[i,][BRCA_mut.s.n[i,] != 0] = 1
  }
  BRCA_mut.s.n = as.data.frame(sapply(BRCA_mut.s.n, as.numeric))
  row.names(BRCA_mut.s.n) = row.names(BRCA_mut.s.ch)
  BRCA_mut.s.n.s = BRCA_mut.s.n[rowSums(BRCA_mut.s.n)/length(s) >= 0.02,]
  #dim(BRCA_mut.s.n.s)
  
  BRCA_CNV.m = as.matrix(BRCA_CNV.s)
  BRCA_CNV.na10 = BRCA_CNV.m[rowSums(is.na(BRCA_CNV.m)) <= 5,]
  BRCA_CNV.na10.sd0.5 = BRCA_CNV.na10[rowSds(BRCA_CNV.na10, na.rm=TRUE)>0.5,] 
  
  BRCA_RNA.m = as.matrix(BRCA_RNA.s)
  BRCA_RNA.na10 = BRCA_RNA.m[rowSums(is.na(BRCA_RNA.m)) <= 5,]
  BRCA_RNA.na10.sd2 = BRCA_RNA.na10[rowSds(BRCA_RNA.na10, na.rm=TRUE)>2,] 
  
  BRCA_PRO.m = as.matrix(BRCA_PRO.s)
  BRCA_PRO.na10 = BRCA_PRO.m[rowSums(is.na(BRCA_PRO.m)) <= 5,]
  BRCA_PRO.na10.sd1 = BRCA_PRO.na10[rowSds(BRCA_PRO.na10, na.rm=TRUE)>1,] 
  
  BRCA_clin.s.s = BRCA_clin.s[c(1:3),]
  
  tBRCA_mut = t(BRCA_mut.s.n.s)
  tBRCA_CNV = t(BRCA_CNV.na10.sd0.5)
  tBRCA_RNA = t(BRCA_RNA.na10.sd2)
  tBRCA_PRO = t(BRCA_PRO.na10.sd1)
  tBRCA_clin = t(BRCA_clin.s.s)
  
  # iCluster can't tolerate NA
  tBRCA_mut[is.na(tBRCA_mut)] = 0
  tBRCA_CNV[is.na(tBRCA_CNV)] = mean(tBRCA_CNV,na.rm=T)
  tBRCA_RNA[is.na(tBRCA_RNA)] = mean(tBRCA_RNA,na.rm=T)
  tBRCA_PRO[is.na(tBRCA_PRO)] = mean(tBRCA_PRO,na.rm=T)
  
  sample_aligned = all(rownames(tBRCA_mut)==rownames(tBRCA_CNV)) && all(rownames(tBRCA_mut)==rownames(tBRCA_RNA)) && all(rownames(tBRCA_mut)==rownames(tBRCA_PRO))
}

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
