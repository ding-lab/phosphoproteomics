##### iCluster.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# cross level clustering between CNV, RNA, and Proteome data

#setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/iCluster")
setwd("/gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer/pan3can_analysis/iCluster")
#source("/Users/khuang/bin/LIB_exp.R")
source("~/bin/LIB_exp.R")

system("mkdir clusterRdata")
#baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"
baseD = "/gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer"

# libraries
library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)

# resources
cancer_genes = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/cgenes_and_druggable.list', header=FALSE, stringsAsFactors = F)
cgenes = as.vector(t(cancer_genes))

# functions
unfactorize = function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

# BRCA data: load
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep=""))
BRCA_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_CNV_formatted_normalized.txt",sep=""))
BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))
BRCA_PRO = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
#BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
# sync samples
s = intersect(colnames(BRCA_mut),colnames(BRCA_CNV))
s = intersect(s,colnames(BRCA_RNA))
s = intersect(s,colnames(BRCA_PRO))

BRCA_mut.s = BRCA_mut[,s]
BRCA_CNV.s = BRCA_CNV[,s]
BRCA_RNA.s = BRCA_RNA[,s]
BRCA_PRO.s = BRCA_PRO[,s]

# feature selection
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
dim(BRCA_mut.s.n.s)

BRCA_CNV.m = as.matrix(BRCA_CNV.s)
BRCA_CNV.na10 = BRCA_CNV.m[rowSums(is.na(BRCA_CNV.m)) <= 5,]
BRCA_CNV.na10.sd0.5 = BRCA_CNV.na10[rowSds(BRCA_CNV.na10, na.rm=TRUE)>0.5,] 

BRCA_RNA.m = as.matrix(BRCA_RNA.s)
BRCA_RNA.na10 = BRCA_RNA.m[rowSums(is.na(BRCA_RNA.m)) <= 5,]
BRCA_RNA.na10.sd2 = BRCA_RNA.na10[rowSds(BRCA_RNA.na10, na.rm=TRUE)>2,] 

BRCA_PRO.m = as.matrix(BRCA_PRO.s)
BRCA_PRO.na10 = BRCA_PRO.m[rowSums(is.na(BRCA_PRO.m)) <= 5,]
BRCA_PRO.na10.sd1 = BRCA_PRO.na10[rowSds(BRCA_PRO.na10, na.rm=TRUE)>1,] 

tBRCA_mut = t(BRCA_mut.s.n.s)
tBRCA_CNV = t(BRCA_CNV.na10.sd0.5)
tBRCA_RNA = t(BRCA_RNA.na10.sd2)
tBRCA_PRO = t(BRCA_PRO.na10.sd1)

tBRCA_mut[is.na(tBRCA_mut)] = 0
tBRCA_CNV[is.na(tBRCA_CNV)] = mean(tBRCA_CNV,na.rm=T)
tBRCA_RNA[is.na(tBRCA_RNA)] = mean(tBRCA_RNA,na.rm=T)
tBRCA_PRO[is.na(tBRCA_PRO)] = mean(tBRCA_PRO,na.rm=T)

sample_aligned = all(rownames(tBRCA_mut)==rownames(tBRCA_CNV)) && all(rownames(tBRCA_mut)==rownames(tBRCA_RNA)) && all(rownames(tBRCA_mut)==rownames(tBRCA_PRO))
if (sample_aligned){
  # looks like it can't tolerate NA...
  fit.single=iClusterPlus(dt1=tBRCA_mut,dt2=tBRCA_CNV,dt3=tBRCA_RNA,dt4=tBRCA_PRO,
                          type=c("binomial","gaussian","gaussian","gaussian"),
                          lambda=c(0.04,0.90,0.90,0.90),K=2,maxiter=10)
  
  bw.col = colorpanel(2,low="white",high="black")
  col.scheme = alist()
  col.scheme[[1]] = bw.col
  col.scheme[[2]] = bluered(256)
  col.scheme[[3]] = bluered(256)
  col.scheme[[4]] = bluered(256)
  fn = paste(pd, "icluster_test.pdf", sep="_")
  pdf(fn, useDingbats=FALSE)
  plotHeatmap(fit=fit.single,datasets=list(tBRCA_mut,tBRCA_CNV,tBRCA_RNA,tBRCA_PRO),
              type=c("binomial","gaussian","gaussian","gaussian"), col.scheme = col.scheme,
              row.order=c(F,T,T,T),chr=chr,plot.chr=c(F,F,F,F),sparse=c(T,T,T,T),cap=c(F,F,F,F))
  dev.off()
  
}

##### model tuning #####
set.seed(123)
date()
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=12,dt1=tBRCA_mut,dt2=tBRCA_CNV,dt3=tBRCA_RNA,dt4=tBRCA_PRO,
                             type=c("binomial","gaussian","gaussian","gaussian"),
                             scale.lambda=c(0.04,0.90,0.90,0.90),n.lambda=307,maxiter=20)
  save(cv.fit, file=paste("clusterRdata/cv.fit.k",k,".Rdata",sep=""))
}
date()

# ##### example data #####
# data(gbm)
# 
# # select mutation features
# dim(gbm.mut)
# mut.rate=apply(gbm.mut,2,mean)
# gbm.mut2 = gbm.mut[,which(mut.rate>0.02)]
# 
# gbm.cn=gbm.cn[order(rownames(gbm.cn)),]
# # check if all the samples are in the same order for the three data sets
# all(rownames(gbm.mut2)==rownames(gbm.exp))
# 
# fit.single=iClusterPlus(dt1=gbm.mut2,dt2=gbm.exp,
#                           type=c("binomial","gaussian"),
#                           lambda=c(0.04,0.90),K=2,maxiter=10)
# 
# # quick plot to check
# bw.col = colorpanel(2,low="white",high="black")
# col.scheme = alist()
# col.scheme[[1]] = bw.col
# col.scheme[[2]] = bluered(256)
# plotHeatmap(fit=fit.single,datasets=list(gbm.mut2,gbm.exp),
#             type=c("binomial","gaussian"), col.scheme = col.scheme,
#             row.order=c(F,T),chr=chr,plot.chr=c(F,F),sparse=c(T,T),cap=c(F,F))
# 
# # model tuning
# set.seed(123)
# date()
# for(k in 1:5){
#   cv.fit = tune.iClusterPlus(cpus=12,dt1=gbm.mut2,dt2=gbm.cn,dt3=gbm.exp,
#                                type=c("binomial","gaussian","gaussian"),K=k,n.lambda=185,
#                                scale.lambda=c(1,1,1),maxiter=20)
#   save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
#   }
# date()
# 
# # model selection
# output=alist()
# files=grep("cv.fit",dir())
# for(i in 1:length(files)){
#   load(dir()[files[i]])
#   output[[i]]=cv.fit
#   }
# nLambda = nrow(output[[1]]$lambda)
# nK = length(output)
# BIC = getBIC(output)
# devR = getDevR(output)
# 
# minBICid = apply(BIC,2,which.min)
# devRatMinBIC = rep(NA,nK)
# for(i in 1:nK){
#   devRatMinBIC[i] = devR[minBICid[i],i]
#   }
# 
# plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
#        ylab="%Explained Variation")
# 
# chr=unlist(strsplit(colnames(gbm.cn),"\\."))
# chr=chr[seq(1,length(chr),by=2)]
# chr=gsub("chr","",chr)
# chr=as.numeric(chr)
# #truncate the values for a better image plot
#   cn.image=gbm.cn
# cn.image[cn.image>1.5]=1.5
# cn.image[cn.image< -1.5]= -1.5
# exp.image=gbm.exp
# exp.image[exp.image>2.5]=2.5
# exp.image[exp.image< -2.5]= -2.5
# 
# # feature selection
# features = alist()
# features[[1]] = colnames(gbm.mut2)
# features[[2]] = colnames(gbm.cn)
# features[[3]] = colnames(gbm.exp)
# sigfeatures=alist()
# for(i in 1:3){
#   rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
#   upper=quantile(rowsum,prob=0.75)
#   sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
#   }
# names(sigfeatures)=c("mutation","copy number","expression")
# #print a few examples of selected features
#   head(sigfeatures[[1]])
# 
# # plot
# bw.col = colorpanel(2,low="white",high="black")
# col.scheme = alist()
# col.scheme[[1]] = bw.col
# col.scheme[[2]] = bluered(256)
# col.scheme[[3]] = bluered(256)
# plotHeatmap(fit=best.fit,datasets=list(gbm.mut2,cn.image,exp.image),
#               type=c("binomial","gaussian","gaussian"), col.scheme = col.scheme,
#               row.order=c(F,F,T),chr=chr,plot.chr=c(F,T,F),sparse=c(T,F,T),cap=c(F,T,F))
