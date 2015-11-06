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

args=commandArgs(TRUE)

k=args[1]
cancer=args[2]

# resources
cancer_genes = read.table(file='/gscmnt/gc2524/dinglab/Proteomics/projects/reference_files/cgenes_and_druggable.list', header=FALSE, stringsAsFactors = F)
cgenes = as.vector(t(cancer_genes))

# functions
unfactorize = function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

# load cancer type data
mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",cacner,"/",cancer,"_SOMATIC_formatted.txt",sep=""))
CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",cacner,"/",cancer,"_CNV_formatted_normalized.txt",sep=""))
RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",cacner,"/",cancer,"_mRNA_formatted_normalized.txt",sep=""))
PRO = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",cacner,"/",cancer,"_PRO_formatted_normalized.txt",sep=""))
#Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",cacner,"/",cancer,"_PHO_formatted_normalized.txt",sep=""))
# sync samples
s = intersect(colnames(mut),colnames(CNV))
s = intersect(s,colnames(RNA))
s = intersect(s,colnames(PRO))

mut.s = mut[,s]
CNV.s = CNV[,s]
RNA.s = RNA[,s]
PRO.s = PRO[,s]

#preProcess = function(mut =, cnv = , rna =, pro = )
# pre-clustering feature selection
mut.s.ch = unfactorize(mut.s[row.names(mut.s) %in% cgenes,])
mut.s.n = mut.s.ch
for (i in 1:nrow(mut.s.ch)){
  mut.s.n[i,][mut.s.n[i,] == "wt"] = 0
  mut.s.n[i,][mut.s.n[i,] == "intronic"] = 0
  mut.s.n[i,][mut.s.n[i,] == "silent"] = 0
  mut.s.n[i,][mut.s.n[i,] == "RNA"] = 0
  mut.s.n[i,][mut.s.n[i,] != 0] = 1
}
mut.s.n = as.data.frame(sapply(mut.s.n, as.numeric))
row.names(mut.s.n) = row.names(mut.s.ch)
mut.s.n.s = mut.s.n[rowSums(mut.s.n)/length(s) >= 0.02,]
dim(mut.s.n.s)

CNV.m = as.matrix(CNV.s)
CNV.na10 = CNV.m[rowSums(is.na(CNV.m)) <= 5,]
CNV.na10.sd0.5 = CNV.na10[rowSds(CNV.na10, na.rm=TRUE)>0.5,] 

RNA.m = as.matrix(RNA.s)
RNA.na10 = RNA.m[rowSums(is.na(RNA.m)) <= 5,]
RNA.na10.sd2 = RNA.na10[rowSds(RNA.na10, na.rm=TRUE)>2,] 

PRO.m = as.matrix(PRO.s)
PRO.na10 = PRO.m[rowSums(is.na(PRO.m)) <= 5,]
PRO.na10.sd1 = PRO.na10[rowSds(PRO.na10, na.rm=TRUE)>1,] 

tmut = t(mut.s.n.s)
tCNV = t(CNV.na10.sd0.5)
tRNA = t(RNA.na10.sd2)
tPRO = t(PRO.na10.sd1)

tmut[is.na(tmut)] = 0
tCNV[is.na(tCNV)] = mean(tCNV,na.rm=T)
tRNA[is.na(tRNA)] = mean(tRNA,na.rm=T)
tPRO[is.na(tPRO)] = mean(tPRO,na.rm=T)

sample_aligned = all(rownames(tmut)==rownames(tCNV)) && all(rownames(tmut)==rownames(tRNA)) && all(rownames(tmut)==rownames(tPRO))

##### model tuning #####
if (sample_aligned){
  #tune_icluster = function(){}
  set.seed(123)
  date()
  
  cv.fit = tune.iClusterPlus(cpus=12,dt1=tBRCA_mut,dt2=tBRCA_CNV,dt3=tBRCA_RNA,dt4=tBRCA_PRO,
                             type=c("binomial","gaussian","gaussian","gaussian"),
                             scale.lambda=c(0.04,0.90,0.90,0.90),n.lambda=307,maxiter=20)
  save(cv.fit, file=paste("clusterRdata/cv.fit.k",cancer, "_",k,".Rdata",sep=""))
  
  date()
}

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
