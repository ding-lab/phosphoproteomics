##### iClusterPlus_preprocessData.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# cross level clustering between CNV, RNA, and Proteome data
# preprocess data for later model tuning and final analysis


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

name = "BRCA"

cat("############################")
cat(date())
cat(paste("Tuning iCluster model in BRCA of cluster number: ",k)) 

# functions
unfactorize = function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

##### MAIN CODE #####

# BRCA data: load
mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",name,"/",name,"_SOMATIC_formatted.txt",sep=""))
CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",name,"/",name,"_CNV_formatted_normalized.txt",sep=""))
RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",name,"/",name,"_mRNA_formatted_normalized.txt",sep=""))
PRO = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",name,"/",name,"_PRO_formatted_normalized.txt",sep=""))
#Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/",name,"/",name,"_PHO_formatted_normalized.txt",sep=""))
clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file=paste(baseD,"pan3can_shared_data/",name,"/",name,"_clinical_summary.txt",sep=""))
# sync samples
s = intersect(colnames(mut),colnames(CNV))
s = intersect(s,colnames(RNA))
s = intersect(s,colnames(PRO))
s = intersect(s,colnames(clin))

mut.s = mut[,s]
CNV.s = CNV[,s]
RNA.s = RNA[,s]
PRO.s = PRO[,s]
clin.s = clin[,s]

#preProcess = function(mut =, cnv = , rna =, pro = )
# pre-clustering
# preliminary feature selection to roughly decide number of features to go into each process

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
#dim(mut.s.n.s)

CNV.m = as.matrix(CNV.s)
CNV.na10 = CNV.m[rowSums(is.na(CNV.m)) <= 5,]
CNV.na10.sd0.5 = CNV.na10[rowSds(CNV.na10, na.rm=TRUE)>0.5,] 

RNA.m = as.matrix(RNA.s)
RNA.na10 = RNA.m[rowSums(is.na(RNA.m)) <= 5,]
RNA.na10.sd2 = RNA.na10[rowSds(RNA.na10, na.rm=TRUE)>2,] 

PRO.m = as.matrix(PRO.s)
PRO.na10 = PRO.m[rowSums(is.na(PRO.m)) <= 5,]
PRO.na10.sd1 = PRO.na10[rowSds(PRO.na10, na.rm=TRUE)>1,] 

clin.s.s = clin.s[c(1:3),]

tmut = t(mut.s.n.s)
tCNV = t(CNV.na10.sd0.5)
tRNA = t(RNA.na10.sd2)
tPRO = t(PRO.na10.sd1)
tclin = t(clin.s.s)

# iCluster can't tolerate NA
tmut[is.na(tmut)] = 0
tCNV[is.na(tCNV)] = mean(tCNV,na.rm=T)
tRNA[is.na(tRNA)] = mean(tRNA,na.rm=T)
tPRO[is.na(tPRO)] = mean(tPRO,na.rm=T)

sample_aligned = all(rownames(tmut)==rownames(tCNV)) && all(rownames(tmut)==rownames(tRNA)) && all(rownames(tmut)==rownames(tPRO))

save.image(file=paste("clusterRdata/",pd,"_",name,".Rdata",sep="")