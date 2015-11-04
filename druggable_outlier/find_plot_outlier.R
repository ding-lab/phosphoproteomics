##### find_plot_outlier.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run outlier analysis for 3 cancer types and plot the result

##### dependencies #####
setwd("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier")
source("/Users/khuang/bin/LIB_exp.R")
source("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier/outlier.R")
system("mkdir logs")
logFile = paste("logs/", date, "_outlier_analysis.log", sep="")
sink(file=logFile)
#sink(file=NULL)
kinaseList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/2014-12-05_CPTAC_Kinase.MATRIX.v3b5_sheet1_genes.list', header=FALSE, stringsAsFactors = F)
drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv_list.txt_list', header=FALSE, stringsAsFactors = F)
kinome = as.vector(t(kinaseList))
druggable = as.vector(t(drugList))

# function to format processed CDAP proteome data 
format_pro = function(Pro.m){
  colnames(Pro.m) = sub(".01A.Unshared.Log.Ratio","",colnames(Pro.m)) 
  colnames(Pro.m) = sub(".01A.1.Unshared.Log.Ratio","",colnames(Pro.m))
  BRCA77 = read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt')
  Pro.m = Pro.m[,colnames(Pro.m) %in% c(colnames(BRCA77),"Gene")]
  Pro.m.d = Pro.m[Pro.m$Gene %in% druggable,]
  row.names(Pro.m.d) = Pro.m.d$Gene
  Pro.m.dn = Pro.m.d[,-1]
  return(Pro.m.dn)
}

format_crc = function(Pro.m){
  row.names(Pro.m) = make.names(Pro.m$Gene, uniq=T)
  Gene = Pro.m$Gene
  Pro.m = as.matrix(Pro.m[,-1])
  colnames(Pro.m) = sub(".Unshared.Spectral.Counts","",colnames(Pro.m)) 
  
  # quantile normalization using function from limma and log2 transformation: Nature CRC proteogenomics 2014
  Pro.mn = normalizeQuantiles(Pro.m,ties=T)
  Pro.mnl = log2(Pro.mn)
  Pro.m.d = Pro.mnl[Gene %in% druggable,]
  
  return(Pro.m.d)
}

# function to normalize CNV based on a log10 scale
normalize_CNV = function(CNV.m){
  colnames(CNV.m) = paste(colnames(CNV.m),".01A", sep="")
  row.names(CNV.m) = sub(" ","", row.names(CNV.m))
  CNV.n.m = as.matrix(CNV.m)
  for (i in 1:nrow(CNV.n.m)){
    CNV.n.m[i,]=log10(CNV.n.m[i,]/mean(CNV.n.m[i,], na.rm=T)) # used log2 in previous versions
  } 
  CNV.n.md = CNV.n.m[row.names(CNV.n.m) %in% druggable,]
  return(CNV.n.md)
}

# function to normalize RSEM
format_RSEM = function(RSEM.m){ # should be normalized using the 75% quantile method already
  RSEM.m.d = RSEM.m[RSEM.m$Hybridization.REF %in% druggable,]
  row.names(RSEM.m.d) = RSEM.m.d$Hybridization.REF
  RSEM.m.d = RSEM.m.d[,-1]
  RSEM.m.d.n = as.matrix(RSEM.m.d)
  # alternative way to normalize that results in less outliers?
#   for (i in 1:nrow(RSEM.m.d.n)){
#     RSEM.m.d.n[i,]=log(RSEM.m.d.n[i,]/mean(RSEM.m.d.n[i,], na.rm=T), base=2)
#   } 
  RSEM.m.d.n = log2(RSEM.m.d.n) # simply log2 transform
  return(RSEM.m.d.n)
}

##### BRCA #####
### CNV ###
BRCA_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_CNV/BRCA_105_CNV.txt")
# "NA " fail, be careful next time about the new space character
BRCA_CNV.d = normalize_CNV(BRCA_CNV)
BRCA_CNV_druggable = find_outlier(BRCA_CNV.d, name = "BRCA druggable CNV")

### RNA ###
BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted.txt",sep=""))
BRCA_RNA.d = BRCA_RNA[row.names(BRCA_RNA) %in% druggable,]
BRCA_RNA_druggable = find_outlier(BRCA_RNA.d, name = "BRCA druggable RNA") # too much outliers

### Proteome ###
BRCA_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Proteome_CDAP.r2/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pro.d = format_pro(BRCA_Pro)
BRCA_Pro_druggable = find_outlier(BRCA_Pro.d, name = "BRCA druggable proteome")

### Phosphoproteome ###
BRCA_Pho = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Phosphoproteome_CDAP.r2/TCGA_Breast_BI_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pho.d = format_pro(BRCA_Pho)
BRCA_Pho_druggable = find_outlier(BRCA_Pho.d, name = "BRCA druggable phosphoproteome")

### RPPA ###
BRCA_RPPA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RPPA/TCGA_Breast_BI_RPPA.tsv")
row.names(BRCA_RPPA) = BRCA_RPPA$Composite.Element.REF
BRCA_RPPA = as.matrix(BRCA_RPPA[,-1])
BRCA_RPPA_outlier = find_outlier(BRCA_RPPA, name = "BRCA RPPA")


### all levels ###
# do this later, may be more benefitial to do the all outlier table and summarize that instead

##### OV #####
### CNV ###
OV_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_CNV/OV_173_CNV.txt")
OV_CNV.d = normalize_CNV(OV_CNV)
OV_CNV_druggable = find_outlier(OV_CNV.d, name = "OV druggable CNV")

### RNA ###
OV_RNA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RNA/TCGA_OV_RSEM.tsv.parsed_hugoified")
OV_RNA.d = format_RSEM(OV_RNA)
OV_RNA_druggable = find_outlier(OV_RNA.d, name = "OV druggable RNA")

### Proteome ###
OV_JHU_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/JHU_Proteome_CDAP.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_PNNL_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/PNNL_Proteome_CDAP.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")

OV_JHU_Pro.d = format_pro(OV_JHU_Pro)
OV_PNNL_Pro.d = format_pro(OV_PNNL_Pro)

OV_JHU_Pro_druggable = find_outlier(OV_JHU_Pro.d, name = "OV JHU druggable proteome")
OV_PNNL_Pro_druggable = find_outlier(OV_PNNL_Pro.d, name = "OV PNNL druggable proteome")

### Phosphoproteome ###
OV_Pho = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/PNNL_Phosphoproteome_CDAP.r2/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_Pho.d = format_pro(OV_Pho)
OV_Pho_druggable = find_outlier(OV_Pho.d, name = "OV PNNL druggable phosphoproteome")

### Glycoproteome ###
OV_Gly = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/JHU_Glycoproteome_CDAP.r2/TCGA_Ovarian_JHU_Glycoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_Gly.d = format_pro(OV_Gly)
OV_Gly_druggable = find_outlier(OV_Gly.d, name = "OV JHU druggable glycoproteome")

### RPPA ###
CRC_RPPA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RPPA/TCGA_OV_RPPA.tsv")
row.names(OV_RPPA) = OV_RPPA$Composite.Element.REF
OV_RPPA = as.matrix(OV_RPPA[,-1])
OV_RPPA_outlier = find_outlier(OV_RPPA, name = "OV RPPA")

### merging the two proteome? ###
### all levels ###

##### CRC #####
### CNV ###
CRC_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_CNV/CRC_88_CNV.txt")
CRC_CNV.d = normalize_CNV(CRC_CNV)
CRC_CNV_druggable = find_outlier(CRC_CNV.d, name = "CRC druggable CNV")

### RNA ###
CRC_RNA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RNA/TCGA_COADREAD_RSEM_combined.tsv.parsed_hugoified")
CRC_RNA.d = format_RSEM(CRC_RNA)
CRC_RNA_druggable = find_outlier(CRC_RNA.d, name = "CRC druggable RNA")

### Proteome ###
CRC_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/CRC/VU_Proteome_CDAP.r2/TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv_hugoified.unshared_spectal_counts.txt")
# spectral count: look into the paper to see how to normalize the count data
CRC_Pro.d = format_crc(CRC_Pro)
CRC_Pro_druggable = find_outlier(CRC_Pro.d, name = "CRC druggable proteome")

### RPPA ###
CRC_RPPA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RPPA/TCGA_COADREAD_RPPA.tsv")
row.names(CRC_RPPA) = CRC_RPPA$Composite.Element.REF
CRC_RPPA = as.matrix(CRC_RPPA[,-1])
CRC_RPPA_outlier = find_outlier(CRC_RPPA, name = "CRC RPPA")
### all levels ###

sink(file=NULL)

##### testing codes#####

if (FALSE){
  m1 = melt(BRCA_Pro_druggable$outlier)
  m2 = melt(BRCA77_druggable$outlier)
  m3 = merge(m1,m2,by=c("Var1","Var2"))
  m3 = m3[complete.cases(m3),]
  nrow(m3[m3$value.x & m3$value.y,])
  nrow(m3[m3$value.x,])
  nrow(m3[m3$value.y,])
  
  m1 = melt(BRCA_Pro_druggable$outlier_zscore)
  m2 = melt(BRCA77_druggable$outlier_zscore)
  m3 = merge(m1,m2,by=c("Var1","Var2"))
  m3 = m3[complete.cases(m3),]
  cor(m3$value.x,m3$value.y)
  nrow(m3[m3$value.x & m3$value.y,])
  nrow(m3[m3$value.x,])
  nrow(m3[m3$value.y,])
  BRCA77_druggable$outlier
}