##### find_plot_pQTL.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run pQTL analysis for 3 cancer types and plot the result

##### dependencies #####
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/mutation_impact/")
source("/Users/khuang/bin/LIB_exp.R")
source("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/mutation_impact/mutation_impact.R")

# system("mkdir logs")
# logFile = paste("logs/", date, "_mutation_impact_analysis.log", sep="")
# sink(file=logFile)
#sink(file=NULL)
kinaseList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/2014-12-05_CPTAC_Kinase.MATRIX.v3b5_sheet1_genes.list', header=FALSE, stringsAsFactors = F)
drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv_list.txt_list', header=FALSE, stringsAsFactors = F)
kinome = as.vector(t(kinaseList))
druggable = as.vector(t(drugList))

# function to format CDAP proteome data processed by Kuan
format_pro = function(Pro.m){
  colnames(Pro.m) = sub(".Unshared.Log.Ratio","",colnames(Pro.m)) 
  Pro.m.d = Pro.m[Pro.m$Gene %in% druggable,]
  row.names(Pro.m.d) = Pro.m.d$Gene
  Pro.m.dn = Pro.m.d[,-1]
  return(Pro.m.dn)
}

format_ov = function(Pro.m){
  colnames(Pro.m) = sub(".Unshared.Log.Ratio","",colnames(Pro.m))
  colnames(Pro.m) = sub("^X","",colnames(Pro.m))
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
  colnames(Pro.m) = sub("\\.\\d\\d$","",colnames(Pro.m)) 
  
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

format_mut = function(mut){
  colnames(mut) = sub("TCGA.","", colnames(mut))
  colnames(mut) = sub(".01",".01A", colnames(mut))
  return(mut)
}
##### BRCA #####
### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201510_pancan_somatic/BRCA_proteomic.maf.matrix.txt")
BRCA_mut_s = format_mut(BRCA_mut)
brcaGenes = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "MLL3")

### Proteome ###
BRCA_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Proteome_CDAP.r2/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pro.d = format_pro(BRCA_Pro)
BRCA_diff_exp = find_diff_exp(BRCA_mut_s,BRCA_Pro.d,name="BRCA Proteome")

### Phosphoproteome ###
BRCA_Pho = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Phosphoproteome_CDAP.r2/TCGA_Breast_BI_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pho.d = format_pro(BRCA_Pho)
BRCA_pho_diff_exp = find_diff_exp(BRCA_mut_s,BRCA_Pho.d,name="BRCA Phosphoproteome")

### all levels ###
# do this later, may be more benefitial to do the all outlier table and summarize that instead

##### OV #####
### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201510_pancan_somatic/OV_proteomic.maf.matrix.txt")
OV_mut_s = format_mut(OV_mut)
ovGenes = c("TP53", "NF1", "KRAS", "BRCA1", "CDK12")

### Proteome ###
OV_JHU_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/JHU_Proteome_CDAP.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_PNNL_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/PNNL_Proteome_CDAP.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")

OV_JHU_Pro.d = format_ov(OV_JHU_Pro)
OV_PNNL_Pro.d = format_ov(OV_PNNL_Pro)

OV_JHU_diff_exp = find_diff_exp(OV_mut_s,OV_JHU_Pro.d,name="OV JHU Proteome")
OV_PNNL_diff_exp = find_diff_exp(OV_mut_s,OV_PNNL_Pro.d,name="OV PNNL Proteome")

### Phosphoproteome ###
OV_Pho = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/PNNL_Phosphoproteome_CDAP.r2/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_Pho.d = format_ov(OV_Pho)
OV_Pho_diff_exp = find_diff_exp(OV_mut_s,OV_Pho.d,name="OV PNNL Phosphoproteome")

### Glycoproteome ###
OV_Gly = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/JHU_Glycoproteome_CDAP.r2/TCGA_Ovarian_JHU_Glycoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_Gly.d = format_ov(OV_Gly)
OV_Gly_diff_exp = find_diff_exp(OV_mut_s,OV_Gly.d,name="OV JHU Glycoproteome")

### merging the two proteome? ###
### all levels ###

##### CRC #####
### Mutation matrix ###
CRC_mut = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201510_pancan_somatic/COADREAD_proteomic.maf.matrix.txt")
CRC_mut_s = format_mut(CRC_mut)
crcGenes = c("TP53","KRAS","APC","PIK3CA","SMAD4")
### Proteome ###
CRC_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/CRC/VU_Proteome_CDAP.r2/TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv_hugoified.unshared_spectal_counts.txt")
# spectral count: look into the paper to see how to normalize the count data
CRC_Pro.d = format_crc(CRC_Pro)

CRC_diff_exp = find_diff_exp(CRC_mut_s,CRC_Pro.d,name="CRC Proteome", geneList =crcGenes)

### all levels ###
