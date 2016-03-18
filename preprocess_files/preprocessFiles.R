##### preprocessFiles.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# preprocess omic files for the pan3can proteomic project
# make sure all the sample names match up
# create normalized files for some inputs

library(limma)
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/"

setwd("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/preprocess_files/")
source("/Users/khuang/bin/LIB_exp.R")
# read the master sample list
brca_s_f = read.table(sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/sample_list/final_list/TCGA_Breast_BI_Proteome_CDAP_sample_77unimodal_u.list')
brca_s = as.vector(t(brca_s_f))
brca_s2 = gsub("-",".",as.vector(t(brca_s_f)))
ov_s_f = read.table(sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/sample_list/final_list/TCGA_Ovarian_Combined_Proteome_CDAP_sample.list')
ov_s = as.vector(t(ov_s_f))
crc_s_f = read.table(sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/sample_list/final_list/TCGA_Colon_VU_Proteome_CDAP_sample_u.list')
crc_s = as.vector(t(crc_s_f))
crc_s2 = gsub("-",".",as.vector(t(crc_s_f)))

pan3can_s = c(brca_s,ov_s,crc_s)

format_brca = function(Pro.m){
  Pro.m = Pro.m[-c(1,2,3),]
  colnames(Pro.m) = sub(".1.Unshared.Log.Ratio","",colnames(Pro.m)) 
  colnames(Pro.m) = sub(".Unshared.Log.Ratio","",colnames(Pro.m)) 
  # extract only the unimodal samples
  BRCA77 = read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt')
  colnames(BRCA77) = paste(colnames(BRCA77),".01A",sep="")
  Pro.m.d = Pro.m[,colnames(Pro.m) %in% c(colnames(BRCA77),"Gene")]
  row.names(Pro.m.d) = make.names(Pro.m.d[,"Gene"], unique=T)
  Pro.m.dn = as.matrix(Pro.m.d[,-1])
  return(Pro.m.dn)
}

format_ov = function(Pro.m){
  Pro.m = Pro.m[-c(1,2,3),]
  colnames(Pro.m) = sub(".Unshared.Log.Ratio","",colnames(Pro.m))
  colnames(Pro.m) = sub("^X","",colnames(Pro.m))
  row.names(Pro.m) = make.names(Pro.m$Gene, unique=T)
  Pro.m = Pro.m[,-1]
  return(Pro.m)
}

normalize_by_sample = function(m){
  m = as.matrix(m)
  m[!is.finite(m)] = NA
  for (i in 1:ncol(m)){
    m[,i] = (m[,i] - mean(m[,i], na.rm=T))/sd(m[,i], na.rm=T)
  }
  return(m)
}

format_crc = function(Pro.m){
  row.names(Pro.m) = make.names(Pro.m$Gene, uniq=T)
  Gene = Pro.m$Gene
  Pro.m = as.matrix(Pro.m[,-1])
  colnames(Pro.m) = sub(".Unshared.Spectral.Counts","",colnames(Pro.m)) 
  colnames(Pro.m) = sub("\\.\\d\\d$","",colnames(Pro.m)) 
  return(Pro.m)
}

normalize_crc = function(m){
  m = as.matrix(m)
  # quantile normalization using function from limma and log2 transformation: Nature CRC proteogenomics 2014
  Pro.mn = normalizeQuantiles(m,ties=T)
  Pro.mnl = log2(Pro.mn+1)
  return(Pro.mnl)
}

format_CNV = function(CNV.m){
  colnames(CNV.m) = paste(colnames(CNV.m),".01A", sep="")
  row.names(CNV.m) = sub(" ","", row.names(CNV.m))
  CNV.n.m = as.matrix(CNV.m)
  return(CNV.n.m)
}

normalize_CNV = function(CNV.m){
  CNV.n.m = as.matrix(CNV.m)
  for (i in 1:nrow(CNV.n.m)){
    CNV.n.m[i,]=log2(CNV.n.m[i,]/mean(CNV.n.m[i,], na.rm=T))
  } 
  return(CNV.n.m)
}

format_RSEM = function(RSEM.m){ # should have been normalized using the 75% quantile across samples already
  row.names(RSEM.m) = make.names(RSEM.m$Hybridization.REF, unique=T)
  RSEM.m = RSEM.m[,-1]
  RSEM.m = as.matrix(RSEM.m)
  #RSEM.m.d.n = log2(RSEM.m.d.n+1) # simple log2 transform
  return(RSEM.m)
}

# format_RPPA = function(m){ #TODO
#   row.names(m) = make.names(RSEM.m$Hybridization.REF, unique=T)
#   RSEM.m = RSEM.m[,-1]
#   RSEM.m = as.matrix(RSEM.m)
#   #RSEM.m.d.n = log2(RSEM.m.d.n+1) # simple log2 transform
#   return(RSEM.m)
# }

format_mut = function(mut){
  colnames(mut) = sub("TCGA.","", colnames(mut))
  colnames(mut) = sub("\\.01",".01A", colnames(mut))
  colnames(mut) = sub("-01","-01A", colnames(mut))
  return(mut)
}

format_clin = function(clin){
  colnames(clin) = sub("TCGA-","", colnames(clin))
  colnames(clin) = paste(colnames(clin),"-01A",sep="")
  return(clin)
}

##### Pre-processing #####

RPPA = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/collaborations/TCPA_2015-10-30/TCGA-PANCAN16-RBN-Trans-Gene_sampleProbe_matrix.tsv")
colnames(RPPA) = sub("TCGA.","",colnames(RPPA))
colnames(RPPA) = paste(colnames(RPPA),".01A",sep="")
##### BRCA #####

### Proteome ###
BRCA_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Proteome_CDAP.r2/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pro.f = format_brca(BRCA_Pro) # retained 77 samples
##write.table(BRCA_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PRO_formatted.txt",sep=""))
BRCA_Pro.n = normalize_by_sample(BRCA_Pro.f)
#write.table(BRCA_Pro.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))

### RPPA ###
BRCA_RPPA = RPPA[,colnames(RPPA) %in% colnames(BRCA_Pro.f)]
#write.table(BRCA_RPPA, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_RPPA_formatted.txt",sep=""))

### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201511_pancan_somatic/BRCA_proteomic_dbFilter.maf.matrix_trv.txt")
BRCA_mut_s = format_mut(BRCA_mut)
BRCA_mut_s = BRCA_mut_s[,colnames(BRCA_mut_s) %in% brca_s2]
#write.table(BRCA_mut_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_SOMATIC_formatted.txt",sep=""))

BRCA_mut_aa = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201511_pancan_somatic/BRCA_proteomic_dbFilter.maf.matrix.txt")
BRCA_mut_aa_s = format_mut(BRCA_mut_aa)
BRCA_mut_aa_s = BRCA_mut_aa_s[,colnames(BRCA_mut_s) %in% brca_s2]
#write.table(BRCA_mut_aa_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_SOMATIC_formatted_amino_acid.txt",sep=""))

### CNV ###
BRCA_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_CNV/BRCA_105_CNV.txt")
# "NA " fail, be careful next time about the new space character
BRCA_CNV.f = format_CNV(BRCA_CNV)
BRCA_CNV.f = BRCA_CNV.f[,colnames(BRCA_CNV.f) %in% colnames(BRCA_Pro.f)]
#write.table(BRCA_CNV.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_CNV_formatted.txt",sep=""))
BRCA_CNV.n = normalize_CNV(BRCA_CNV.f)
#write.table(BRCA_CNV.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_CNV_formatted_normalized.txt",sep=""))

### RNA ###
BRCA_RNA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_RNA/TCGA_Breast_BI_RSEM.tsv.parsed_hugoified")
BRCA_RNA.f = format_RSEM(BRCA_RNA)
BRCA_RNA.f = BRCA_RNA.f[,colnames(BRCA_RNA.f) %in% colnames(BRCA_Pro.f)]
#write.table(BRCA_RNA.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_mRNA_formatted.txt",sep=""))
BRCA_RNA.n = log2(BRCA_RNA.f+1)
#write.table(BRCA_RNA.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))

### Phosphoproteome ###
BRCA_Pho = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Phosphoproteome_CDAP.r2/TCGA_Breast_BI_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pho.f = format_brca(BRCA_Pho)
#write.table(BRCA_Pho.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PHO_formatted.txt",sep=""))
BRCA_Pho.n = normalize_by_sample(BRCA_Pho.f)
#write.table(BRCA_Pho.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))

### Normalize Phospho by proteome ### 
# match genes and samples
b_genes = row.names(BRCA_Pho.f)[row.names(BRCA_Pho.f) %in% row.names(BRCA_Pro.f)]
BRCA_Pro.f.g = BRCA_Pro.f[b_genes,]
BRCA_Pho.f.g = BRCA_Pho.f[b_genes,]
BRCA_Pho_by_Pro.f = BRCA_Pho.f.g - BRCA_Pro.f.g
#write.table(BRCA_Pho_by_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PHO_by_PRO_formatted.txt",sep=""))
BRCA_Pho_by_Pro.n = normalize_by_sample(BRCA_Pho_by_Pro.f)
#write.table(BRCA_Pho_by_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PHO_by_PRO_formatted_normalized.txt",sep=""))

### Clinical ###
BRCA_Clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file="/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/clinical_association/20150916_20150821v2_stddate_clinical_all_cancer_type_pan7943/BRCA_clinical.txt")
tBRCA_Clin = t(BRCA_Clin)
tBRCA_Clin_f = format_clin(tBRCA_Clin)
tBRCA_Clin_f_s = tBRCA_Clin_f[,colnames(tBRCA_Clin_f) %in% brca_s]
#write.table(tBRCA_Clin_f_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_clinical.txt",sep=""))

BRCA_Clin_short = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_brca_201506/80\ Tumors\ Data/clin.info.CPTAC2.qc.miss.25.perct.short.txt")
tBRCA_Clin_short = t(BRCA_Clin_short)
colnames(tBRCA_Clin_short) = toupper(colnames(tBRCA_Clin_short))
tBRCA_Clin_short_f = format_clin(tBRCA_Clin_short)
tBRCA_Clin_short_f_s = tBRCA_Clin_short_f[,colnames(tBRCA_Clin_short_f) %in% brca_s]
#write.table(tBRCA_Clin_short_f_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_clinical_summary.txt",sep=""))


##### OV #####

### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201511_pancan_somatic/OV_proteomic_dbFiletred.maf.matrix_trv.txt")
OV_mut_s = format_mut(OV_mut)
##OV_mut_s = OV_mut_s[,-ncol(OV_mut_s)]
#write.table(OV_mut_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_SOMATIC_formatted.txt",sep=""))

OV_mut_aa = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201511_pancan_somatic/OV_proteomic_dbFiletred.maf.matrix.txt")
OV_mut_aa_s = format_mut(OV_mut_aa)
##OV_mut_aa_s = OV_mut_aa_s[,-ncol(OV_mut_aa_s)]
#write.table(OV_mut_aa_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_SOMATIC_formatted_amino_acid.txt",sep=""))

### CNV ###
OV_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_CNV/OV_173_CNV.txt")
OV_CNV.f = format_CNV(OV_CNV)
#write.table(OV_CNV.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_CNV_formatted.txt",sep=""))
OV_CNV.n = normalize_CNV(OV_CNV.f)
#write.table(OV_CNV.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_CNV_formatted_normalized.txt",sep=""))

### RPPA ###
RPPA2 = RPPA
colnames(RPPA2) = gsub("\\.","-",colnames(RPPA2))
OV_RPPA = RPPA2[,colnames(RPPA2) %in% ov_s]
#write.table(OV_RPPA, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_RPPA_formatted.txt",sep=""))

### RNA ###
OV_RNA = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_RNA/TCGA_OV_RSEM.tsv.parsed_hugoified")
colnames(OV_RNA)[1]="Hybridization.REF"
OV_RNA.f = format_RSEM(OV_RNA)
#write.table(OV_RNA.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_mRNA_formatted.txt",sep=""))
OV_RNA.n = log2(OV_RNA.f+1)
#write.table(OV_RNA.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_mRNA_formatted_normalized.txt",sep=""))

### Proteome ###
OV_JHU_Pro = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_proteome_CDAP_r2/OV/JHU_Proteome_CDAP.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_PNNL_Pro = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_proteome_CDAP_r2/OV/PNNL_Proteome_CDAP.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")

OV_JHU_Pro.f = format_ov(OV_JHU_Pro) 
#write.table(OV_JHU_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_PRO_formatted.txt",sep=""))
OV_JHU_Pro.n = normalize_by_sample(OV_JHU_Pro.f)
#write.table(OV_JHU_Pro.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_PRO_formatted_normalized.txt",sep=""))

# exclude control
OV_JHU_Pro.f_u = OV_JHU_Pro.f[,-c(grep(".*CONTROL.*",colnames(OV_JHU_Pro.f)))]
#write.table(OV_JHU_Pro.f_u, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_PRO_formatted_noControl.txt",sep=""))
OV_JHU_Pro.n_u = OV_JHU_Pro.n[,-c(grep(".*CONTROL.*",colnames(OV_JHU_Pro.n)))]
#write.table(OV_JHU_Pro.n_u, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_PRO_formatted_normalized_noControl.txt",sep=""))

OV_PNNL_Pro.f = format_ov(OV_PNNL_Pro) 
#write.table(OV_PNNL_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PRO_formatted.txt",sep=""))
OV_PNNL_Pro.n = normalize_by_sample(OV_PNNL_Pro.f)
#write.table(OV_PNNL_Pro.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))

OV_PNNL_Pro.n_noJHU = OV_PNNL_Pro.n[,!(colnames(OV_PNNL_Pro.n) %in% colnames(OV_JHU_Pro.n_u))]
OV_Pro_no_overlap = merge(OV_JHU_Pro.n_u, OV_PNNL_Pro.n_noJHU, by="row.names",all=F)
row.names(OV_Pro_no_overlap) = OV_Pro_no_overlap$Row.names
OV_Pro_no_overlap = OV_Pro_no_overlap[,-1]
#write.table(OV_Pro_no_overlap, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_merged_PRO_noOverlap_formatted_normalized.txt",sep=""))

colnames(OV_JHU_Pro.n_u) = paste(colnames(OV_JHU_Pro.n_u),"JHU",sep="_")
colnames(OV_PNNL_Pro.n) = paste(colnames(OV_PNNL_Pro.n),"PNNL",sep="_")
OV_Pro.n = merge(OV_JHU_Pro.n_u, OV_PNNL_Pro.n, by="row.names")
row.names(OV_Pro.n) = OV_Pro.n$Row.names
OV_Pro.n = OV_Pro.n[,-1]
#write.table(OV_Pro.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_merged_PRO_formatted_normalized.txt",sep=""))


# not an issue for PNNL set
# OV_PNNL_Pro.f_u = OV_PNNL_Pro.f[,-c(grep(".*CONTROL.*",colnames(OV_PNNL_Pro.f)))]
# #write.table(OV_JHU_Pro.f_u, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_PRO_formatted_noControl.txt",sep=""))
# OV_PNNL_Pro.n_u = OV_PNNL_Pro.n[,-c(grep(".*CONTROL.*",colnames(OV_JHU_Pro.n)))]
# #write.table(OV_JHU_Pro.n_u, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_PRO_formatted_normalized_noControl.txt",sep=""))

### Phosphoproteome ###
OV_PNNL_Pho = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_proteome_CDAP_r2/OV/PNNL_Phosphoproteome_CDAP.r2/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_PNNL_Pho.f = format_ov(OV_PNNL_Pho)
colnames(OV_JHU_Pro.f) %in% ov_s
#write.table(OV_PNNL_Pho.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PHO_formatted.txt",sep=""))
OV_PNNL_Pho.n = normalize_by_sample(OV_PNNL_Pho.f)
#write.table(OV_PNNL_Pho.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))

### Normalize Phospho by proteome ### 
# match genes and samples
o_genes = row.names(OV_PNNL_Pho.f)[row.names(OV_PNNL_Pho.f) %in% row.names(OV_PNNL_Pro.f)]
o_samples = colnames(OV_PNNL_Pho.f)[colnames(OV_PNNL_Pho.f) %in% colnames(OV_PNNL_Pro.f)]
OV_PNNL_Pro.f.g = OV_PNNL_Pro.f[o_genes,o_samples]
OV_PNNL_Pho.f.g = OV_PNNL_Pho.f[o_genes,o_samples]
OV_PNNL_Pho_by_Pro.f = OV_PNNL_Pho.f.g - OV_PNNL_Pro.f.g
#write.table(OV_PNNL_Pho_by_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PHO_by_PRO_formatted.txt",sep=""))
OV_PNNL_Pho_by_Pro.n = normalize_by_sample(OV_PNNL_Pho_by_Pro.f)
#write.table(OV_PNNL_Pho_by_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PHO_by_PRO_formatted_normalized.txt",sep=""))

### Glycoproteome ###
OV_JHU_Gly = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_proteome_CDAP_r2/OV/JHU_Glycoproteome_CDAP.r2/TCGA_Ovarian_JHU_Glycoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_JHU_Gly.f = format_ov(OV_JHU_Gly) 
#write.table(OV_JHU_Gly.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_GLY_formatted.txt",sep=""))
OV_JHU_Gly.n = normalize_by_sample(OV_JHU_Gly.f)
#write.table(OV_JHU_Gly.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_GLY_formatted_normalized.txt",sep=""))

### Clinical ###
# if (FALSE){
# OV_Clin = read.table(row.names=1, quote = "", header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/clinical_association/20150916_20150821v2_stddate_clinical_all_cancer_type_pan7943/OV_clinical.txt")
# tOV_Clin = t(OV_Clin)
# tOV_Clin_f = format_clin(tOV_Clin)
# tOV_Clin_f_s = tOV_Clin_f[,colnames(tOV_Clin_f) %in% ov_s]
# #write.table(tOV_Clin_f_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_clinical.txt",sep=""))
# }

OV_Clin2 = read.table(row.names=1, quote = "", header=TRUE, sep="\t",file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201603_clin_downloads/nationwidechildrens.org_OV.bio.Level_2.0.42.0/nationwidechildrens.org_clinical_patient_ov.txt")
row.names(OV_Clin2) = OV_Clin2[,1]
OV_Clin2.r = OV_Clin2[-(1:2),-(1:2)]
tOV_Clin2 = t(OV_Clin2.r)
tOV_Clin_f2 = format_clin(tOV_Clin2)
tOV_Clin_f_s2 = tOV_Clin_f2[,colnames(tOV_Clin_f2) %in% ov_s]
#write.table(tOV_Clin_f_s2, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_clinical.txt",sep=""))

## not working for some reason
# OV_Clin_short = read.table(row.names=1, header=TRUE, sep=",", quote = "", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/References/TCGA_downloads/TCGA_2011Nat_OV/2010-09-11380C-Table_S1.2.csv")
# tOV_Clin_short = t(OV_Clin_short)
# colnames(tOV_Clin_short) = toupper(colnames(tOV_Clin_short))
# tOV_Clin_short_f = format_clin(tOV_Clin_short)
# tOV_Clin_short_f_s = tOV_Clin_short_f[,colnames(tOV_Clin_short_f) %in% OV_s]

### merging the two proteome? ###

##### CRC #####
### Mutation matrix ###
CRC_mut = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201511_pancan_somatic/COADREAD_proteomic_dbFilter.maf.matrix_trv.txt")
CRC_mut_s = format_mut(CRC_mut)
##CRC_mut_s = CRC_mut_s[,-ncol(CRC_mut_s)] # tumor sample barcode
write.table(CRC_mut_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_SOMATIC_formatted.txt",sep=""))

CRC_mut_aa = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201511_pancan_somatic/COADREAD_proteomic_dbFilter.maf.matrix.txt")
CRC_mut_aa_s = format_mut(CRC_mut_aa)
##CRC_mut_aa_s = CRC_mut_aa_s[,colnames(CRC_mut_s) %in% CRC_s2]
write.table(CRC_mut_aa_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_SOMATIC_formatted_amino_acid.txt",sep=""))

### CNV ###
CRC_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_CNV/CRC_88_CNV.txt")
CRC_CNV.f = format_CNV(CRC_CNV)
#write.table(CRC_CNV.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_CNV_formatted.txt",sep=""))
CRC_CNV.n = normalize_CNV(CRC_CNV.f)
#write.table(CRC_CNV.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_CNV_formatted_normalized.txt",sep=""))

### RNA ###
CRC_RNA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_RNA/TCGA_COADREAD_RSEM_combined.tsv.parsed_hugoified")
CRC_RNA.f = format_RSEM(CRC_RNA)
CRC_RNA.f = CRC_RNA.f[,colnames(CRC_RNA.f) %in% crc_s2]
#write.table(CRC_RNA.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_mRNA_formatted.txt",sep=""))
CRC_RNA.n = log2(CRC_RNA.f+1)
#write.table(CRC_RNA.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_mRNA_formatted_normalized.txt",sep=""))

### Proteome ###
CRC_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201507_pancan_proteome_CDAP_r2/CRC/VU_Proteome_CDAP.r2/TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv_hugoified.unshared_spectal_counts.txt")
CRC_Pro.f = format_crc(CRC_Pro)
CRC_Pro.f = CRC_Pro.f[,-c(ncol(CRC_Pro.f),which(duplicated(colnames(CRC_Pro.f))))]
#write.table(CRC_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_PRO_formatted.txt",sep=""))
CRC_Pro.f = normalize_crc(CRC_Pro.f)
#write.table(CRC_Pro.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_PRO_formatted_normalized.txt",sep=""))

### RPPA ###
CRC_RPPA = RPPA[,colnames(RPPA) %in% colnames(CRC_Pro.f)]
#write.table(CRC_RPPA, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_RPPA_formatted.txt",sep=""))

### Clinical ###
# if (FALSE){
# COAD_Clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file="/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/clinical_association/20150916_20150821v2_stddate_clinical_all_cancer_type_pan7943/COAD_clinical.txt")
# READ_Clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file="/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/clinical_association/20150916_20150821v2_stddate_clinical_all_cancer_type_pan7943/READ_clinical.txt")
# COAD_Clin = COAD_Clin[,intersect(colnames(COAD_Clin),colnames(READ_Clin))]
# READ_Clin = READ_Clin[,intersect(colnames(COAD_Clin),colnames(READ_Clin))]
# CRC_Clin = rbind(COAD_Clin,READ_Clin)
# tCRC_Clin = t(CRC_Clin)
# tCRC_Clin_f = format_clin(tCRC_Clin)
# tCRC_Clin_f_s = tCRC_Clin_f[,colnames(tCRC_Clin_f) %in% crc_s] # no samples extracted for now...
# #write.table(tCRC_Clin_f_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_clinical.txt",sep="")) 
# }

COAD_Clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201603_clin_downloads/nationwidechildrens.org_COAD.bio.Level_2.0.43.0/nationwidechildrens.org_clinical_patient_coad.txt")
READ_Clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201603_clin_downloads/nationwidechildrens.org_READ.bio.Level_2.0.42.0/nationwidechildrens.org_clinical_patient_read.txt")
row.names(COAD_Clin) = COAD_Clin[,1]
COAD_Clin = COAD_Clin[-(1:2),-(1:2)]
row.names(READ_Clin) = READ_Clin[,1]
READ_Clin = READ_Clin[-(1:2),-(1:2)]
COAD_Clin = COAD_Clin[,intersect(colnames(COAD_Clin),colnames(READ_Clin))]
READ_Clin = READ_Clin[,intersect(colnames(COAD_Clin),colnames(READ_Clin))]
CRC_Clin = rbind(COAD_Clin,READ_Clin)
tCRC_Clin = t(CRC_Clin)
tCRC_Clin_f = format_clin(tCRC_Clin)
tCRC_Clin_f_s = tCRC_Clin_f[,colnames(tCRC_Clin_f) %in% crc_s]
#write.table(tCRC_Clin_f_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_clinical.txt",sep="")) 

##### pan3can #####
### Somatic ###
Somatic_2 = merge(BRCA_mut_s,OV_mut_s,by="row.names", all=T)
row.names(Somatic_2) = Somatic_2[,1]
Somatic_2 = Somatic_2[,-1]
pan3can_somatic = merge(Somatic_2, CRC_mut_s,by="row.names", all=T)
row.names(pan3can_somatic) = pan3can_somatic[,1]
pan3can_somatic = pan3can_somatic[,-1]
pan3can_somatic[is.na(pan3can_somatic)]="wt"
#write.table(pan3can_somatic, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"pan3can/pan3can_SOMATIC_formatted.txt",sep="")) 

### CNV ###
CNV_2 = merge(BRCA_CNV.f,OV_CNV.f,by="row.names", all=F)
row.names(CNV_2) = CNV_2[,1]
CNV_2 = CNV_2[,-1]
pan3can_CNV = merge(CNV_2, CRC_CNV.f,by="row.names", all=F)
row.names(pan3can_CNV) = pan3can_CNV[,1]
pan3can_CNV = pan3can_CNV[,-1]
#write.table(pan3can_CNV, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"pan3can/pan3can_CNV_formatted.txt",sep="")) 
pan3can_CNV.n = normalize_CNV(pan3can_CNV)
#write.table(pan3can_CNV.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"pan3can/pan3can_CNV_formatted_normalized.txt",sep="")) 

### RNA ###
RNA_2 = merge(BRCA_RNA.f,OV_RNA.f,by="row.names", all=F)
row.names(RNA_2) = RNA_2[,1]
RNA_2 = RNA_2[,-1]
pan3can_RNA = merge(RNA_2, CRC_RNA.f,by="row.names", all=F)
row.names(pan3can_RNA) = pan3can_RNA[,1]
pan3can_RNA = pan3can_RNA[,-1]
#write.table(pan3can_RNA, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"pan3can/pan3can_RNA_formatted.txt",sep="")) 
pan3can_RNA.n = log2(pan3can_RNA+1)
#write.table(pan3can_RNA.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"pan3can/pan3can_RNA_formatted_normalized.txt",sep="")) 

### RPPA ###
RPPA_2 = merge(BRCA_RPPA,OV_RPPA,by="row.names", all=T)
row.names(RPPA_2) = RPPA_2[,1]
RPPA_2 = RPPA_2[,-1]
pan3can_RPPA = merge(RPPA_2, CRC_RPPA,by="row.names", all=T)
row.names(pan3can_RPPA) = pan3can_RPPA[,1]
pan3can_RPPA = pan3can_RPPA[,-1]
#write.table(pan3can_RPPA, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"pan3can/pan3can_RPPA_formatted.txt",sep="")) 

### PRO ###
## experimental ##
## do a simple sample-by-sample normalization ##
PRO_2 = merge(BRCA_Pro.n,OV_Pro_no_overlap,by="row.names", all=F)
row.names(PRO_2) = PRO_2[,1]
PRO_2 = PRO_2[,-1]
pan3can_PRO = merge(PRO_2, CRC_Pro.n,by="row.names", all=F)
row.names(pan3can_PRO) = pan3can_PRO[,1]
pan3can_PRO = pan3can_PRO[,-1]
pan3can_PRO.n = normalize_by_sample(pan3can_PRO)
#write.table(pan3can_PRO.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"pan3can/pan3can_PRO_formatted_normalized.txt",sep="")) 

### Clinical ###
brca_s_f$cancer = "BRCA"
brca_s_f$V1 = gsub("-",".",brca_s_f$V1)
ov_s_f$cancer = "OV"
crc_s_f$cancer = "CRC"
crc_s_f$V1 = gsub("-",".",crc_s_f$V1)
ctype1 = rbind(brca_s_f,ov_s_f)
ctype = rbind(ctype1, crc_s_f)
t_ctype = t(ctype)
#write.table(t_ctype, col.names=F, quote=F, sep = '\t', file=paste(baseD,"pan3can/pan3can_ctype.txt",sep="")) 


##### find distributions of RNA exp data #####

# not normalized RNA
pan3can_RNA.m = melt(as.matrix(pan3can_RNA[,c(1:10)]))

p = ggplot(pan3can_RNA.m,aes(as.numeric(value), fill=Var2))
p = p + geom_density(alpha=0.6) + theme_bw() + theme_nogrid() + xlim(0,25)
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
# use 5 as a cut-off for genes not expressed in the sample

# not normalized CRC Pro
