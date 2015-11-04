##### preprocessFiles.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# preprocess omic files for the pan3can proteomic project
# make sure all the sample names match up
# create normalized files for some inputs

library(limma)
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/"
setwd(baseD)

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
  Pro.mnl = log2(Pro.mn)
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
  colnames(mut) = sub(".01",".01A", colnames(mut))
  return(mut)
}


##### Pre-processing #####

RPPA = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/collaborations/TCPA_2015-10-30/TCGA-PANCAN16-RBN-Trans-Gene_sampleProbe_matrix.tsv")
colnames(RPPA) = sub("TCGA.","",colnames(RPPA))
colnames(RPPA) = paste(colnames(RPPA),".01A",sep="")
##### BRCA #####

### Proteome ###
BRCA_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Proteome_CDAP.r2/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pro.f = format_brca(BRCA_Pro) # retained 77 samples
write.table(BRCA_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PRO_formatted.txt",sep=""))
BRCA_Pro.n = normalize_by_sample(BRCA_Pro.f)
write.table(BRCA_Pro.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))

### RPPA ###
BRCA_RPPA = RPPA[,colnames(RPPA) %in% colnames(BRCA_Pro.f)]
write.table(BRCA_RPPA, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_RPPA_formatted.txt",sep=""))

### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201510_pancan_somatic/BRCA_proteomic.maf.matrix.txt")
BRCA_mut_s = format_mut(BRCA_mut)
BRCA_mut_s = BRCA_mut_s[,colnames(BRCA_mut_s) %in% colnames(BRCA_Pro.f)]
write.table(BRCA_mut_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_SOMATIC_formatted.txt",sep=""))

### CNV ###
BRCA_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_CNV/BRCA_105_CNV.txt")
# "NA " fail, be careful next time about the new space character
BRCA_CNV.f = format_CNV(BRCA_CNV)
BRCA_CNV.f = BRCA_CNV.f[,colnames(BRCA_CNV.f) %in% colnames(BRCA_Pro.f)]
write.table(BRCA_CNV.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_CNV_formatted.txt",sep=""))
BRCA_CNV.n = normalize_CNV(BRCA_CNV.f)
write.table(BRCA_CNV.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_CNV_formatted_normalized.txt",sep=""))

### RNA ###
BRCA_RNA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RNA/TCGA_Breast_BI_RSEM.tsv.parsed_hugoified")
BRCA_RNA.f = format_RSEM(BRCA_RNA)
BRCA_RNA.f = BRCA_RNA.f[,colnames(BRCA_RNA.f) %in% colnames(BRCA_Pro.f)]
write.table(BRCA_RNA.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_mRNA_formatted.txt",sep=""))
BRCA_RNA.n = log2(BRCA_RNA.f+1)
write.table(BRCA_RNA.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))

### Phosphoproteome ###
BRCA_Pho = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Phosphoproteome_CDAP.r2/TCGA_Breast_BI_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pho.f = format_brca(BRCA_Pho)
write.table(BRCA_Pho.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PHO_formatted.txt",sep=""))
BRCA_Pho.n = normalize_by_sample(BRCA_Pho.f)
write.table(BRCA_Pho.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))

##### OV #####

### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201510_pancan_somatic/OV_proteomic.maf.matrix.txt")
OV_mut_s = format_mut(OV_mut)
write.table(OV_mut_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_SOMATIC_formatted.txt",sep=""))

### CNV ###
OV_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_CNV/OV_173_CNV.txt")
OV_CNV.f = format_CNV(OV_CNV)
write.table(OV_CNV.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_CNV_formatted.txt",sep=""))
OV_CNV.n = normalize_CNV(OV_CNV.f)
write.table(OV_CNV.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_CNV_formatted_normalized.txt",sep=""))

### RPPA ###
RPPA2 = RPPA
colnames(RPPA2) = sub("\\.","-",colnames(RPPA2))
OV_RPPA = RPPA2[,colnames(RPPA2) %in% colnames(OV_CNV.f)]
write.table(OV_RPPA, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_RPPA_formatted.txt",sep=""))

### RNA ###
OV_RNA = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RNA/TCGA_OV_RSEM.tsv.parsed_hugoified")
colnames(OV_RNA)[1]="Hybridization.REF"
OV_RNA.f = format_RSEM(OV_RNA)
write.table(OV_RNA.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_mRNA_formatted.txt",sep=""))
OV_RNA.n = log2(OV_RNA.f+1)
write.table(OV_RNA.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_mRNA_formatted_normalized.txt",sep=""))

### Proteome ###
OV_JHU_Pro = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/JHU_Proteome_CDAP.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_PNNL_Pro = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/PNNL_Proteome_CDAP.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")

OV_JHU_Pro.f = format_ov(OV_JHU_Pro) 
write.table(OV_JHU_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_PRO_formatted.txt",sep=""))
OV_JHU_Pro.n = normalize_by_sample(OV_JHU_Pro.f)
write.table(OV_JHU_Pro.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_PRO_formatted_normalized.txt",sep=""))

OV_PNNL_Pro.f = format_ov(OV_PNNL_Pro) 
write.table(OV_PNNL_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PRO_formatted.txt",sep=""))
OV_PNNL_Pro.n = normalize_by_sample(OV_PNNL_Pro.f)
write.table(OV_PNNL_Pro.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))

### Phosphoproteome ###
OV_PNNL_Pho = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/PNNL_Phosphoproteome_CDAP.r2/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_PNNL_Pho.f = format_ov(OV_PNNL_Pho)
write.table(OV_PNNL_Pho.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PHO_formatted.txt",sep=""))
OV_PNNL_Pho.n = normalize_by_sample(OV_PNNL_Pho.f)
write.table(OV_PNNL_Pho.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))

### Glycoproteome ###
OV_JHU_Gly = read.table(header=TRUE, sep="\t", check.names=F, file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/JHU_Glycoproteome_CDAP.r2/TCGA_Ovarian_JHU_Glycoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_JHU_Gly.f = format_ov(OV_JHU_Gly) 
write.table(OV_JHU_Gly.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_GLY_formatted.txt",sep=""))
OV_JHU_Gly.n = normalize_by_sample(OV_JHU_Gly.f)
write.table(OV_JHU_Gly.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"OV/OV_JHU_GLY_formatted_normalized.txt",sep=""))

### merging the two proteome? ###

##### CRC #####
### Mutation matrix ###
CRC_mut = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201510_pancan_somatic/COADREAD_proteomic.maf.matrix.txt")
CRC_mut_s = format_mut(CRC_mut)
write.table(CRC_mut_s, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_SOMATIC_formatted.txt",sep=""))

### CNV ###
CRC_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_CNV/CRC_88_CNV.txt")
CRC_CNV.f = format_CNV(CRC_CNV)
write.table(CRC_CNV.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_CNV_formatted.txt",sep=""))
CRC_CNV.n = normalize_CNV(CRC_CNV.f)
write.table(CRC_CNV.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_CNV_formatted_normalized.txt",sep=""))

### RNA ###
CRC_RNA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RNA/TCGA_COADREAD_RSEM_combined.tsv.parsed_hugoified")
CRC_RNA.f = format_RSEM(CRC_RNA)
write.table(CRC_RNA.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_mRNA_formatted.txt",sep=""))
CRC_RNA.n = log2(CRC_RNA.f+1)
write.table(CRC_RNA.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_mRNA_formatted_normalized.txt",sep=""))

### Proteome ###
CRC_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/CRC/VU_Proteome_CDAP.r2/TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv_hugoified.unshared_spectal_counts.txt")
CRC_Pro.f = format_crc(CRC_Pro)
write.table(CRC_Pro.f, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_PRO_formatted.txt",sep=""))
CRC_Pro.n = normalize_crc(CRC_Pro.f)
write.table(CRC_Pro.n, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_PRO_formatted_normalized.txt",sep=""))

### RPPA ###
CRC_RPPA = RPPA[,colnames(RPPA) %in% colnames(CRC_Pro.f)]
write.table(CRC_RPPA, col.names=NA, quote=F, sep = '\t', file=paste(baseD,"CRC/CRC_RPPA_formatted.txt",sep=""))

##### TO-DO #####
##### RPPA #####
