##### preprocessPhosphosites.R #####
# Kuan-lin Huang @ WashU 2016 Feb
# preprocess phosphosite data of BI BRCA and PNNL OV

library(limma)
library(biomaRt)
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/"

normalize_by_sample = function(m){
  m = as.matrix(m)
  m[!is.finite(m)] = NA
  for (i in 1:ncol(m)){
    m[,i] = (m[,i] - mean(m[,i], na.rm=T))/sd(m[,i], na.rm=T)
  }
  return(m)
}

brca_s_f = read.table(sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/sample_list/final_list/TCGA_Breast_BI_Proteome_CDAP_sample_77unimodal_u.list')
brca_s = as.vector(t(brca_s_f))
##### convert the refseq protein name to hgnc ID #####
#ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")

BRCA_Psite = read.table(header=T,na.strings="",sep="\t","/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv")
#sites = sub("(.*):(.*)","\\2",BRCA_Psite$Phosphosite)
gene_sites = paste(BRCA_Psite$Gene,BRCA_Psite$Phosphosite,sep=":")
row.names(BRCA_Psite) = gene_sites#make.names(gene_sites, unique=T)
BRCA_P_abbr = BRCA_Psite[, !(colnames(BRCA_Psite) %in% c("Phosphosite","Peptide","Gene","Organism"))]
colnames(BRCA_P_abbr) = gsub("\\.Log.Ratio","",colnames(BRCA_P_abbr))
colnames(BRCA_P_abbr) = gsub("\\.","-",colnames(BRCA_P_abbr))
colnames(BRCA_P_abbr) = gsub("-01A-.","-01A",colnames(BRCA_P_abbr))
BRCA_P_abbr = BRCA_P_abbr[,brca_s]
#write.table(BRCA_P_abbr,col.names=NA, quote=F, sep = '\t', file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev.tsv")
BRCA_P_abbr_n = normalize_by_sample(BRCA_P_abbr)
#write.table(BRCA_P_abbr_n,col.names=NA, quote=F, sep = '\t', file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv")
BRCA_P_abbr_n_10NA = BRCA_P_abbr_n[rowSums(!is.na(BRCA_P_abbr_n))>=10,]
#write.table(BRCA_P_abbr_n_10NA, col.names=NA, quote=F, sep = '\t', file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv")

OV_Psite = read.table(header=T,na.strings="",sep="\t","/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Ovarian_PNNL_Phosphoproteome/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv")
sites = sub("(.*):(.*)","\\2",OV_Psite$Phosphosite)
gene_sites = paste(OV_Psite$Gene,OV_Psite$Phosphosite,sep=":")#sites,sep=".")
row.names(OV_Psite) = gene_sites#make.names(gene_sites, unique=T)
OV_P_abbr = OV_Psite[, !(colnames(OV_Psite) %in% c("Phosphosite","Peptide","Gene","Organism"))]
colnames(OV_P_abbr) = gsub("\\.Log.Ratio","",colnames(OV_P_abbr))
#write.table(OV_P_abbr,col.names=NA, quote=F, sep = '\t', file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Ovarian_PNNL_Phosphoproteome/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev.tsv")
OV_P_abbr_n = normalize_by_sample(OV_P_abbr)
#write.table(OV_P_abbr_n,col.names=NA, quote=F, sep = '\t', file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Ovarian_PNNL_Phosphoproteome/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv")
OV_P_abbr_n_10NA = OV_P_abbr_n[rowSums(!is.na(OV_P_abbr_n))>=10,]
#write.table(OV_P_abbr_n_10NA, col.names=NA, quote=F, sep = '\t', file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Ovarian_PNNL_Phosphoproteome/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA.tsv")

num_BRCA_sites = nrow(BRCA_P_abbr_n )
cat("Number of BRCA phosphosites:", num_BRCA_sites, "\n")
num_BRCA_sites = nrow(BRCA_P_abbr_n_10NA )
cat("Number of BRCA phosphosites with more than 10 observations:", num_BRCA_sites, "\n")

num_OV_sites = nrow(OV_P_abbr_n )
cat("Number of OV phosphosites:", num_OV_sites, "\n")
num_OV_sites = nrow(OV_P_abbr_n_10NA )
cat("Number of OV phosphosites with more than 10 observations:", num_OV_sites, "\n")

num_shared_sites = sum(row.names(BRCA_P_abbr_n_10NA) %in% row.names(OV_P_abbr_n_10NA))
cat("Number of shared phosphosites with more than 10 observations each:", num_OV_sites, "\n")

#BRCA_Psite$refseq_peptide = gsub("\\..*","",BRCA_Psite$Phosphosite)
# #getBM(filters="refseq_peptide", attributes="external_gene_id", values=refseq, mart=ensembl)
# mapTab = getBM(attributes = c("refseq_peptide","hgnc_symbol"), filters = "refseq_peptide", values = BRCA_Psite$refseq_peptide, mart = ensembl, uniqueRows=FALSE)
# #dupRows = which(duplicated(mapTab[,2]))#union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
# BRCA_Psite_m = merge(BRCA_Psite, mapTab, by="refseq_peptide", all.x=T)
# BRCA_Psite_m[BRCA_Psite_m$Gene != BRCA_Psite_m$hgnc_symbol,]
# # didn't make much difference, don't use
