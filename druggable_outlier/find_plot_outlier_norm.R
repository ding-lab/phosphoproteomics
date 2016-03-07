##### find_plot_outlier.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# do another smad normalization for all data type first
# run outlier analysis for 3 cancer types and plot the result

##### dependencies #####
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"
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

extract_exp_genes = function(expFile){
  RNA.c = read.table(row.names=1, header=TRUE, sep="\t", file=expFile)
  RNA.c.d = RNA.c[row.names(RNA.c) %in% druggable,]
  #expressed at higher than 10 counts in at least 70% of samples in the specific tissue
  exp_gene = row.names(RNA.c.d)[rowSums(RNA.c.d >= 10) > dim(RNA.c)[2]*0.7] 
  return(exp_gene)
}
##### BRCA #####

### extract tissue-specific expressed gene list from RNA data before normalization ###
BRCA_RNA.c = paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted.txt",sep="")
BRCA_exp_gene = extract_exp_genes(BRCA_RNA.c)

### CNV ###
BRCA_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_CNV_formatted_normalized.txt",sep=""))
# "NA " fail, be careful next time about the new space character
BRCA_CNV.d = BRCA_CNV[row.names(BRCA_CNV) %in% BRCA_exp_gene,]
BRCA_CNV_druggable = find_outlier(BRCA_CNV.d, name = "BRCA druggable CNV normalized")

### RNA ###
BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))
BRCA_RNA.d = BRCA_RNA[row.names(BRCA_RNA) %in% BRCA_exp_gene,]
BRCA_RNA_druggable = find_outlier(BRCA_RNA.d, name = "BRCA druggable RNA normalized")

### Proteome ###
BRCA_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
BRCA_Pro.d = BRCA_Pro[row.names(BRCA_Pro) %in% BRCA_exp_gene,]
BRCA_Pro_druggable = find_outlier(BRCA_Pro.d, name = "BRCA druggable proteome normalized")

### Phosphoproteome ###
BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
BRCA_Pho.d = BRCA_Pho[row.names(BRCA_Pho) %in% BRCA_exp_gene,]
BRCA_Pho_druggable = find_outlier(BRCA_Pho.d, name = "BRCA druggable phosphoproteome normalized")

### RPPA ### 
BRCA_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_RPPA_formatted.txt",sep=""))
gene = sub("_.*","",row.names(BRCA_RPPA))
gene = sub(" .*","",gene)
BRCA_RPPA.d = BRCA_RPPA[gene %in% BRCA_exp_gene,]
BRCA_RPPA_druggable = find_outlier(BRCA_RPPA.d, name = "BRCA druggable RPPA")

### all levels ###
# do this later, may be more benefitial to do the all outlier table and summarize that instead

##### OV #####
### extract tissue-specific expressed gene list from RNA data before normalization ###
OV_RNA.c = paste(baseD,"pan3can_shared_data/OV/OV_mRNA_formatted.txt",sep="")
OV_exp_gene = extract_exp_genes(OV_RNA.c)

### CNV ###
OV_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_CNV_formatted_normalized.txt",sep=""))
OV_CNV.d = OV_CNV[row.names(OV_CNV) %in% OV_exp_gene,]
OV_CNV_druggable = find_outlier(OV_CNV.d, name = "OV druggable CNV normalized")

### RNA ###
OV_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_mRNA_formatted_normalized.txt",sep=""))
OV_RNA.d = OV_RNA[row.names(OV_RNA) %in% OV_exp_gene,]
OV_RNA_druggable = find_outlier(OV_RNA.d, name = "OV druggable RNA normalized")

### Proteome ###
OV_JHU_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_PRO_formatted_normalized.txt",sep=""))
OV_PNNL_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))

OV_JHU_Pro.d = OV_JHU_Pro[row.names(OV_JHU_Pro) %in% OV_exp_gene,]
OV_PNNL_Pro.d = OV_PNNL_Pro[row.names(OV_PNNL_Pro) %in% OV_exp_gene,]

OV_JHU_Pro_druggable = find_outlier(OV_JHU_Pro.d, name = "OV JHU druggable proteome normalized")
OV_PNNL_Pro_druggable = find_outlier(OV_PNNL_Pro.d, name = "OV PNNL druggable proteome normalized")

### Phosphoproteome ###
OV_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))
OV_Pho.d = OV_Pho[row.names(OV_Pho) %in% OV_exp_gene,]
OV_Pho_druggable = find_outlier(OV_Pho.d, name = "OV PNNL druggable phosphoproteome normalized")

### Glycoproteome ###
OV_Gly = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_GLY_formatted_normalized.txt",sep=""))
OV_Gly.d = OV_Gly[row.names(OV_Gly) %in% OV_exp_gene,]
#OV_Gly_druggable = find_outlier(OV_Gly.d, name = "OV JHU druggable glycoproteome normalized")
## glyco throwing errors for now

### RPPA ### 
OV_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_RPPA_formatted.txt",sep=""))
gene = sub("_.*","",row.names(OV_RPPA))
gene = sub(" .*","",gene)
OV_RPPA.d = OV_RPPA[gene %in% OV_exp_gene,]
OV_RPPA_druggable = find_outlier(OV_RPPA.d, name = "OV druggable RPPA")

### merging the two proteome? ###
### all levels ###

##### CRC #####
CRC_RNA.c = paste(baseD,"pan3can_shared_data/CRC/CRC_mRNA_formatted.txt",sep="")
CRC_exp_gene = extract_exp_genes(CRC_RNA.c)

# CRC_pro.c = paste(baseD,"pan3can_shared_data/CRC/CRC_PRO_formatted.txt",sep="")
# CRC_exp_pro = extract_exp_genes(CRC_pro.c) # only 7 proteins...

### CNV ###
CRC_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_CNV_formatted_normalized.txt",sep=""))
CRC_CNV.d = CRC_CNV[row.names(CRC_CNV) %in% CRC_exp_gene,]
CRC_CNV_druggable = find_outlier(CRC_CNV.d, name = "CRC druggable CNV normalized")

### RNA ###
CRC_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_mRNA_formatted_normalized.txt",sep=""))
CRC_RNA.d = CRC_RNA[row.names(CRC_RNA) %in% CRC_exp_gene,]
CRC_RNA_druggable = find_outlier(CRC_RNA.d, name = "CRC druggable RNA normalized")

### Proteome ###
CRC_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_PRO_formatted_normalized.txt",sep=""))
# spectral count: look into the paper to see how to normalize the count data
CRC_Pro.d = CRC_Pro[row.names(CRC_Pro) %in% CRC_exp_gene,]
CRC_Pro_druggable = find_outlier(CRC_Pro.d, name = "CRC druggable proteome normalized")

### RPPA ### 
CRC_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_RPPA_formatted.txt",sep=""))
gene = sub("_.*","",row.names(CRC_RPPA))
gene = sub(" .*","",gene)
CRC_RPPA.d = CRC_RPPA[gene %in% CRC_exp_gene,]
CRC_RPPA_druggable = find_outlier(CRC_RPPA.d, name = "CRC druggable RPPA")

### all levels ###

sink(file=NULL)
