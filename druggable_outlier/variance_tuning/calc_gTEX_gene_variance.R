##### plot_variance.R #####
# Kuan-lin Huang @ WashU 2016 Jan
# use gTEX expression data as normal background distribution
# find and plot relative variance in normal samples across genes
# command run: Rscript plot_gene_variance.R > gTEX_rpkm_tissue_gene_stat.txt

##### dependencies #####
setwd("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier/variance_tuning")
source("/Users/khuang/bin/LIB_exp.R")
# library(biomaRt)
# ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
sink("gTEX_rpkm_tissue_gene_stat.txt")
sink("gTEX_rpkm_tissue_gene_stat_log2.txt")


rpkmColClass = c("character", "character", rep("numeric", 8555) )
rpkmFile = read.table(header=TRUE, sep="\t", quote="", row.names=1, skip=2, colClasses=rpkmColClass, file="/Users/khuang/Box Sync/PhD/proteogenomics/gTEX/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct")
phenoTable = read.table(header=TRUE, sep="\t", quote="", row.names=1, colClasses=c("character", "character", "factor", "character", "character", "character", "character"), file="/Users/khuang/Box\ Sync/PhD/proteogenomics/gTEX/GTEx_Data_V6_Annotations_SampleAttributesDS.txt")[,c(1:8)]

# replace row names using gene names instead of 
row.names(rpkmFile) = make.names(rpkmFile[,1], unique =T) # it seems like the multiple entries problem does not affect canonical genes much; some examples: U3, 803 Y_RNA
tissues = unique(as.vector(unlist(phenoTable$SMTS)))

cat(paste("tissue", "gene", "mean", "median","IQR", "SD", "\n", sep="\t"))

# having trouble running from command line can check later
for (tissue in tissues){ # subset sample using each tissue
  samples = row.names(phenoTable[phenoTable$SMTS==tissue,])
  samples2 = gsub("-",".",samples)
  samples_indexes = which(colnames(rpkmFile) %in% samples2)
  data = rpkmFile[,samples_indexes]
  data = log2(data+1)
  # going through gene-sample (row-col) table
  # loop through each gene to find stats
  for (i in 1:nrow(data)){
    data_mean = mean(as.numeric(data[i,]), na.rm=T)
    data_median = median(as.numeric(data[i,]), na.rm=T)
    IQR = quantile(as.numeric(data[i,]), probs=0.75, na.rm=T) - quantile(as.numeric(data[i,]), probs=0.25, na.rm=T) 
    SD = sd(as.numeric(data[i,]), na.rm=T)
    cat(paste(tissue, row.names(data[i,]), data_mean, data_median, IQR, SD, "\n", sep="\t"))
  }
}

sink()
