##### phosphosite_summary.R #####
# Kuan-lin Huang @ WashU 2016 Feb
# cross-correlation between different phosphosites

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/data_summary/")
source("/Users/khuang/bin/LIB_exp.R")

brca_file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
ov_file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
brca_pho = read.table(row.names=1, header=TRUE, sep="\t", file= brca_file)
ov_pho = read.table(row.names=1, header=TRUE, sep="\t", file= ov_file)

num_sites = length(unique(c(row.names(brca_pho),row.names(ov_pho))))
brca_genes = gsub("(.*?)\\.(.).*","\\1",row.names(brca_pho))
ov_genes = gsub("(.*?)\\.(.).*","\\1",row.names(ov_pho))
num_genes = length(unique(c(brca_genes,ov_genes)))
cat("number of unique genes characterized:", num_genes)
cat("number of unique sites characterized:", num_sites)

# BRCA
non_unique = gsub("\\.[0-9]+","",row.names(brca_pho))
brca_residue = gsub("(.*)\\.(.).*","\\2",non_unique)
BRCA = table(brca_residue)
BRCA

brca_pho_10NA = brca_pho[rowSums(!is.na(brca_pho))>=10,]
non_unique = gsub("\\.[0-9]+","",row.names(brca_pho_10NA))
brca_residue_10NA = gsub("(.*)\\.(.).*","\\2",non_unique)
BRCA10NA = table(brca_residue_10NA)

# OV
non_unique = gsub("\\.[0-9]+","",row.names(ov_pho))
OV_residue = gsub("(.*)\\.(.).*","\\2",non_unique)
OV = table(OV_residue)
OV

OV_pho_10NA = ov_pho[rowSums(!is.na(ov_pho))>=10,]
non_unique = gsub("\\.[0-9]+","",row.names(OV_pho_10NA))
OV_residue_10NA = gsub("(.*)\\.(.).*","\\2",non_unique)
OV10NA = table(OV_residue_10NA)

# both
residue_table = rbind(BRCA,BRCA10NA,OV,OV10NA)
residue_table_m = melt(residue_table)

fn = paste(pd, 'phosphosite_brca.ov_sty_summary.pdf',sep ="_")
p = ggplot(residue_table_m,aes(x=Var1, y=value, fill=Var2))
p = p + geom_bar(binwidth = 0.05, stat="identity") + theme_bw() + theme_nogrid()
p = p + labs(x = "Data type", y="# of phosphosites")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
ggsave(file=fn, height=5, width=8, useDingbats=FALSE)