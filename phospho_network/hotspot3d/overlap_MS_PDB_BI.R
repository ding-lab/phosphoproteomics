# Yige Wu @ WashU 2017 Jan
# overlap regression result with PDB active sites

# directory and library ---------------------------------------------------
# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"

library(stringr)
library(ggplot2)
library(readr)

setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# choose kinase or phosphotase, significance level, outlier threshold and least sample number-------------------------
sig <- 0.05 # significance level
out_thres <- 1.5 #threshold for outlier
least_samples <- 5# least number of samples with complete data for each model
protein <- "kinase"


# input PDB file and processed regression result ---------------------------------------
PDB_CPTAC_sites = read.table(paste(baseD,"pan3can_shared_data/analysis_results/tables/PDB_BRCA77_sites.pairwise",sep = ""), sep="\t", header=F, quote="")
table_HUMAN_cis_sig = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_cis_sig_fam.txt",sep = ""))
table_HUMAN_trans_sig = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_trans_sig_fam.txt",sep = ""))
# colnames(PDB_CPTAC_sites) = c("gene1","phosphosite1","feature1","other1","note1",
#                               "gene2","phosphosite2","feature2","other2","note2",
#                               "distance","PDB","P")
colnames(PDB_CPTAC_sites) = c("gene1","chr1","start1","stop1","phosphosite1","other1","note1","feature1","noteother1",
                              "gene2","chr2","start2","stop2","phosphosite2","other2","note2","feature2","noteother2",
                              "note","PDB")

# cis overlapping significantly regulated site to pairwise features
PDB_CPTAC_sites$pair = paste(PDB_CPTAC_sites$gene1,PDB_CPTAC_sites$gene1,gsub("p.","",PDB_CPTAC_sites$phosphosite1),sep=":")
table_HUMAN_cis_sig_overlap1 = merge(table_HUMAN_cis_sig,PDB_CPTAC_sites,by="pair")
PDB_CPTAC_sites$pair = paste(PDB_CPTAC_sites$gene2,PDB_CPTAC_sites$gene2,gsub("p.","",PDB_CPTAC_sites$phosphosite2),sep=":")
table_HUMAN_cis_sig_overlap2 = merge(table_HUMAN_cis_sig,PDB_CPTAC_sites,by="pair")
table_HUMAN_cis_sig_overlap = rbind(table_HUMAN_cis_sig_overlap1,table_HUMAN_cis_sig_overlap2)
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_regression_cis_sig_fam_overlapPDB.txt",sep = "")
write.table(table_HUMAN_cis_sig_overlap, file=tn, quote=F, sep = '\t', row.names = FALSE)


# trans overlapping significantly regulated site to pairwise features
table_HUMAN_trans_sig$sub_site = paste(table_HUMAN_trans_sig$SUBSTRATE,gsub("(.*:)","",table_HUMAN_trans_sig$pair),sep=":")
PDB_CPTAC_sites$sub_site = paste(PDB_CPTAC_sites$gene1,gsub("p.","",PDB_CPTAC_sites$phosphosite1),sep=":")
table_HUMAN_trans_sig_overlap1 = merge(table_HUMAN_trans_sig,PDB_CPTAC_sites,by="sub_site")
PDB_CPTAC_sites$sub_site = paste(PDB_CPTAC_sites$gene2,gsub("p.","",PDB_CPTAC_sites$phosphosite2),sep=":")
table_HUMAN_trans_sig_overlap2 = merge(table_HUMAN_trans_sig,PDB_CPTAC_sites,by="sub_site")
table_HUMAN_trans_sig_overlap = rbind(table_HUMAN_trans_sig_overlap1,table_HUMAN_trans_sig_overlap2)
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_regression_trans_sig_fam_overlapPDB.txt",sep = "")
write.table(table_HUMAN_trans_sig_overlap, file=tn, quote=F, sep = '\t', row.names = FALSE)