# Kuan @ WashU March 2017
# analysis on regulated sites
# keep only the sites with significant kinase regulators after multiple kinase model

# choose kinase/phosphotase, cancer , significance level, model -----------------------------------------------
protein <- "kinase"
# protein <- "phosphotase"
sig <- 0.05
cancer <- "BRCA"
out_thres <- 1.5

# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)

# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

# input regression processed data -----------------------------------------
multi_reg <- read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_kinase_substrate_regression_multi.txt",sep = ""))

multi_reg$sub_res = paste(multi_reg$SUBSTRATE,multi_reg$SUB_MOD_RSD,sep="_")
cat("Number of tested substrate sites: ",length(unique(multi_reg$sub_res)),"\n")
cat("Number of tested substrate-kinase in multi-kinase model: ",dim(multi_reg)[1],"\n")
multi_reg_sig_sites = multi_reg$sub_res[multi_reg$FDR_pho_kin < sig]

multi_reg_sig = multi_reg[multi_reg$sub_res %in% multi_reg_sig_sites,]
cat("Number of remaining significant substrate sites: ",length(unique(multi_reg_sig_sites)),"\n")
cat("Number of remaining tested substrate-kinase in multi-kinase model: ",dim(multi_reg_sig)[1],"\n")

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_kinase_substrate_regression_multi_sig_residue.txt",sep = "")
write.table(multi_reg_sig, file=tn, quote=F, sep = '\t', row.names = FALSE)