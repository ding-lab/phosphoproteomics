## output kinase/phosphotase/substrates in significant regression results, then do pathway enrichment analysis
## divide kinase and phosphotase and substrate, divide cis and trans, divide BRCA and OV

# choose kinase/phosphotase, cancer , significance level, model -----------------------------------------------
protein <- "kinase"
sig <- 0.05

# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# load data for kinase/phosphotase ------------------------------------
# table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",protein,"_substrate_regression_trans_edited.txt",sep = ""))
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))

table <- table_2can
table$protein <- "kinase"

# output gene list and score --------------------------------------------------------
for(cancer in c("BRCA","OV")) {
  for (iscis in c("cis","trans") ) {
    if (iscis == "cis") {
      var <- "pro_kin"
    }
    if (iscis == "trans") {
      var <- "pho_kin"
    }
    
    valid_table <- data.frame()
    temp <- table[table$Cancer==cancer & table$SELF==iscis & table$protein=="kinase",]
    for( g in unique(union(temp$KINASE,temp$SUBSTRATE)) ) {
      # no matter g is kinase or substrate, as long as it's in a record, and coef > 0 and FDR <= .., it's counted validated
      table_g <- temp[temp$KINASE==g | temp$SUBSTRATE==g,]
      all_count <- nrow(table_g)
      valid_count <- length(which(table_g[,paste("FDR_",var,sep = "")] <=sig & table_g[,paste("coef_",var,sep = "")] > 0))
      valid_table[g,"gene"] <- g
      valid_table[g,"valid_count"] <- valid_count
      valid_table[g,"valid_ratio"] <- valid_count/all_count
    }
    valid_table <- valid_table[order(valid_table$valid_ratio, decreasing = T),]
    tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/enrichment/",cancer,"_",iscis,"_kinase_result_valid_ratio_sig_",sig,".txt", sep="")
    write.table(valid_table[valid_table$valid_ratio>0,c("gene","valid_ratio")], file=tn, quote=F, sep = '\t', row.names = F, col.names = F)
    
    tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/enrichment/",cancer,"_",iscis,"_kinase_result_ratio_ranked_background.txt", sep="")
    write.table(valid_table$gene, file=tn, quote=F, sep = '\t', row.names = F, col.names = F)
    
    valid_table <- valid_table[order(valid_table$valid_count, decreasing = T),]
    tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/enrichment/",cancer,"_",iscis,"_kinase_result_valid_count_sig_",sig,".txt", sep="")
    write.table(valid_table[valid_table$valid_count>0,c("gene","valid_count")], file=tn, quote=F, sep = '\t', row.names = F, col.names = F)
    
    tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/enrichment/",cancer,"_",iscis,"_kinase_result_count_ranked_background.txt", sep="")
    write.table(valid_table$gene, file=tn, quote=F, sep = '\t', row.names = F, col.names = F)
  }
}


