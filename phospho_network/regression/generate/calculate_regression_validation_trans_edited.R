# Yige Wu @ WashU 2016 Dec
# calculate the validation statistics for regression results


# choose kinase/phosphotase, cancer , significance level -----------------------------------------------
protein <- "kinase"
sig <- 0.05
cancer <- "BRCA"
# cancer <- "OV"

# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)
library(readr)

# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory




# input k_s_table-------------------------------------------------------------------
if ( protein == "kinase" ) {
  ### read in the kinase/substrate table/ phosphorylation data ###
  k_s_table = read.delim(paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
}

if ( protein == "phosphotase" ) {
  ### read in the phosphotase/substrate table/ phosphorylation data ### 
  k_s_table <- read.csv(paste(baseD,"pan3can_shared_data/Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.csv",sep = ""))
  colnames(k_s_table) <- c("Phosphatase_UniProtAC_human","GENE","Substrate_UniProtAC_ref","SUB_GENE","Substrate_Type","DephosphoSite","BioassayType", "PubMed_ID_rev")
}

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

unique_kinase <- unique(k_s_table$GENE)



# different standard of validation for kinase and phosphotase -------------
if ( protein == "kinase" ) {
  valid_type <- "pos_sig"
}
if ( protein == "phosphotase" ) {
  valid_type <- "neg_sig"
}


# input regression processed data -----------------------------------------
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",protein,"_substrate_regression_trans_edited.txt",sep = ""))
table1 <- table_2can[table_2can$Cancer== cancer,]

# initiate for gene-level validation ----------------------------------------------------------------
x <- vector(mode = "numeric", length = length(unique_kinase) )
temp <- data.frame(matrix(rep(x,4), ncol=4, byrow=T))
valid_trans <- cbind(unique_kinase,temp)
colnames(valid_trans) <- c("kinase","pos_sig","pos_insig","neg_sig","neg_insig")
rownames(valid_trans) <- unique_kinase
valid_cis <- valid_trans


# make tables for trans pairs ---------------------------------------------
for( kinase in unique_kinase) {
  temp <- table1[table1$KINASE==kinase,]
  valid_trans[kinase,"pos_sig"] <- length(which(temp$coef_pho_kin>0 & temp$FDR_pho_kin <= sig & !temp$self))
  valid_trans[kinase,"pos_insig"] <- length(which(temp$coef_pho_kin>0 & temp$FDR_pho_kin > sig & !temp$self))
  valid_trans[kinase,"neg_sig"] <- length(which(temp$coef_pho_kin<0 & temp$FDR_pho_kin <= sig & !temp$self))
  valid_trans[kinase,"neg_insig"] <- length(which(temp$coef_pho_kin<0 & temp$FDR_pho_kin > sig & !temp$self))
}
if ( protein == "phosphotase") {
  valid_trans <- valid_trans[,c("kinase","neg_sig","neg_insig","pos_sig","pos_insig")]
}
valid_trans$all_count <- rowSums(valid_trans[,c("pos_sig","pos_insig","neg_sig","neg_insig")])
valid_trans$valid_count <- valid_trans[,valid_type]
valid_trans$valid_ratio <- valid_trans$valid_count/valid_trans$all_count



# write out table for trans ---------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",cancer,'_',protein,"_trans_regression_validation_statistics.txt", sep="")
write.table(valid_trans, file=tn, quote=F, sep = '\t', row.names = FALSE)





# make tables for cis pairs -----------------------------------------------
for( kinase in unique_kinase) {
  temp <- table1[table1$KINASE==kinase,]
  valid_cis[kinase,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & temp$self))
  valid_cis[kinase,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & temp$self))
  valid_cis[kinase,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & temp$self))
  valid_cis[kinase,"neg_insig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & temp$self))
}
if ( protein == "phosphotase") {
  valid_cis <- valid_cis[,c("kinase","neg_sig","neg_insig","pos_sig","pos_insig")]
}
valid_cis$all_count <- rowSums(valid_cis[,c("pos_sig","pos_insig","neg_sig","neg_insig")])
valid_cis$valid_count <- valid_cis[,valid_type]
valid_cis$valid_ratio <- valid_cis$valid_count/valid_cis$all_count

# write out table for cis ---------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",cancer,'_',protein,"_cis_regression_validation_statistics.txt", sep="")
write.table(valid_cis, file=tn, quote=F, sep = '\t', row.names = FALSE)


# initiate for residual-level validation ----------------------------------------
x <- vector(mode = "numeric", length = length(3) )
temp <- data.frame(matrix(rep(x,4), ncol=4, byrow=T))
valid_rsd_trans <- cbind(c("S","T","Y"),temp)
colnames(valid_rsd_trans) <- c("residue","pos_sig","pos_insig","neg_sig","neg_insig")
rownames(valid_rsd_trans) <- c("S","T","Y")
valid_rsd_cis <- valid_rsd_trans
# make table for trans-residue -----------------------------------------------------------
for(rsd in c("S","T","Y") ) {
  temp <- table_2can[grepl(rsd,table_2can$SUB_MOD_RSD),]
  valid_rsd_trans[rsd,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & !temp$self))
  valid_rsd_trans[rsd,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & !temp$self))
  valid_rsd_trans[rsd,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & !temp$self))
  valid_rsd_trans[rsd,"neg_insig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & !temp$self))
}
if ( protein == "phosphotase") {
  valid_rsd_trans <- valid_rsd_trans[,c("residue","neg_sig","neg_insig","pos_sig","pos_insig")]
}
valid_rsd_trans$all_count <- rowSums(valid_rsd_trans[,c("pos_sig","pos_insig","neg_sig","neg_insig")])
valid_rsd_trans$valid_count <- valid_rsd_trans[,valid_type]
valid_rsd_trans$valid_ratio <- valid_rsd_trans$valid_count/valid_rsd_trans$all_count

# write out table for trans-residual --------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",cancer,'_',protein,"_trans_regression_validation_statistics_for_residual.txt", sep="")
write.table(valid_rsd_trans, file=tn, quote=F, sep = '\t', row.names = FALSE)

# make table for cis-residue -------------------------------------------------------------
## cis-residue
for(rsd in c("S","T","Y") ) {
  temp <- table_2can[grepl(rsd,table_2can$SUB_MOD_RSD),]
  valid_rsd_cis[rsd,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & temp$self))
  valid_rsd_cis[rsd,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & temp$self))
  valid_rsd_cis[rsd,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & temp$self))
  valid_rsd_cis[rsd,"neg_insig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & temp$self))
}
if ( protein == "phosphotase") {
  valid_rsd_cis <- valid_rsd_cis[,c("residue","neg_sig","neg_insig","pos_sig","pos_insig")]
}
valid_rsd_cis$all_count <- rowSums(valid_rsd_cis[,c("pos_sig","pos_insig","neg_sig","neg_insig")])
valid_rsd_cis$valid_count <- valid_rsd_cis[,valid_type]
valid_rsd_cis$valid_ratio <- valid_rsd_cis$valid_count/valid_rsd_cis$all_count

# write out table for trans-residual --------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",cancer,'_',protein,"_cis_regression_validation_statistics_for_residual.txt", sep="")
write.table(valid_rsd_cis, file=tn, quote=F, sep = '\t', row.names = FALSE)


# output validated kinase-substrate pairs ---------------------------------
table <- table_2can
table$protein <- "kinase"

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
    
    tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/enrichment/",cancer,"_",iscis,"_kinase_result_background.txt", sep="")
    write.table(valid_table$gene, file=tn, quote=F, sep = '\t', row.names = F, col.names = F)
    
    valid_table <- valid_table[order(valid_table$valid_count, decreasing = T),]
    tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/enrichment/",cancer,"_",iscis,"_kinase_result_valid_count_sig_",sig,".txt", sep="")
    write.table(valid_table[valid_table$valid_count>0,c("gene","valid_count")], file=tn, quote=F, sep = '\t', row.names = F, col.names = F)
    
  }
}




