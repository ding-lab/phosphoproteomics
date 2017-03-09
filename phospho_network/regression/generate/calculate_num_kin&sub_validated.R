# Yige Wu @ WashU 2016 Dec
# calculate the number of the validated kinases and substrates

# choose kinase/phosphotase, cancer , significance level -----------------------------------------------
protein <- "kinase"
sig <- 0.05

# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)
library(readr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# input regression result -------------------------------------------------
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))

# print result statistics -------------------------------------------------
## input row numbers of k_s_fit, output how many substrates, phosphosites and kinases, etc
KSP_count <- function(result) {
  row_count <- nrow(result)
  sub_count <- length(unique(result$SUBSTRATE))
  kin_count <- length(unique(result$KINASE))
  kin_brca <- length(unique(result$KINASE[result$Cancer=="BRCA"]))
  kin_ov <- length(unique(result$KINASE[result$Cancer=="OV"]))
  ksp_count <- nrow(unique(result[,c("KINASE","SUBSTRATE","SUB_MOD_RSD")]))
  ksp_brca <- nrow(unique(result[result$Cancer=="BRCA",c("KINASE","SUBSTRATE","SUB_MOD_RSD")]))
  ksp_ov <- nrow(unique(result[result$Cancer=="OV",c("KINASE","SUBSTRATE","SUB_MOD_RSD")]))
  ks_count <- nrow(unique(result[,c("KINASE","SUBSTRATE")]))
  p_count <- nrow(unique(result[,c("SUBSTRATE","SUB_MOD_RSD")]))
  p_brca <- nrow(unique(result[result$Cancer=="BRCA",c("SUBSTRATE","SUB_MOD_RSD")]))
  p_ov <- nrow(unique(result[result$Cancer=="OV",c("SUBSTRATE","SUB_MOD_RSD")]))
  all_count <- data.frame(row_count,ksp_count,ksp_brca,ksp_ov,ks_count,sub_count,kin_count,p_count,kin_brca,kin_ov, p_brca, p_ov)
  return(all_count)
}

time <- Sys.time()
sink(file = paste(baseD,"pan3can_shared_data/analysis_results/tables/No.validated_kinases_substrates_fdr_",sig,"_",time,".txt",sep = ""))
temp <- KSP_count(table_2can)
cat(paste("all regression result has: ", temp$row_count,"lines\n", temp$kin_count," kinases(BRCA: ",temp$kin_brca,";OV: ",temp$kin_ov, ")\n", 
          temp$sub_count," substrates", "\n",
          temp$p_count," phosphosites(BRCA: ",temp$p_brca,";OV: ",temp$p_ov, ")\n", 
          temp$ks_count," substrate:kinase pairs", "\n",
          temp$ksp_count," substrate:phosphosite:kinase pairs(BRCA: ",temp$ksp_brca,";OV: ",temp$ksp_ov, ")\n",  "\n", sep=" "))

temp <- KSP_count(table_2can[table_2can$self,])
cat(paste("cis regression result has: ", temp$row_count,"lines\n", temp$kin_count," kinases(BRCA: ",temp$kin_brca,";OV: ",temp$kin_ov, ")\n", 
          temp$sub_count," substrates", "\n",
          temp$p_count," phosphosites(BRCA: ",temp$p_brca,";OV: ",temp$p_ov, ")\n", 
          temp$ks_count," substrate:kinase pairs", "\n",
          temp$ksp_count," substrate:phosphosite:kinase pairs(BRCA: ",temp$ksp_brca,";OV: ",temp$ksp_ov, ")\n",  "\n", sep=" "))

temp <- KSP_count(table_2can[!table_2can$self,])
cat(paste("trans regression result has: ", temp$row_count,"lines\n", temp$kin_count," kinases(BRCA: ",temp$kin_brca,";OV: ",temp$kin_ov, ")\n", 
          temp$sub_count," substrates", "\n",
          temp$p_count," phosphosites(BRCA: ",temp$p_brca,";OV: ",temp$p_ov, ")\n", 
          temp$ks_count," substrate:kinase pairs", "\n",
          temp$ksp_count," substrate:phosphosite:kinase pairs(BRCA: ",temp$ksp_brca,";OV: ",temp$ksp_ov, ")\n",  "\n", sep=" "))

temp <- KSP_count(table_2can[table_2can$self & table_2can$FDR_pro_kin <= sig,])
cat(paste("cis significant regression result has: ", temp$row_count,"lines\n", temp$kin_count," kinases(BRCA: ",temp$kin_brca,";OV: ",temp$kin_ov, ")\n", 
          temp$sub_count," substrates", "\n",
          temp$p_count," phosphosites(BRCA: ",temp$p_brca,";OV: ",temp$p_ov, ")\n", 
          temp$ks_count," substrate:kinase pairs", "\n",
          temp$ksp_count," substrate:phosphosite:kinase pairs(BRCA: ",temp$ksp_brca,";OV: ",temp$ksp_ov, ")\n",  "\n", sep=" "))

temp <- KSP_count(table_2can[!table_2can$self & table_2can$FDR_pho_kin <= sig,])
cat(paste("trans significant regression result has: ", temp$row_count,"lines\n", temp$kin_count," kinases(BRCA: ",temp$kin_brca,";OV: ",temp$kin_ov, ")\n", 
          temp$sub_count," substrates", "\n",
          temp$p_count," phosphosites(BRCA: ",temp$p_brca,";OV: ",temp$p_ov, ")\n", 
          temp$ks_count," substrate:kinase pairs", "\n",
          temp$ksp_count," substrate:phosphosite:kinase pairs(BRCA: ",temp$ksp_brca,";OV: ",temp$ksp_ov, ")\n",  "\n", sep=" "))

temp <- KSP_count(table_2can[table_2can$self & table_2can$FDR_pro_kin <= sig & table_2can$coef_pro_kin>0,])
cat(paste("cis regression result(FDR<=sig, coef>0) has: ", temp$row_count,"lines\n", temp$kin_count," kinases(BRCA: ",temp$kin_brca,";OV: ",temp$kin_ov, ")\n", 
          temp$sub_count," substrates", "\n",
          temp$p_count," phosphosites(BRCA: ",temp$p_brca,";OV: ",temp$p_ov, ")\n", 
          temp$ks_count," substrate:kinase pairs", "\n",
          temp$ksp_count," substrate:phosphosite:kinase pairs(BRCA: ",temp$ksp_brca,";OV: ",temp$ksp_ov, ")\n",  "\n", sep=" "))

temp <- KSP_count(table_2can[!table_2can$self & table_2can$FDR_pho_kin <= sig & table_2can$coef_pho_kin>0,])
cat(paste("trans regression result(FDR<=sig, coef>0) has: ", temp$row_count,"lines\n", temp$kin_count," kinases(BRCA: ",temp$kin_brca,";OV: ",temp$kin_ov, ")\n", 
          temp$sub_count," substrates", "\n",
          temp$p_count," phosphosites(BRCA: ",temp$p_brca,";OV: ",temp$p_ov, ")\n", 
          temp$ks_count," substrate:kinase pairs", "\n",
          temp$ksp_count," substrate:phosphosite:kinase pairs(BRCA: ",temp$ksp_brca,";OV: ",temp$ksp_ov, ")\n",  "\n", sep=" "))

sink()



# examine important cis-regulated proteins --------------------------------

gene_list <- c("ERBB2","PAK1","RIPK2","BRAF","AKT1","MARK2","PAK2","WNK1")
time <- Sys.time()
sink(file = paste(baseD,"pan3can_shared_data/analysis_results/tables/No.cis_validated_phosphosites_given_genes",time,".txt",sep = ""))
for(gene in gene_list) {
  table_cis_vad <- table_2can[table_2can$self & table_2can$KINASE== gene & table_2can$FDR_pro_kin<=sig,]
  cat(paste("cis validated result for ", gene, ":\n", 
            nrow(table_cis_vad)," phosphosites(BRCA: ",length(which(table_cis_vad$Cancer=="BRCA")),";OV: ",length(which(table_cis_vad$Cancer=="OV")),")\n", sep=" "))
  
}
sink()




