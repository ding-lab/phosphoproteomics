# Yige Wu @ WashU 2017 Feb

# choose kinase/phosphotase, significance level, outlier threshold and least sample number-------------------------
protein <- "kinase"
# protein <- "phosphotase"


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


# input processed data for single-model and multi-model -----------------------------------
table_single <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))
table_multi <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/multi_",protein,"_substrate_regression.txt",sep = ""))


# qqplot ------------------------------------------------------------------
multi_single_overlap = function(c) {
  table_multi_temp <- table_multi[table_multi$Cancer==c & !table_multi$self]
  table_can_temp <- table_single[table_single$Cancer == c & !table_single$self,]
  
  merge_id <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","pair")
  merge_cols <- c(merge_id,"scale_FDR","FDR_pho_kin","coef_pho_kin","model","Size")
  table_overlap <- merge(table_multi_temp[,merge_cols],table_can_temp[,merge_cols], by = merge_id)
  return(table_overlap)
}


#for (cancer in c("BRCA","OV")) {
for (cancer in "BRCA") {
  for (iscis in "trans") { # loop around cis and trans
    if (iscis == "trans") {
      var <- "pho_kin"
    }
    fdr_var <- paste("FDR_",var,sep = "")
    coef_var <- paste("coef_",var,sep = "")
    
    table_multi_temp <- table_multi[table_multi$Cancer==cancer & table_multi$SELF == iscis,]
    table_can_temp <- table_single[table_single$Cancer == cancer & table_single$SELF == iscis,]
    table_multi_temp$scale_FDR <- -log10(table_multi_temp$FDR_pho_kin)
    table_can_temp$scale_FDR <- -log10(table_can_temp$FDR_pho_kin)
    
    merge_id <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","pair")
    merge_cols <- c(merge_id,"scale_FDR","FDR_pho_kin","coef_pho_kin","model","Size")
    table_overlap <- merge(table_multi_temp[,merge_cols],table_can_temp[,merge_cols], by = merge_id)
    
    cat(paste("multi-model overlap with single-model ",cancer," ",iscis," regression result:\n",
              nrow(table_overlap)," kinase:substrate:phosphosite pairs (",
              nrow(table_multi_temp)," PDX pairs vs ",
              nrow(table_can_temp)," BRCA pairs)\n\n",
              sep = ""))
    
    p = ggplot(table_overlap,aes(x=scale_FDR.x, y=scale_FDR.y))
    p = p + geom_point(alpha=0.05)
    #p = p + geom_text(aes(label= ifelse( ( scale_FDR.x >= -log10(0.01) & scale_FDR.y >= -log10(0.01)) , pair, NA)),size=2,alpha=0.5)
    p = p + geom_text(aes(label= ifelse( ( scale_FDR.y >= -log10(0.05) & scale_FDR.x < -log10(0.05) ) , pair, NA)),size=2,alpha=0.2)
    #p = p + geom_text(aes(label= ifelse( top & scale_FDR.x >= -log10(0.05) , pair, NA)),size=2,alpha=0.5)
    p = p + theme_bw() #+ theme_nogrid()
    p = p + geom_abline(slope=1)
    p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
    p = p + labs(x = "-log(FDR) in multi-model", y="-log(FDR) in single-model")
    p = p + xlim(0,20) + ylim(0,20)
    p
  }
}

# substrate <- "APC";rsd <- "S2129"
# substrate <- "RAF1";rsd <- "S29"
substrate <- "NPM1";rsd <- "S70"

old <- table_single[table_single$SUBSTRATE==substrate & table_single$SUB_MOD_RSD==rsd & !table_single$self,]
new <- table_multi[table_multi$SUBSTRATE==substrate & table_multi$SUB_MOD_RSD==rsd,]

change10 <- table_overlap[table_overlap$FDR_pho_kin.y<sig & table_overlap$FDR_pho_kin.x>=sig,]
phosphosites10 <- unique(change10[,c("SUBSTRATE","SUB_MOD_RSD")])
substrates <- as.vector(phosphosites10$SUBSTRATE)
rsds <- as.vector(phosphosites10$SUB_MOD_RSD)
for (i in 1:nrow(phosphosites10)) {
  substrate <- substrates[i]
  rsd <- rsds[i]
  #table_old <- table_single[table_single$SUBSTRATE==substrate & table_single$SUB_MOD_RSD==rsd & !table_single$self,]
  table_new <- table_multi[table_multi$SUBSTRATE==substrate & table_multi$SUB_MOD_RSD==rsd & table_multi$Cancer==cancer,]
  kin_sig <- table_new$KINASE[table_new$FDR_pho_kin<sig]
  if (length(kin_sig) > 0) {
    # print(phosphosites10[i,])
    for (k in kin_sig) {
      k_subs <- unique(table_single$SUBSTRATE[table_single$KINASE==k & !table_single$self])
      print(k_subs)
      # if (length(k_subs)==1) {
      #   cat(k,substrate,rsd,sep = ":")
      # }
    }
  }
}

#for (cancer in c("BRCA","OV")) {
for (cancer in "BRCA") {
    
  for (iscis in "trans") { # loop around cis and trans
    if (iscis == "cis") {
      var <- "pro_kin"
    }
    if (iscis == "trans") {
      var <- "pho_kin"
    }
    fdr_var <- paste("FDR_",var,sep = "")
    coef_var <- paste("coef_",var,sep = "")
    
    table_multi_temp <- table_multi[table_multi$Cancer==cancer & table_multi$SELF == iscis,]
    table_can_temp <- table_single[table_single$Cancer == cancer & table_single$SELF == iscis,]
    table_multi_temp$scale_FDR <- -log10(table_multi_temp$FDR_pho_kin)
    table_can_temp$scale_FDR <- -log10(table_can_temp$FDR_pho_kin)
    
    merge_id <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","pair")
    merge_cols <- c(merge_id,"scale_FDR","FDR_pho_kin","coef_pho_kin","model","Size")
    table_overlap <- merge(table_multi_temp[,merge_cols],table_can_temp[,merge_cols], by = merge_id)
    
    cat(paste("multi-model overlap with single-model ",cancer," ",iscis," regression result:\n",
              nrow(table_overlap)," kinase:substrate:phosphosite pairs (",
              nrow(table_multi_temp)," PDX pairs vs ",
              nrow(table_can_temp)," BRCA pairs)\n\n",
              sep = ""))
    
    p = ggplot(table_overlap,aes(x=scale_FDR.x, y=scale_FDR.y))
    p = p + geom_point(alpha=0.05)
    #p = p + geom_text(aes(label= ifelse( ( scale_FDR.x >= -log10(0.01) & scale_FDR.y >= -log10(0.01)) , pair, NA)),size=2,alpha=0.5)
    p = p + geom_text(aes(label= ifelse( ( scale_FDR.y >= -log10(0.05) & scale_FDR.x < -log10(0.05) ) , pair, NA)),size=2,alpha=0.2)
    
    #p = p + geom_text(aes(label= ifelse( top & scale_FDR.x >= -log10(0.05) , pair, NA)),size=2,alpha=0.5)
    p = p + theme_bw() #+ theme_nogrid()
    p = p + geom_abline(slope=1)
    p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
    p = p + labs(x = "-log(FDR) in multi-model", y="-log(FDR) in single-model")
    p = p + xlim(0,20) + ylim(0,20)
    p
    
    
  }
}

# substrate <- "APC";rsd <- "S2129"
# substrate <- "RAF1";rsd <- "S29"
substrate <- "NPM1";rsd <- "S70"

old <- table_single[table_single$SUBSTRATE==substrate & table_single$SUB_MOD_RSD==rsd & !table_single$self,]
new <- table_multi[table_multi$SUBSTRATE==substrate & table_multi$SUB_MOD_RSD==rsd,]

change10 <- table_overlap[table_overlap$FDR_pho_kin.y<sig & table_overlap$FDR_pho_kin.x>=sig,]
phosphosites10 <- unique(change10[,c("SUBSTRATE","SUB_MOD_RSD")])
substrates <- as.vector(phosphosites10$SUBSTRATE)
rsds <- as.vector(phosphosites10$SUB_MOD_RSD)
for (i in 1:nrow(phosphosites10)) {
  substrate <- substrates[i]
  rsd <- rsds[i]
  #table_old <- table_single[table_single$SUBSTRATE==substrate & table_single$SUB_MOD_RSD==rsd & !table_single$self,]
  table_new <- table_multi[table_multi$SUBSTRATE==substrate & table_multi$SUB_MOD_RSD==rsd & table_multi$Cancer==cancer,]
  kin_sig <- table_new$KINASE[table_new$FDR_pho_kin<sig]
  if (length(kin_sig) > 0) {
    # print(phosphosites10[i,])
    for (k in kin_sig) {
      k_subs <- unique(table_single$SUBSTRATE[table_single$KINASE==k & !table_single$self])
      print(k_subs)
      # if (length(k_subs)==1) {
      #   cat(k,substrate,rsd,sep = ":")
      # }
    }
  }
}