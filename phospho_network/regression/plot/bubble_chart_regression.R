# Yige Wu @ WashU 2016 Dec
# draw bubble charts for regression result

# choose kinase/phosphotase, cancer , significance level, model -----------------------------------------------
protein <- "kinase"
sig <- 0.05
mode <- "pho_sub~pro_kin"
top <- 50 # choose top n rows for FDR_pro_kin for model1
c <- "BRCA" #choose in which cancer extract the top n rows

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

# input regression processed data -----------------------------------------
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",protein,"_substrate_regression.txt",sep = ""))
table_sig <- table_2can[table_2can$model==mode & table_2can$FDR_pro_kin<=sig,]

# bubble chart ---------------------
for (cis in c(TRUE,FALSE)) { # loop around cis and trans
  t0 <- table_sig[table_sig$self==cis,]
  # sort by FDR_pro_kin
  t1 <- t0[t0$Cancer==c,]
  t1 <- t1[order(t1$FDR_pro_kin),] 
  
  ## corresponding results for other two models are extracted and ordered
  rows <- c()
  for(i in 1:top){
    r <- unlist(which(t0$pair==t1$pair[i]))
    rows <- c(rows,r)
  }
  table <- t0[rows,]
  #table_2can_brca_top = table_2can[table_2can$pair %in% table$pair,]
  
  ## actual plotting
  lim = max(max(table$coef_pro_kin),min(table$coef_pro_kin))
  p = ggplot(table,aes(x=model, y=pair))# make this the original ethni
  p = p + facet_grid(KINASE~Cancer, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
  p = p + geom_point(aes(fill=coef_pro_kin, size =-log10(FDR_pro_kin), color=ifelse(sig, "black",NA)),pch=21) 
  p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
  p = p + scale_colour_manual(values=c("black",NA))
  p = p + theme_bw() #+ theme_nogrid()
  p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p
  fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/',c,'_cis_',cis,'_fdr',sig,'_top_',top,'.pdf',sep ="")
  ggsave(file=fn, height=10, width=5, useDingbats=FALSE)
}

