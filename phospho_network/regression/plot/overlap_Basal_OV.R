# Yige Wu @ WashU 2017 Jan
# overlap the BRCA Basal regression result and OV result

# directory and library ---------------------------------------------------
# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"

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


# input processed basal and OV data ---------------------------------------
table_Basal <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/Basal_kinase_substrate_regression.txt",sep = ""), 
                          +     "\t", escape_double = FALSE, trim_ws = TRUE)
table_Basal <- table_Basal[table_Basal$model == "pho_sub~pro_kin",c("KINASE","SUBSTRATE","SUB_MOD_RSD","FDR_pro_kin","coef_pro_kin","size","pair","Cancer","SELF")]

table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/kinase_substrate_regression.txt",sep = ""))
table_OV <- table_2can[table_2can$Cancer=="OV" & table_2can$model=="pho_sub~pro_kin",colnames(table_Basal)]

table_2can <- rbind(table_Basal,table_OV)

table1 <- table_2can[table_2can$FDR_pro_kin <= sig,]

# bubble chart, appoint gene list ---------------------
# choose the parameters for bubble chart
# gene_list <- c("PIK3CA", "PTEN", "AKT1", "TP53", "GATA3", "CDH1", "RB1", "MLL3", "MAP3K1", "CDKN1B","TBX3", "RUNX1", "CBFB", "AFF2", "PIK3R1", "PTPN22", "PTPRD", "NF1", "SF3B1", "CCND3","NF1", "BRCA1", "BRCA2", "RB1", "CDK12","RB1", "NF1", "FAT3", "CSMD3", "GABRA6", "CDK12")
# gene_list <- unique(gene_list)
# 
# RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
# gene_list <- as.vector(t(RTK_file))
gene_list <- genes_sub;

for (cis in c("cis","trans")) { # loop around cis and trans
  t0 <- table1[table1$SELF==cis,]
  rows <- c()
  for(i in gene_list){
    r <- unlist(which(t0$SUBSTRATE==i))
    rows <- c(rows,r)
  }
  table <- t0[rows,]
  if (nrow(table) > 0){
    table$sig <- (table$FDR_pro_kin <= sig)
    ## actual plotting
    lim = max(max(table$coef_pro_kin),min(table$coef_pro_kin))
    p = ggplot(table,aes(x=SELF, y=pair))# make this the original ethni
    #  p = p + facet_grid(KINASE~., drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
    p = p + facet_grid(SUBSTRATE~Cancer, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
    p = p + geom_point(aes(fill=coef_pro_kin, size =-log10(FDR_pro_kin), color=ifelse(sig, "black",NA)),pch=21) 
    p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
    p = p + scale_colour_manual(values=c("black",NA))
    p = p + theme_bw() #+ theme_nogrid()
    p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
    p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
    p
    fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/Basal_OV_overlap_',cis,'_sig_',sig,'.pdf',sep ="")
    ggsave(file=fn, height=12, width=5.5, useDingbats=FALSE)
  }
}


# qqplot ------------------------------------------------------------------
overlap <- merge(table_Basal, table_OV, by = c("KINASE","SUBSTRATE","SUB_MOD_RSD","pair","SELF"))
overlap_pos <- overlap[overlap$coef_pro_kin.x > 0 & overlap$coef_pro_kin.y > 0,]
p = ggplot(overlap_pos,aes(x=-log10(FDR_pro_kin.x), y=-log10(FDR_pro_kin.y)))
p = p + facet_grid(.~SELF,scales = "free_y")
p = p + geom_point(alpha=0.05)
p = p + geom_text(aes(label= ifelse( FDR_pro_kin.x <= sig & FDR_pro_kin.y <= sig , pair, NA)),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "-logFDR_Basal", y="-logFDR_OV")
#p = p + xlim(0,10) + ylim(0,10)
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/kinase/overlap/Basal_OV_regression_overlap_coef_positive_both_sig_',sig,'.pdf',sep ="")
ggsave(file=fn, height=8, width=8, useDingbats=FALSE)

p = ggplot(overlap_pos,aes(x=-log10(FDR_pro_kin.x), y=-log10(FDR_pro_kin.y)))
p = p + facet_grid(.~SELF,scales = "free_y")
p = p + geom_point(alpha=0.05)
p = p + geom_text(aes(label= ifelse( FDR_pro_kin.x <= sig & FDR_pro_kin.y > sig , pair, NA)),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "-logFDR_Basal", y="-logFDR_OV")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/kinase/overlap/Basal_OV_regression_overlap_coef_positive_Basal_sig_',sig,'.pdf',sep ="")
ggsave(file=fn, height=8, width=8, useDingbats=FALSE)

p = ggplot(overlap_pos,aes(x=-log10(FDR_pro_kin.x), y=-log10(FDR_pro_kin.y)))
p = p + facet_grid(.~SELF,scales = "free_y")
p = p + geom_point(alpha=0.05)
p = p + geom_text(aes(label= ifelse( FDR_pro_kin.x > sig & FDR_pro_kin.y <= sig , pair, NA)),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "-logFDR_Basal", y="-logFDR_OV")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/kinase/overlap/Basal_OV_regression_overlap_coef_positive_OV_sig_',sig,'.pdf',sep ="")
ggsave(file=fn, height=8, width=8, useDingbats=FALSE)

# bubble chart ---------------------
# choose the parameters for bubble chart
top <- 50 # choose top n rows for FDR_pro_kin for model1
c <- "Basal" #choose in which cancer extract the top n rows
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
  if(nrow(table) >0) {
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
    ggsave(file=fn, height=10, width=6, useDingbats=FALSE)
  }
}


