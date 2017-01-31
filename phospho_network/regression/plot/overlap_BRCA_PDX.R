# Yige Wu @ WashU 2016 Nov
# look at overlap between BRCA data and PDX data, positive coefficient

# library -----------------------------------------------------------------
library(ggplot2)
library(grid)
library(dplyr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# integrate BRCA & FDX processed data -------------------------------------
load(paste(baseD,"pan3can_shared_data/analysis_results/PDX/Rdata/regression.RData", sep = ""))
table_2can$scale_FDR <- -log10(table_2can$FDR_pro_kin)
table_PDX <- table_2can[table_2can$model=="pho_sub~pro_kin",c("KINASE","SUBSTRATE","SUB_MOD_RSD","size","coef_pro_kin","pair","Cancer","SELF","scale_FDR")]
table_PDX$cohort <- "PDX"

load(paste(baseD,"pan3can_shared_data/analysis_results/regression/Rdata/kinase.RData",sep = ""))
table_2can$SELF <- "trans"; table_2can$SELF[table_2can$self] <- "cis"
table_2can$scale_FDR <- -log10(table_2can$FDR_pro_kin)
table_BRCA <- table_2can[table_2can$Cancer=="BRCA" & table_2can$model=="pho_sub~pro_kin", c("KINASE","SUBSTRATE","SUB_MOD_RSD","size","coef_pro_kin","pair","Cancer","SELF","scale_FDR")]
table_BRCA$cohort <- "BRCA"

table_OV <- table_2can[table_2can$Cancer=="OV" & table_2can$model=="pho_sub~pro_kin", c("KINASE","SUBSTRATE","SUB_MOD_RSD","size","coef_pro_kin","pair","Cancer","SELF","scale_FDR")]
table_OV$cohort <- "OV"

pdx_brca <- merge(table_BRCA, table_PDX, by = c("KINASE","SUBSTRATE","SUB_MOD_RSD","pair","SELF"))
pdx_brca_pos <- pdx_brca[pdx_brca$coef_pro_kin.x>0 & pdx_brca$coef_pro_kin.y>0,]

pdx_ov <- merge(table_OV, table_PDX, by = c("KINASE","SUBSTRATE","SUB_MOD_RSD","pair","SELF"))
pdx_ov_pos <- pdx_ov[pdx_ov$coef_pro_kin.x>0 & pdx_ov$coef_pro_kin.y>0,]

overlap_pos <- rbind(pdx_brca_pos, pdx_ov_pos)

# qqplot: -log10(FDR) of BRCA and PDX -------------------------------------
# plot_fdr_scale <- -log10(0.05)
overlap_pos$top <- FALSE
temp <- overlap_pos[order(overlap_pos$scale_FDR.y, decreasing = TRUE),]
top <- 100
for (self in c("cis","trans")) {
  for (cohort in c("BRCA","OV")) {
    fdr_order <- temp[temp$SELF==self & temp$cohort.x==cohort,]
    rows <- overlap_pos$SELF==self & overlap_pos$cohort.x==cohort
    overlap_pos$top[rows] <- overlap_pos$scale_FDR.y[rows] >= fdr_order$scale_FDR.y[min(top,nrow(fdr_order))]
  }
}

p = ggplot(overlap_pos,aes(x=scale_FDR.x, y=scale_FDR.y))
p = p + facet_grid(SELF~cohort.x,scales = "free_y")
p = p + geom_point(alpha=0.05)
p = p + geom_text(aes(label= ifelse( top & scale_FDR.x >= -log10(0.05) , pair, NA)),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "-logFDR_BRCA/OV", y="-logFDR_PDX")
p = p + xlim(0,10) + ylim(0,10)
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/PDX/overlap/BRCA_OV_regression_overlap_PDX_both_coef_positive_top_',top,'PDX.pdf',sep ="")
ggsave(file=fn, height=8, width=8, useDingbats=FALSE)

# bubble chart ---------------------
temp <- rbind(table_BRCA,table_PDX)
table_sig <- temp[temp$scale_FDR <= -log10(sig),]
# choose the parameters for bubble chart
top <- 50 # choose top n rows for FDR_pro_kin for model1
c <- "BRCA" #choose in which cancer extract the top n rows
for (cis in c(TRUE,FALSE)) { # loop around cis and trans
  t0 <- table_sig[table_sig$self==cis,]
  # sort by FDR_pro_kin
  t1 <- t0[t0$cohort==c,]
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
  fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/PDX_',c,'_cis_',cis,'_fdr',sig,'_top_',top,'.pdf',sep ="")
  #ggsave(file=fn, height=10, width=5, useDingbats=FALSE)
}

