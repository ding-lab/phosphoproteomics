# Yige Wu @ WashU 2016 Dec
# draw volcano plots for correlation between kinase/substrate phosphoryaltion level and protein phosphorylation level


# choose kinase or phosphotase, outlier threshold and least sample number-------------------------
protein <- "kinase"
cancer <- "BRCA"
sig <- 0.05

# library -----------------------------------------------------------------
library(reshape)
library(stringr)
library(readr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))
source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# input and melt correlation table -------------------------------------------------
ks_phos_corr <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/",protein,"_",cancer,"_kin:sub-based_phos_score_correlate_protein_phos_score.txt", sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
ks_phos_corr_s <- ks_phos_corr[,c("gene","s_fdr","s_estimate")]; colnames(ks_phos_corr_s) <- c("gene","fdr","cor"); ks_phos_corr_s$score_type <- "substrate-based"
ks_phos_corr_k <- ks_phos_corr[,c("gene","k_fdr","k_estimate")]; colnames(ks_phos_corr_k) <- c("gene","fdr","cor"); ks_phos_corr_k$score_type <- "kinase-based"
ks_phos_corr_table <- rbind(ks_phos_corr_s,ks_phos_corr_k)
ks_phos_corr_table$Cancer <- cancer

# plot --------------------------------------------------------------------
p = ggplot(ks_phos_corr_table,aes(x=cor, y=-log10(fdr)))
p = p + facet_grid(Cancer~score_type,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.1 ,  stroke = 0 )
p = p + geom_text(aes(label= ifelse(-log10(fdr)> 3 , as.character(gene), NA)),size=1.5,alpha=0.5,vjust=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Correlation coefficients", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/',protein,"_",cancer,'4corr_kin:sub_phos_score.pdf',sep ="")
ggsave(file=fn, height=3, width=6, useDingbats=FALSE)
