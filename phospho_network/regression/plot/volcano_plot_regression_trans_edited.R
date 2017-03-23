# Yige Wu @ WashU 2016 Dec
# draw volcano plots for regression result

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

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# function -------------------------------------------------------------
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


# input regression processed data -----------------------------------------
# table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",protein,"_substrate_regression_trans_edited.txt",sep = ""))
table_2can_prev <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited_phosphonetwork.txt",sep = ""))
table_2can = table_2can[!(table_2can$pair %in% table_2can_prev$pair),] # only keep the ones not in previous

table_2can$coef_pho_kin_filtered = remove_outliers(table_2can$coef_pho_kin)
table_2can$coef_pro_kin_filtered = remove_outliers(table_2can$coef_pro_kin)

if ( protein == "kinase") {
  plot_fdr_scale <- 3
}
if ( protein == "phosphotase") {
  plot_fdr_scale <- 2
}

table_cis_outlier_removed_m = table_2can[table_2can$SELF=="cis" & !is.na(table_2can$coef_pro_kin_filtered),]
table_trans_outlier_removed_m = table_2can[table_2can$SELF=="trans" & !is.na(table_2can$coef_pho_kin_filtered),]

# cis volcano plotting module -------------------------------------------------
p = ggplot(table_cis_outlier_removed_m,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin)))
p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.05 ,  stroke = 0 )
p = p + geom_text(aes(label= ifelse(-log10(FDR_pro_kin)>plot_fdr_scale | (Cancer=="OV" & -log10(FDR_pro_kin)>2), as.character(pair), NA)),size=1.5,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase protein expression", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/Cis_',protein,'_substrate_volcano_trans_edited_phosphonetwork.pdf',sep ="")
ggsave(file=fn, height=3, width=6, useDingbats=FALSE)

# trans volcano plotting module -------------------------------------------------
p = ggplot(table_trans_outlier_removed_m,aes(x=coef_pho_kin, y=-log10(FDR_pho_kin)))
p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.05 ,  stroke = 0 )
p = p + geom_text(aes(label= ifelse(-log10(FDR_pho_kin)>plot_fdr_scale | (Cancer=="OV" & -log10(FDR_pho_kin)>2), as.character(pair), NA)),size=1.5,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase phosphorylation level", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/Trans_',protein,'_substrate_volcano_trans_edited_phosphonetwork.pdf',sep ="")
ggsave(file=fn, height=3, width=6, useDingbats=FALSE)


# Size~signficance --------------------------------------------------------
table_PDX <- read_delim("~/Box Sync/pan3can_shared_data/analysis_results/tables/kinase_PDX_substrate_regression_trans_edited.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
table_2can_pdx <- rbind(table_2can,table_PDX)
# table_2can$sig[table_2can$self] <- table_2can$FDR_pro_kin[table_2can$self] < 0.05
# table_2can$sig[!table_2can$self] <- table_2can$FDR_pho_kin[!table_2can$self] < 0.05
table_2can_pdx$sig[table_2can_pdx$self] <- table_2can_pdx$FDR_pro_kin[table_2can_pdx$self] < 0.05
table_2can_pdx$sig[!table_2can_pdx$self] <- table_2can_pdx$FDR_pho_kin[!table_2can_pdx$self] < 0.05

p <- ggplot(table_2can_pdx, aes(x = sig , y = Size))
p = p + facet_grid(SELF~Cancer,scales = "free_y", space = "free_y", drop=T)#, , space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p <- p + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
# p = p + geom_text(aes(label= ifelse( (Size > 60 & !sig ), as.character(pair), NA)),size=1.5,alpha=0.5)
p = p + theme_bw()
p = p + theme_nogrid()
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/Size~sig_violin_FDR_',sig,'.pdf',sep ="")
ggsave(file=fn, height=5, width=6)

