# Yige Wu @ WashU 2016 Dec
# draw volcano plots for regression result

# choose kinase/phosphotase, cancer , significance level, model -----------------------------------------------
protein <- "kinase"
sig <- 0.05
cancer <- "BRCA"
mode <- "pho_sub~pro_kin"

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
#table <- table_2can[table_2can$Cancer==c & table_2can$model==mode & table_2can$self==cis,]
table_2can$coef_pro_kin_filtered = remove_outliers(table_2can$coef_pro_kin)
table_2can_outlier_removed_m = table_2can[!is.na(table_2can$coef_pro_kin_filtered) & table_2can$model == mode,]

# function -------------------------------------------------------------
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


# volcano plotting module -------------------------------------------------
if ( protein == "kinase") {
  plot_fdr_scale <- 5
}
if ( protein == "phosphotase") {
  plot_fdr_scale <- 2
}
p = ggplot(table_2can_outlier_removed_m,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin)))
p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.05)
p = p + geom_text(aes(label= ifelse(-log10(FDR_pro_kin)>plot_fdr_scale, pair, NA)),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/',protein,'_substrate_volcano.pdf',sep ="")
ggsave(file=fn, height=6, width=8, useDingbats=FALSE)

