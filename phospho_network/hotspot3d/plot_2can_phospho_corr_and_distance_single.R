# Yige Wu @ WashU 2017 Jan
# plot 3D/linear distance and co-phosphorylation correlation FDRs and coefficients

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


# choose cohort and significance level ------------------------------------
sig <- 0.05

# input within protein pairwise processed file ----------------------------
pairwise_brca <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/BRCA_phosphosite_within_protein_distance_and_correlation.txt", sep=""),"\t", escape_double = FALSE, trim_ws = TRUE)
pairwise_ov <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/OV_phosphosite_within_protein_distance_and_correlation.txt", sep=""),"\t", escape_double = FALSE, trim_ws = TRUE)
pairwise_brca$cancer <- "BRCA"
pairwise_ov$cancer <- "OV"
pairwise <- rbind(pairwise_brca, pairwise_ov)

# plot correlation between coef_corr and distances ------------------------
p = ggplot(data = pairwise, aes(x = dis_3d , y = coef_corr, color = fdr_corr < 0.05))
p = p + facet_grid(.~cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x)
p = p + geom_point(alpha=0.3, stroke=0) #+ scale_color_gradientn(name= "FDR", na.value=NA)
p = p + geom_text(aes(label= ifelse((dis_3d < 4 & coef_corr <0.4) | (dis_3d > 9 & coef_corr > 0.95), pair, NA ), vjust = -1, hjust = 1),size=2,alpha=0.5)
p = p + theme_bw()
p = p + labs(x="3D distance (ångström)", y = "Correlation coefficient")
#p = p + expand_limits(x = 0) 
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/hotspot3d/within_protein_2can_distance_3d_and_coef_corr_correlation.pdf',sep ="")
ggsave(file=fn, height=5, width=10)



limx <- 100
if ( cohort == "OV" ) {
  limx <- 75
}
p = ggplot(data = pairwise[pairwise$dis_lin<limx,], aes(x = dis_lin , y = coef_corr, color = fdr_corr < 0.05))
p = p + facet_grid(.~cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x)
p = p + geom_point(alpha=0.3, stroke = 0)
p = p + geom_text(aes(label= ifelse((dis_lin < 5 & coef_corr < 0.1) | (dis_lin > 0.5*limx & coef_corr > 0.90), pair, NA ), vjust = 1, hjust = -0.2 ),size=2,alpha=0.5)
p = p + theme_bw()
p = p + labs(x="linear distance", y = "Correlation coefficient")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/hotspot3d/within_protein_2can_distance_linear_and_coef_corr_correlation_max',limx,'.pdf',sep ="")
ggsave(file=fn, height=5, width=10)
