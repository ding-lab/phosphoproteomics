# Yige Wu @ WashU 2017 Jan
# adopted by Kuan @ WashU 2017 March
# plot 3D/linear distance and co-phosphorylation correlation FDRs and coefficients

# directory and library ---------------------------------------------------
# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"

library(stringr)
library(ggplot2)
library(readr)


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# choose cohort and significance level ------------------------------------
sig <- 0.05

# input within protein pairwise processed file ----------------------------
pairwise_brca <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_phosphosite_within_protein_distance_and_correlation_broad.txt", sep=""),sep="\t",header=T)
brca_f = paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_phosphosite_within_protein_distance_and_correlation_broad.txt", sep="")
pairwise_brca = read.table(header=T,sep="\t", file = brca_f)

cor.test(pairwise_brca$coef_corr,pairwise_brca$dis_lin, method = "spearman")
cor.test(pairwise_brca$coef_corr,pairwise_brca$dis_3d, method = "spearman")

cor.test(pairwise_brca$coef_corr,pairwise_brca$dis_lin, method = "pearson")
cor.test(pairwise_brca$coef_corr,pairwise_brca$dis_3d, method = "pearson")

cor.test(pairwise_brca$coef_corr[pairwise_brca$dis_lin>5],pairwise_brca$dis_lin[pairwise_brca$dis_lin>5], method = "spearman")
cor.test(pairwise_brca$coef_corr[pairwise_brca$dis_lin>5],pairwise_brca$dis_3d[pairwise_brca$dis_lin>5], method = "spearman")
# pairwise_ov <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/OV_phosphosite_within_protein_distance_and_correlation.txt", sep=""),"\t", escape_double = FALSE, trim_ws = TRUE)
# pairwise_brca$cancer <- "BRCA"
# pairwise_ov$cancer <- "OV"
# pairwise <- rbind(pairwise_brca, pairwise_ov)

# # plot correlation between coef_corr and distances ------------------------
# p = ggplot(data = pairwise, aes(x = dis_3d , y = coef_corr, color = fdr_corr < 0.05))
# p = p + facet_grid(.~cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
# p = p + geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x)
# p = p + geom_point(alpha=0.3, stroke=0) #+ scale_color_gradientn(name= "FDR", na.value=NA)
# p = p + geom_text(aes(label= ifelse((dis_3d < 4 & coef_corr <0.4) | (dis_3d > 9 & coef_corr > 0.95), pair, NA ), vjust = -1, hjust = 1),size=2,alpha=0.5)
# p = p + theme_bw()
# p = p + labs(x="3D distance (ångström)", y = "Correlation coefficient")
# #p = p + expand_limits(x = 0) 
# p
# fn = paste(baseD,'pan3can_shared_data/analysis_results/hotspot3d/within_protein_2can_distance_3d_and_coef_corr_correlation.pdf',sep ="")
# ggsave(file=fn, height=5, width=10)

# plot correlation between coef_corr and distances ------------------------
p = ggplot(data = pairwise_brca, aes(x = dis_3d , y = coef_corr))
#p = p + facet_grid(.~cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_smooth(method = "glm", color="black", formula = y ~ x)
p = p + geom_point(alpha=0.3, stroke=0) #+ scale_color_gradientn(name= "FDR", na.value=NA)
p = p + theme_bw()
p = p + labs(x="3D distance (ångström)", y = "Correlation coefficient")
#p = p + expand_limits(x = 0) 
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/phosphopairs_correlation_3d_distance.pdf',sep ="")
ggsave(file=fn, height=5, width=5,useDingbats=F)

# plot correlation between coef_corr and distances ------------------------
p = ggplot(data = pairwise_brca, aes(x = dis_lin , y = coef_corr))
#p = p + facet_grid(.~cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_smooth(method = "glm", color="black", formula = y ~ x)
p = p + geom_point(alpha=0.3, stroke=0) #+ scale_color_gradientn(name= "FDR", na.value=NA)
p = p + theme_bw() + xlim(0,25)
p = p + labs(x="Linear distance (aa)", y = "Correlation coefficient")
#p = p + expand_limits(x = 0) 
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/phosphopairs_correlation_linear_distance.pdf',sep ="")
ggsave(file=fn, height=5, width=5,useDingbats=F)