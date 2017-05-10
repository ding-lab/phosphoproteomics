# find_discordant_sites.R
# Kuan Huang @ WashU 2017 Feb

# choose kinase/phosphotase, cancer , significance level, model -----------------------------------------------
protein <- "kinase"
sig = 0.05
strict_sig = 0.01

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

plot_trans_landscape = function(pair){
  thres = 0
  all_gene_sites = table_HUMAN_trans[table_HUMAN_trans$KS==pair,]
  all_gene_sites_count = all_gene_sites[all_gene_sites$Size>thres,]
  all_gene_sites_count$pho_pos = as.numeric(sub("[A-Z]([0-9]+).*","\\1",all_gene_sites_count$SUB_MOD_RSD))
  
  p = ggplot(all_gene_sites_count,aes(x=pho_pos, y=coef_pho_kin))# make this the original ethni
  #p = p + facet_grid(Cancer~., drop=T, space = "free_y")#, space = "free", scales = "free")
  p = p + geom_point(aes(fill=-log10(FDR_pho_kin), color=ifelse(sig, "black","grey"), size = Size/20),pch=21)
  p = p + scale_fill_gradientn(name= "-log(FDR)", na.value=NA, colours=col_paletteR(100))
  p = p + scale_size_continuous(breaks = c(1,2,4))
  p = p + geom_text(aes(label=SUB_MOD_RSD),size=2)
  p = p + scale_colour_manual(values=c("black","grey"))
  p = p + theme_bw() #+ theme_nogrid()
  p = p + geom_hline(yintercept = 0, alpha=0.5) + expand_limits(x=0)
  p = p + labs(x = paste(pair,"phosphosite position"), y="Regression coefficient")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p
  fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/landscape/PDX_',pair,'_site_trans_coef.pdf',sep ="")
  ggsave(file=fn, height = 4, width = 7, useDingbats=FALSE)
}

plot_cis_landscape = function(gene){
  thres = 0
  all_gene_sites = table_HUMAN_cis[table_HUMAN_cis$KINASE==gene,]
  all_gene_sites_count = all_gene_sites[all_gene_sites$Size>thres,]
  all_gene_sites_count$pho_pos = as.numeric(sub("[A-Z]([0-9]+).*","\\1",all_gene_sites_count$SUB_MOD_RSD))

  p = ggplot(all_gene_sites_count,aes(x=pho_pos, y=coef_pro_kin))# make this the original ethni
  #p = p + facet_grid(Cancer~., drop=T, space = "free_y")#, space = "free", scales = "free")
  p = p + geom_point(aes(fill=-log10(FDR_pro_kin), color=ifelse(sig, "black","grey"), size = Size/20),pch=21, alpha = 0.8)
  p = p + scale_fill_gradientn(name= "-log(FDR)", na.value=NA, colours=col_paletteR(100))
  p = p + geom_text(aes(label=SUB_MOD_RSD),size=2)
  p = p + scale_colour_manual(values=c("black","grey"))
  p = p + scale_size_continuous(breaks = c(1,2,4))
  p = p + theme_bw() #+ theme_nogrid()
  p = p + geom_hline(yintercept = 0, alpha=0.5) + expand_limits(x=0)
  p = p + labs(x = paste(gene,"phosphosite position"), y="Regression coefficient")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p
  fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/landscape/PDX_',gene,'_site_cis_coef.pdf',sep ="")
  ggsave(file=fn, height = 4, width = 7, useDingbats=FALSE)
}

# input regression processed data -----------------------------------------
table_HUMAN_cis = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_PDX_", protein,"_substrate_regression_cis.txt",sep = ""))
table_HUMAN_trans = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_PDX_", protein,"_substrate_regression_trans.txt",sep = ""))

# cis
table_HUMAN_cis$sig = table_HUMAN_cis$FDR_pro_kin < sig
cis_genes = table_HUMAN_cis$KINASE[table_HUMAN_cis$FDR_pro_kin < sig]
cis_genes_3 = names(table(cis_genes)[table(cis_genes)>2])

# # some dry runs
# gene = "RAF1" # weird the color for none-significant ones not showing up now
# for (gene in cis_genes_3){
#   thres = 0
#   all_gene_sites = table_HUMAN_cis[table_HUMAN_cis$KINASE==gene,]
#   all_gene_sites_count = all_gene_sites[all_gene_sites$Size>thres,]
#   all_gene_sites_count$pho_pos = as.numeric(sub("[A-Z]([0-9]+).*","\\1",all_gene_sites_count$SUB_MOD_RSD))
#   
#   p = ggplot(all_gene_sites_count,aes(x=pho_pos, y=coef_pro_kin))# make this the original ethni
#   #p = p + facet_grid(Cancer~., drop=T, space = "free_y")#, space = "free", scales = "free")
#   p = p + geom_point(aes(fill=-log10(FDR_pro_kin), color=ifelse(sig, "black","grey"), size = Size/20),pch=21, alpha = 0.8) 
#   p = p + scale_fill_gradientn(name= "-log(FDR)", na.value=NA, colours=col_paletteR(100))
#   p = p + geom_text(aes(label=SUB_MOD_RSD),size=2)
#   p = p + scale_colour_manual(values=c("black","grey"))
#   p = p + scale_size_continuous(breaks = c(1,2,4))
#   p = p + theme_bw() #+ theme_nogrid()
#   p = p + geom_hline(yintercept = 0, alpha=0.5) + expand_limits(x=0)
#   p = p + labs(x = paste(gene,"phosphosite position"), y="Regression coefficient")
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p
#   fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/landscape/PDX_',gene,'_site_cis_coef.pdf',sep ="")
#   ggsave(file=fn, height = 4, width = 7, useDingbats=FALSE)
# }

# trans
table_HUMAN_trans$sig = table_HUMAN_trans$FDR_pho_kin < sig
table_HUMAN_trans$KS = paste(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE,sep="-")
trans_pairs = table_HUMAN_trans$KS[table_HUMAN_trans$FDR_pho_kin < sig]
trans_pairs_3 = names(table(trans_pairs)[table(trans_pairs)>2])
# # trans with a bunch of significant sites
# for (pair in trans_pairs_3){
#   thres = 0
#   all_gene_sites = table_HUMAN_trans[table_HUMAN_trans$KS==pair,]
#   all_gene_sites_count = all_gene_sites[all_gene_sites$Size>thres,]
#   all_gene_sites_count$pho_pos = as.numeric(sub("[A-Z]([0-9]+).*","\\1",all_gene_sites_count$SUB_MOD_RSD))
# 
#   p = ggplot(all_gene_sites_count,aes(x=pho_pos, y=coef_pho_kin))# make this the original ethni
#   #p = p + facet_grid(Cancer~., drop=T, space = "free_y")#, space = "free", scales = "free")
#   p = p + geom_point(aes(fill=-log10(FDR_pho_kin), color=ifelse(sig, "black","grey"), size = Size/20),pch=21)
#   p = p + scale_fill_gradientn(name= "-log(FDR)", na.value=NA, colours=col_paletteR(100))
#   p = p + scale_size_continuous(breaks = c(1,2,4))
#   p = p + geom_text(aes(label=SUB_MOD_RSD),size=2)
#   p = p + scale_colour_manual(values=c("black","grey"))
#   p = p + theme_bw() #+ theme_nogrid()
#   p = p + geom_hline(yintercept = 0, alpha=0.5) + expand_limits(x=0)
#   p = p + labs(x = paste(pair,"phosphosite position"), y="Regression coefficient")
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p
#   fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/landscape/PDX_',pair,'_site_trans_coef.pdf',sep ="")
#   ggsave(file=fn, height = 4, width = 7, useDingbats=FALSE)
# }
# 
# # trans with mostly non-significant sites but specific significant sites
# trans_pairs = unique(table_HUMAN_trans$KS[table_HUMAN_trans$FDR_pho_kin < strict_sig])
# thres=40
# trans_pairs_count = table(table_HUMAN_trans$KS[all_gene_sites$Size>thres])
# trans_pairs_3_001 = names(trans_pairs_count)[(trans_pairs_count > 2) & (names(trans_pairs_count) %in% trans_pairs)]
# for (pair in trans_pairs_3_001){
#   thres = 0
#   all_gene_sites = table_HUMAN_trans[table_HUMAN_trans$KS==pair,]
#   all_gene_sites_count = all_gene_sites[all_gene_sites$Size>thres,]
#   all_gene_sites_count$pho_pos = as.numeric(sub("[A-Z]([0-9]+).*","\\1",all_gene_sites_count$SUB_MOD_RSD))
#   
#   p = ggplot(all_gene_sites_count,aes(x=pho_pos, y=coef_pho_kin))# make this the original ethni
#   #p = p + facet_grid(Cancer~., drop=T, space = "free_y")#, space = "free", scales = "free")
#   p = p + geom_point(aes(fill=-log10(FDR_pho_kin), color=ifelse(sig, "black","grey"), size = Size/20),pch=21) 
#   p = p + scale_fill_gradientn(name= "-log(FDR)", na.value=NA, colours=col_paletteR(100))
#   p = p + geom_text(aes(label=SUB_MOD_RSD),size=2)
#   p = p + scale_colour_manual(values=c("black","grey"))
#   p = p + theme_bw() #+ theme_nogrid()
#   p = p + geom_hline(yintercept = 0, alpha=0.5) + expand_limits(x=0)
#   p = p + labs(x = paste(pair,"phosphosite position"), y="Regression coefficient")
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p
#   fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/landscape/PDX_',pair,'_site_trans_disc_coef.pdf',sep ="")
#   ggsave(file=fn, height = 4, width = 7, useDingbats=FALSE)
# }
plot_trans_landscape("MAPK3-RAF1")
plot_cis_landscape("RAF1")
plot_cis_landscape("AKT1")
plot_cis_landscape("ERBB2")
# plot_trans_landscape("EGFR-GAB1")
# plot_trans_landscape("ATM-TP53BP1")
# plot_trans_landscape("CDK1-NUP98")
# plot_trans_landscape("MAP3K5-MAP2K6")
# plot_trans_landscape("MAPK3-MAP2K1")
