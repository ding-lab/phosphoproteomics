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

# input regression processed data -----------------------------------------
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/manuscript/supplementary_tables/kinase_substrate_regression_trans_edited.txt",sep = ""))

# cis
table_cis = table_2can[table_2can$self,]
table_cis$sig = table_cis$FDR_pro_kin < sig
cis_genes = table_cis$KINASE[table_cis$FDR_pro_kin < sig]
cis_genes_3 = names(table(cis_genes)[table(cis_genes)>2])
for (gene in cis_genes_3){
  thres = 0
  all_gene_sites = table_cis[table_cis$KINASE==gene,]
  all_gene_sites_count = all_gene_sites[all_gene_sites$Size>thres,]
  all_gene_sites_count$pho_pos = as.numeric(sub("[A-Z]([0-9]+).*","\\1",all_gene_sites_count$SUB_MOD_RSD))
  
  p = ggplot(all_gene_sites_count,aes(x=pho_pos, y=coef_pro_kin))# make this the original ethni
  p = p + facet_grid(Cancer~., drop=T, space = "free_y")#, space = "free", scales = "free")
  p = p + geom_point(aes(fill=-log10(FDR_pro_kin), color=ifelse(sig, "black",NA), size = Size/20),pch=21) 
  p = p + scale_fill_gradientn(name= "-log(FDR)", na.value=NA, colours=col_paletteR(100))
  p = p + geom_text(aes(label=SUB_MOD_RSD),size=2)
  p = p + scale_colour_manual(values=c("black",NA))
  p = p + theme_bw() #+ theme_nogrid()
  p = p + geom_hline(yintercept = 0, alpha=0.5)
  p = p + labs(x = paste(gene,"phosphosite position"), y="Regression coefficient")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p
  fn = paste(pd,gene,'site_cis_coef.pdf',sep ="_")
  ggsave(file=fn, useDingbats=FALSE)
}

# trans
table_trans = table_2can[!table_2can$self,]
table_trans$sig = table_trans$FDR_pho_kin < sig
table_trans$KS = paste(table_trans$KINASE,table_trans$SUBSTRATE,sep=":")
trans_pairs = table_trans$KS[table_trans$FDR_pho_kin < sig]
trans_pairs_3 = names(table(trans_pairs)[table(trans_pairs)>2])
# trans with a bunch of significant sites
for (pair in trans_pairs_3){
  thres = 0
  all_gene_sites = table_trans[table_trans$KS==pair,]
  all_gene_sites_count = all_gene_sites[all_gene_sites$Size>thres,]
  all_gene_sites_count$pho_pos = as.numeric(sub("[A-Z]([0-9]+).*","\\1",all_gene_sites_count$SUB_MOD_RSD))
  
  p = ggplot(all_gene_sites_count,aes(x=pho_pos, y=coef_pho_kin))# make this the original ethni
  p = p + facet_grid(Cancer~., drop=T, space = "free_y")#, space = "free", scales = "free")
  p = p + geom_point(aes(fill=-log10(FDR_pho_kin), color=ifelse(sig, "black",NA), size = Size/20),pch=21) 
  p = p + scale_fill_gradientn(name= "-log(FDR)", na.value=NA, colours=col_paletteR(100))
  p = p + geom_text(aes(label=SUB_MOD_RSD),size=2)
  p = p + scale_colour_manual(values=c("black",NA))
  p = p + theme_bw() #+ theme_nogrid()
  p = p + geom_hline(yintercept = 0, alpha=0.5)
  p = p + labs(x = paste(pair,"phosphosite position"), y="Regression coefficient")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p
  fn = paste(pd,pair,'site_trans_coef.pdf',sep ="_")
  ggsave(file=fn, useDingbats=FALSE)
}

# trans with mostly non-significant sites but specific significant sites
trans_pairs = unique(table_trans$KS[table_trans$FDR_pho_kin < strict_sig])
thres=40
trans_pairs_count = table(table_trans$KS[all_gene_sites$Size>thres])
trans_pairs_3_001 = names(trans_pairs_count)[(trans_pairs_count > 2) & (names(trans_pairs_count) %in% trans_pairs)]
for (pair in trans_pairs_3_001){
  thres = 0
  all_gene_sites = table_trans[table_trans$KS==pair,]
  all_gene_sites_count = all_gene_sites[all_gene_sites$Size>thres,]
  all_gene_sites_count$pho_pos = as.numeric(sub("[A-Z]([0-9]+).*","\\1",all_gene_sites_count$SUB_MOD_RSD))
  
  p = ggplot(all_gene_sites_count,aes(x=pho_pos, y=coef_pho_kin))# make this the original ethni
  p = p + facet_grid(Cancer~., drop=T, space = "free_y")#, space = "free", scales = "free")
  p = p + geom_point(aes(fill=-log10(FDR_pho_kin), color=ifelse(sig, "black",NA), size = Size/20),pch=21) 
  p = p + scale_fill_gradientn(name= "-log(FDR)", na.value=NA, colours=col_paletteR(100))
  p = p + geom_text(aes(label=SUB_MOD_RSD),size=2)
  p = p + scale_colour_manual(values=c("black",NA))
  p = p + theme_bw() #+ theme_nogrid()
  p = p + geom_hline(yintercept = 0, alpha=0.5)
  p = p + labs(x = paste(pair,"phosphosite position"), y="Regression coefficient")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p
  fn = paste(pd,pair,'site_trans_coef.pdf',sep ="_")
  ggsave(file=fn, useDingbats=FALSE)
}