# Kuan @ WashU March 2018
# compare quantitatively regulated sites vs. predicted kinase-binding sites
# compare cis in RTK/non RTK
# compare in vivo/in vitro in PhosphositePlus

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
library(ggrepel)

# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"

setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory
# load("/Users/khuang/bin/2015-08-01_Gene_Set.RData")
# KEGG_signaling = KEGG[c(grep("signaling", names(KEGG)),grep("Cell cycle", names(KEGG)))]

# manning kinome allow family-wise analysis
manning_kinome = read.table(header=TRUE, quote = "", sep="\t",file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Gene_family/knome_render_peerj-01-126-s001_wgene.txt")
#manning_kinome = manning_kinome[,-c(12:14)]
manning_kinome_map = manning_kinome[,c(12,9,10)]
colnames(manning_kinome_map)[1] = "KINASE"
manning_kinome_map = manning_kinome_map[!duplicated(manning_kinome_map$KINASE),]

# phosphosite plus data
if ( protein == "kinase" ) {
  ### read in the kinase/substrate table/ phosphorylation data ###
  k_s_table_phosphosite = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
  k_s_table_network = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphoNetworks/comKSI.csv",sep=""))
  k_s_table_network_site = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphoNetworks/highResolutionNetwork.table.csv",sep=""))
}
k_s_table_RXN = k_s_table_phosphosite[,c("GENE","SUB_GENE","IN_VIVO_RXN","IN_VITRO_RXN")]
colnames(k_s_table_RXN) = c("KINASE","SUBSTRATE","IN_VIVO_RXN","IN_VITRO_RXN")
k_s_table_RXN = k_s_table_RXN[!duplicated(paste(k_s_table_RXN$KINASE,k_s_table_RXN$SUBSTRATE)),] # should do an overlap instead, but mostly consistent

networkinf = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/NetworKIN/networkin_human_predictions_3.1.tsv_hugoified.gz"
networkin = read.table(header=F, quote = "", sep="\t",file=gzfile(networkinf))
colnames(networkin) = c("SubstrateString","Position","Name","Score1","KIN","KIN_fam","Score2","ENSP","Score3","Substrate","Sequence","String")
networkin$Group = gsub("_group","",networkin$KIN_fam)

# NOT A PERFECT SOLUTION BUT DOES MOST OF IT
manning_kinome_lookup = manning_kinome[,c("Name","hgnc_symbol")]
networkin_wname = merge(networkin,manning_kinome_lookup,by="Name",all.x=T)
networkin_wname$Name = toupper(networkin_wname$Name)
networkin_wname$Name = gsub("ALPHA","A",networkin_wname$Name)
networkin_wname$Name = gsub("BETA","B",networkin_wname$Name)
networkin_wname$Name = gsub("DELTA","D",networkin_wname$Name)
networkin_wname$Name = gsub("GAMMA","G",networkin_wname$Name)
networkin_wname$Name = gsub("EPSILON","E",networkin_wname$Name)
networkin_wname$Name = gsub("IOTA","I",networkin_wname$Name)
networkin_wname$Name = gsub("THETA","Q",networkin_wname$Name)
networkin_wname$Name = gsub("PKC","PRKC",networkin_wname$Name)
networkin_wname$Name = gsub("AURORA","AURK",networkin_wname$Name)
networkin_wname$Name = gsub("PDHK","PDK",networkin_wname$Name)
networkin_wname$hgnc_symbol[is.na(networkin_wname$hgnc_symbol)] = networkin_wname$Name
# 
# table(networkin_wname$Name[is.na(networkin_wname$hgnc_symbol)])[table(networkin_wname$Name[is.na(networkin_wname$hgnc_symbol)])>3]
# table(networkin_wname$Name[!(networkin_wname$Name %in% manning_kinome_lookup$hgnc_symbol)])
networkin_wname$SUB_MOD_RSD = paste(toupper(gsub("[A-Z]","",networkin_wname$Sequence)),networkin_wname$Position,sep="")
networkin_wname = networkin_wname[!is.na(networkin_wname$hgnc_symbol) & !is.na(networkin_wname$Substrate) & !is.na(networkin_wname$SUB_MOD_RSD),]
networkin_wname$pair = paste(networkin_wname$hgnc_symbol,networkin_wname$Substrate,networkin_wname$SUB_MOD_RSD,sep=":")
networkin_wname_KS = networkin_wname[networkin_wname$hgnc_symbol %in% c(k_s_table_network$Kinase,k_s_table_phosphosite$GENE),]

if ( protein == "kinase") {
  plot_fdr_scale <- 3
}
if ( protein == "phosphotase") {
  plot_fdr_scale <- 2
}


# input regression processed data -----------------------------------------
table_HUMAN_cis = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_cis.txt",sep = ""))
table_HUMAN_trans = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_trans.txt",sep = ""))

table_HUMAN_trans$dataset = "PhosphoNetwork"
table_HUMAN_trans$dataset[ paste(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE) %in% paste(k_s_table_phosphosite$GENE,k_s_table_phosphosite$SUB_GENE)] = "PhosphositePlus"
cat("Number of unique kinase-substrate pairs in trans regulation: \n")
table(table_HUMAN_trans$dataset[!duplicated(paste(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE))])
num_cis = dim(table_HUMAN_cis[!duplicated(table_HUMAN_cis$KINASE),])[1]
cat("Number of unique tested kinases in cis regulation:", num_cis, "\n")
num_trans = length(table_HUMAN_trans$KINASE[!duplicated(table_HUMAN_cis$KINASE)])
cat("Number of unique tested kinases in trans regulation:", num_trans, "\n")

table_HUMAN_cis = merge(table_HUMAN_cis,manning_kinome_map,by="KINASE",all.x=T,all.y=F)
table_HUMAN_trans = merge(table_HUMAN_trans,manning_kinome_map,by="KINASE",all.x=T)
table_HUMAN_cis_sig = table_HUMAN_cis[table_HUMAN_cis$FDR_pro_kin < sig & table_HUMAN_cis$coef_pro_kin > 0,]
table_HUMAN_trans_sig = table_HUMAN_trans[table_HUMAN_trans$FDR_pho_kin < sig & table_HUMAN_trans$coef_pho_kin > 0,]

networkin_wname_cis = networkin_wname_KS[as.character(networkin_wname_KS$hgnc_symbol)==as.character(networkin_wname_KS$Substrate),]
networkin_wname_trans = networkin_wname_KS[as.character(networkin_wname_KS$hgnc_symbol)!=as.character(networkin_wname_KS$Substrate),]

k_s_table_phosphosite_cis = k_s_table_phosphosite[k_s_table_phosphosite$GENE==k_s_table_phosphosite$SUB_GENE,]
k_s_table_network_site_cis = k_s_table_network_site[k_s_table_network_site$kinase==k_s_table_network_site$substrate,]
k_s_table_phosphosite_trans = k_s_table_phosphosite[k_s_table_phosphosite$GENE!=k_s_table_phosphosite$SUB_GENE,]
k_s_table_network_site_trans = k_s_table_network_site[k_s_table_network_site$kinase!=k_s_table_network_site$substrate,]

exp_cis = c(paste(sep=":", k_s_table_phosphosite_cis$GENE,k_s_table_phosphosite_cis$SUB_GENE,k_s_table_phosphosite_cis$SUB_MOD_RSD),
            paste(sep=":", k_s_table_network_site_cis$kinase,k_s_table_network_site_cis$substrate,k_s_table_network_site_cis$site))

exp_trans = c(paste(sep=":", k_s_table_phosphosite_trans$GENE,k_s_table_phosphosite_trans$SUB_GENE,k_s_table_phosphosite_trans$SUB_MOD_RSD),
            paste(sep=":", k_s_table_network_site_trans$kinase,k_s_table_network_site_trans$substrate,k_s_table_network_site_trans$site))


paste(k_s_table_phosphosite$GENE,k_s_table_phosphosite$SUB_GENE,k_s_table_phosphosite$SUB_MOD_RSD)
paste(k_s_table_network_site$kinase,k_s_table_network_site$substrate,k_s_table_network_site$site)

##### Venn diagrams #####
library(VennDiagram)
compare_cis = list(experimental = exp_cis, motif = networkin_wname_cis$pair, quantitative = table_HUMAN_cis_sig$pair)
venn.plot = venn.diagram(
  x = compare_cis, 
  filename = NULL,
  lwd = 4,
  alpha = 0.2,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  margin=0.2
)
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/Cis_kinase_substrate_exp_vs_motif_vs_quant.pdf',sep ="")
pdf(file=fn)
grid.draw(venn.plot)
dev.off()

compare_trans = list(experimental = exp_trans, motif = networkin_wname_trans$pair, quantitative = table_HUMAN_trans_sig$pair)
venn.plot = venn.diagram(
  x = compare_trans, 
  filename = NULL,
  lwd = 4,
  alpha = 0.2,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  margin=0.2
)
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/Trans_kinase_substrate_exp_vs_motif_vs_quant.pdf',sep ="")

pdf(file=fn)
grid.draw(venn.plot)
dev.off()

### limit to observed/tested relationship in our dataset ###


### limit to networKIN kinases ###
compare_cis = list(experimental = exp_cis, motif = networkin_wname_cis$pair, quantitative = table_HUMAN_cis_sig$pair)
venn.plot = venn.diagram(
  x = compare_cis, 
  filename = NULL,
  lwd = 4,
  alpha = 0.2,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  margin=0.2
)
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/Cis_kinase_substrate_exp_vs_motif_vs_quant_networKIN.pdf',sep ="")
pdf(file=fn)
grid.draw(venn.plot)
dev.off()

compare_trans = list(experimental = exp_trans, motif = networkin_wname_trans$pair, quantitative = table_HUMAN_trans_sig$pair)
venn.plot = venn.diagram(
  x = compare_trans, 
  filename = NULL,
  lwd = 4,
  alpha = 0.2,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  margin=0.2
)
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/Trans_kinase_substrate_exp_vs_motif_vs_quant_networKIN.pdf',sep ="")

pdf(file=fn)
grid.draw(venn.plot)
dev.off()


#http://eulerr.co/
three_component_venn = function(a,b,C){

  fit = c(length(a), length(b), length(C), 
          length(intersect(a,b)), length(intersect(a,C)), length(intersect(a,C)),
           length(intersect(intersect(a,b),intersect(a,C))))
  return(fit)
}
three_component_venn(c(1:3),c(2:5),c(1:9))
