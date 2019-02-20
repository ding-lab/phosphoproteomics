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

networkinf = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/NetworKIN/networkin_human_predictions_3.1.tsv.gz"
networkin = read.table(header=F, quote = "", sep="\t",file=gzfile(networkinf))
colnames(networkin) = c("Substrate","Position","Name","Score1","KIN","KIN_fam","Score2","ENSP","Score3","Substrate","Sequence","String")
networkin$Group = gsub("_group","",networkin$KIN_fam)

manning_kinome_lookup = manning_kinome[,c("Name","hgnc_symbol")]
networkin_wname = merge(networkin,manning_kinome_lookup,by="Name")

# function -------------------------------------------------------------
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

run_fisher = function(universe, pathway_genes, sig_genes){
  p = NA; OR = NA
  
  fisher_elements = c(sum(!universe %in% c(pathway_genes,sig_genes)),sum((universe %in% pathway_genes) & !(universe %in% sig_genes)),
                      sum(!(universe %in% pathway_genes) & (universe %in% sig_genes)),sum((universe %in% pathway_genes) & (universe %in% sig_genes)))
  if (fisher_elements[1] > 0 && fisher_elements[3] > 0 && fisher_elements[2] >= 0 && fisher_elements[4] >= 0){
    test.table = matrix(as.numeric(fisher_elements), nrow=2)
    f.test = fisher.test(test.table,alternative="greater")
    OR = f.test$estimate
    p = f.test$p.value
  }
  
  return(p)
  #return(list("p"=p, "OR"=OR))
}

if ( protein == "kinase") {
  plot_fdr_scale <- 3
}
if ( protein == "phosphotase") {
  plot_fdr_scale <- 2
}

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

table_HUMAN_cis$PresentInPhosphositeplus = paste(table_HUMAN_cis$KINASE,table_HUMAN_cis$SUBSTRATE,table_HUMAN_cis$SUB_MOD_RSD) %in% paste(k_s_table_phosphosite$GENE,k_s_table_phosphosite$SUB_GENE,k_s_table_phosphosite$SUB_MOD_RSD)
table_HUMAN_trans$PresentInPhosphositeplus = paste(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE,table_HUMAN_trans$SUB_MOD_RSD) %in% paste(k_s_table_phosphosite$GENE,k_s_table_phosphosite$SUB_GENE,k_s_table_phosphosite$SUB_MOD_RSD)
table_HUMAN_cis$PresentInPhosphonetwork = paste(table_HUMAN_cis$KINASE,table_HUMAN_cis$SUBSTRATE,table_HUMAN_cis$SUB_MOD_RSD) %in% paste(k_s_table_network_site$kinase,k_s_table_network_site$substrate,k_s_table_network_site$site)
table_HUMAN_trans$PresentInPhosphonetwork = paste(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE,table_HUMAN_trans$SUB_MOD_RSD) %in% paste(k_s_table_network_site$kinase,k_s_table_network_site$substrate,k_s_table_network_site$site)
table_HUMAN_cis$PresentInExpDatabase = table_HUMAN_cis$PresentInPhosphositeplus | table_HUMAN_cis$PresentInPhosphonetwork
table_HUMAN_trans$PresentInExpDatabase = table_HUMAN_trans$PresentInPhosphositeplus | table_HUMAN_trans$PresentInPhosphonetwork

##table_HUMAN_cis$PresentInNetworkin = paste(table_HUMAN_cis$KINASE,table_HUMAN_cis$SUBSTRATE,table_HUMAN_cis$SUB_MOD_RSD) %in% paste(k_s_table_phosphosite$GENE,k_s_table_phosphosite$SUB_GENE,k_s_table_phosphosite$SUB_MOD_RSD)

# save results for supplementary tables
table_HUMAN_cis_sig = table_HUMAN_cis[table_HUMAN_cis$FDR_pro_kin < sig & table_HUMAN_cis$coef_pro_kin > 0,]
table_HUMAN_trans_sig = table_HUMAN_trans[table_HUMAN_trans$FDR_pho_kin < sig & table_HUMAN_trans$coef_pho_kin > 0,]

# specific sites found in databases
cat("Number of significant KS sites found previously:\n")
fisher.test(table(table_HUMAN_cis$FDR_pro_kin < sig & table_HUMAN_cis$coef_pro_kin > 0,table_HUMAN_cis$PresentInExpDatabase))
CisSigVsExpData = data.frame(table(table_HUMAN_cis$FDR_pro_kin < sig & table_HUMAN_cis$coef_pro_kin > 0,table_HUMAN_cis$PresentInExpDatabase))
colnames(SigVsPhoPlus) = c("Sig","InExpDatabase","Count")

table(table_HUMAN_cis_sig$PresentInPhosphositeplus)
length(unique(table_HUMAN_cis_sig$KINASE[!table_HUMAN_cis_sig$PresentInPhosphonetwork]))

table(table_HUMAN_trans_sig$PresentInPhosphositeplus,table_HUMAN_trans_sig$PresentInPhosphonetwork)

fisher.test(table(table_HUMAN_trans$FDR_pho_kin < sig & table_HUMAN_trans$coef_pho_kin > 0,table_HUMAN_trans$PresentInExpDatabase))
TransSigVsExpData = data.frame(table(table_HUMAN_trans$FDR_pho_kin < sig & table_HUMAN_trans$coef_pho_kin > 0,table_HUMAN_trans$PresentInExpDatabase))
colnames(TransSigVsExpData) = c("Sig","InExpDatabase","Count")

table(table_HUMAN_trans_sig$dataset,table_HUMAN_trans_sig$PresentInPhosphositeplus)
sum(table(table_HUMAN_trans_sig$KINASE[!table_HUMAN_trans_sig$PresentInPhosphositeplus])[c(7,8,19:27,72:80)])
length(unique(table_HUMAN_trans_sig$SUBSTRATE[!table_HUMAN_trans_sig$PresentInPhosphositeplus]))
cat("\n")

num_cis = dim(table_HUMAN_cis_sig[!duplicated(table_HUMAN_cis_sig$KINASE),])[1]
cat("Number of unique tested kinases with significant cis regulation:", num_cis, "\n")
num_trans = length(table_HUMAN_trans_sig$KINASE[!duplicated(table_HUMAN_trans_sig$KINASE)])
table(table_HUMAN_cis_sig$KINASE)[order(table(table_HUMAN_cis_sig$KINASE),decreasing = T)][1:30]

cat("Number of unique tested kinases with significant trans regulation:", num_trans, "\n")
table(table_HUMAN_trans_sig$KINASE)[order(table(table_HUMAN_trans_sig$KINASE),decreasing = T)][1:30]



##### compare in vitro vs. in vivo #####
# cis
table_HUMAN_cis = merge(table_HUMAN_cis,k_s_table_RXN,by=c("KINASE","SUBSTRATE"),all.x=T,all.y=F)
table_HUMAN_cis$evidence="none"
table_HUMAN_cis$evidence[table_HUMAN_cis$IN_VITRO_RXN=="X"]="in vitro"
table_HUMAN_cis$evidence[table_HUMAN_cis$IN_VIVO_RXN=="X"] = "in vivo"
table_HUMAN_cis$evidence[table_HUMAN_cis$IN_VIVO_RXN=="X" & table_HUMAN_cis$IN_VITRO_RXN=="X"] = "both"

cat("Cis evidence vs. validaton","\n")
table(table_HUMAN_cis[table_HUMAN_cis$SELF=="cis",]$FDR_pro_kin<0.05 & table_HUMAN_cis[table_HUMAN_cis$SELF=="cis",]$coef_pro_kin>0,table_HUMAN_cis[table_HUMAN_cis$SELF=="cis",]$evidence)
           
p = ggplot(table_HUMAN_cis,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin), color=evidence))
#p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.2,  stroke = 0 )
#p = p + geom_text(aes(label= ifelse(-log10(FDR_pro_kin)>5 | (Cancer=="OV" & -log10(FDR_pro_kin)>2), as.character(pair), NA)),size=1.5,alpha=0.5)
p = p + theme_bw()
p = p + geom_vline(xintercept = 0, alpha=0.3) + xlim(-2,2)
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase protein expression", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/Cis_kinase_substrate_volcano_in_vivo_vitro.pdf',sep ="")
ggsave(file=fn, height=4, width=4.5, useDingbats=FALSE)

# trans
table_HUMAN_trans = merge(table_HUMAN_trans,k_s_table_RXN,by=c("KINASE","SUBSTRATE"),all.x=T,all.y=F)
table_HUMAN_trans$evidence="none"
table_HUMAN_trans$evidence[table_HUMAN_trans$dataset=="PhosphoNetwork"] = "PhosphoNetwork"
table_HUMAN_trans$evidence[table_HUMAN_trans$IN_VITRO_RXN=="X"]="in vitro"
table_HUMAN_trans$evidence[table_HUMAN_trans$IN_VIVO_RXN=="X"] = "in vivo"
table_HUMAN_trans$evidence[table_HUMAN_trans$IN_VIVO_RXN=="X" & table_HUMAN_trans$IN_VITRO_RXN=="X"] = "both"

cat("Trans evidence vs. validaton","\n")
table(table_HUMAN_trans[table_HUMAN_trans$SELF=="trans",]$FDR_pho_kin<0.05 & table_HUMAN_trans[table_HUMAN_trans$SELF=="trans",]$coef_pho_kin>0,table_HUMAN_trans[table_HUMAN_trans$SELF=="trans",]$evidence)

p = ggplot(table_HUMAN_trans,aes(x=coef_pho_kin, y=-log10(FDR_pho_kin), color =evidence))
p = p + geom_point(alpha=0.05 ,  stroke = 0 , size=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + geom_vline(xintercept = 0, alpha=0.3) + xlim(-2,2)
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase phosphorylation level", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/trans_kinase_substrate_volcano_in_vivo_vitro.pdf',sep ="")
ggsave(file=fn, height=4, width=4.5, useDingbats=FALSE)