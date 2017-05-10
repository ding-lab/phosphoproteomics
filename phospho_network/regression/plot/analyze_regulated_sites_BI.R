# Kuan @ WashU March 2017
# analysis on regulated sites
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

# save results for supplementary tables
table_HUMAN_cis_sig = table_HUMAN_cis[table_HUMAN_cis$FDR_pro_kin < sig & table_HUMAN_cis$coef_pro_kin > 0,]
table_HUMAN_trans_sig = table_HUMAN_trans[table_HUMAN_trans$FDR_pho_kin < sig & table_HUMAN_trans$coef_pho_kin > 0,]

num_cis = dim(table_HUMAN_cis_sig[!duplicated(table_HUMAN_cis_sig$KINASE),])[1]
cat("Number of unique tested kinases with significant cis regulation:", num_cis, "\n")
num_trans = length(table_HUMAN_trans_sig$KINASE[!duplicated(table_HUMAN_trans_sig$KINASE)])
table(table_HUMAN_cis_sig$KINASE)[order(table(table_HUMAN_cis_sig$KINASE),decreasing = T)][1:30]

cat("Number of unique tested kinases with significant trans regulation:", num_trans, "\n")
table(table_HUMAN_trans_sig$KINASE)[order(table(table_HUMAN_trans_sig$KINASE),decreasing = T)][1:30]

# output stats 
cat("Number of total cis K-S relations investigated:",nrow(table_HUMAN_cis),"\n")
cat("Number of significant cis K-S regulations observed:",nrow(table_HUMAN_cis_sig),"\n")
cat("Percentage of significant cis K-S regulations observed:",nrow(table_HUMAN_cis_sig)/nrow(table_HUMAN_cis),"\n")
cat("Number of significant cis K-S regulations observed with coef > 1:",nrow(table_HUMAN_cis_sig[table_HUMAN_cis_sig$coef_pro_kin>1,]),"\n\n")

cat("Number of total trans K-S relations investigated:",nrow(table_HUMAN_trans),"\n")
cat("Number of significant trans K-S regulations observed:",nrow(table_HUMAN_trans_sig),"\n")
cat("Percentage of significant trans K-S regulations observed:",nrow(table_HUMAN_trans_sig)/nrow(table_HUMAN_trans),"\n\n")

cat("Number of total trans K-S relations investigated in PhosphositePlus:",nrow(table_HUMAN_trans[table_HUMAN_trans$dataset=="PhosphositePlus",]),"\n")
cat("Number of significant trans K-S regulations observed:",nrow(table_HUMAN_trans_sig[table_HUMAN_trans_sig$dataset=="PhosphositePlus",]),"\n")
cat("Percentage of significant trans K-S regulations observed:",nrow(table_HUMAN_trans_sig[table_HUMAN_trans_sig$dataset=="PhosphositePlus",])/nrow(table_HUMAN_trans[table_HUMAN_trans$dataset=="PhosphositePlus",]),"\n")
cat("Number of total trans K-S relations investigated in PhosphoNetwork:",nrow(table_HUMAN_trans[table_HUMAN_trans$dataset=="PhosphoNetwork",]),"\n")
cat("Number of significant trans K-S regulations observed:",nrow(table_HUMAN_trans_sig[table_HUMAN_trans_sig$dataset=="PhosphoNetwork",]),"\n")
cat("Percentage of significant trans K-S regulations observed:",nrow(table_HUMAN_trans_sig[table_HUMAN_trans_sig$dataset=="PhosphoNetwork",])/nrow(table_HUMAN_trans[table_HUMAN_trans$dataset=="PhosphoNetwork",]),"\n\n")

trans_sig_site_c = length(unique(paste(table_HUMAN_trans_sig$SUBSTRATE,table_HUMAN_trans_sig$SUB_MOD_RSD)))
trans_all_site_c = length(unique(paste(table_HUMAN_trans$SUBSTRATE,table_HUMAN_trans$SUB_MOD_RSD)))

cat("Number of kinases in significantly associated cis pairs:",length(unique(table_HUMAN_cis$KINASE)),"\n")
cat("Number of kinases in significantly associated trans pairs:",length(unique(table_HUMAN_trans$KINASE)),"\n")
cat("Number of all genes in significantly associated trans pairs:",length(unique(c(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE))),"\n")

cat("Number of sig phosphosites in trans analysis:",trans_sig_site_c,"\n")
cat("Number of all phosphosites in trans analysis:",trans_all_site_c,"\n")
cat("Percentage of sig phosphosites in trans analysis:",trans_sig_site_c/trans_all_site_c,"\n\n")

# # write output table
# tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_cis_sig_fam.txt",sep = "")
# write.table(table_HUMAN_cis_sig, file=tn, quote=F, sep = '\t', row.names = FALSE)
# tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_trans_sig_fam.txt",sep = "")
# write.table(table_HUMAN_trans_sig, file=tn, quote=F, sep = '\t', row.names = FALSE)

# ##### scan for pathway enrichment
# 
# stats = matrix(,,ncol=9)
# colnames(stats) = c("hsa","Pathway_name","Num_pathway_genes","Num_genes_in_cis_analysis",
#                     "Num_sig_genes_in_cis_analysis", "Cis_enrichment_P","Num_genes_in_trans_analysis", "Num_sig_genes_in_trans_analysis","Trans_enrichment_P")
# 
# for (pathway in names(KEGG_signaling)){
#   hsa = strsplit(pathway, split = "\t")[[1]][1]
#   pathway_name=strsplit(pathway, split = "\t")[[1]][2]
#   pathway_genes = KEGG[[pathway]]
# 
#   numPathwayGene = length(pathway_genes)
#   inCisSet = sum(pathway_genes %in% table_HUMAN_cis$KINASE)
#   inTransSet = sum(pathway_genes %in% c(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE))
#   sigInCisSet = sum(pathway_genes %in% table_HUMAN_cis_sig$KINASE)
#   sigInTransSet = sum(pathway_genes %in% c(table_HUMAN_trans_sig$KINASE,table_HUMAN_trans_sig$SUBSTRATE))
# 
#   cis_p = run_fisher(universe = unique(table_HUMAN_cis$KINASE),pathway_genes,unique(table_HUMAN_cis_sig$KINASE))
#   trans_p = run_fisher(universe = unique(c(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE)),pathway_genes,unique(c(table_HUMAN_trans_sig$KINASE,table_HUMAN_trans_sig$SUBSTRATE)))
# 
#   a = c(hsa, pathway_name,numPathwayGene,inCisSet, sigInCisSet, cis_p, inTransSet, sigInTransSet, trans_p)
#   stats=rbind(stats,a)
#   row.names(stats)=NULL
# }
# stats = data.frame(stats)
# stats$Cis_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Cis_enrichment_P)), method="BH")
# stats$Trans_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Trans_enrichment_P)), method="BH")
# stats = stats[order(stats$Trans_enrichment_P),]
# 
# # tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/sig_pairs_KEGG_enrichment.txt", sep="")
# # write.table(stats, file=tn, quote=F, sep = '\t', row.names = FALSE)

# kinase family
stats = matrix(,,ncol=8)
colnames(stats) = c("family_name","Num_Family_genes","Num_family_in_cis_analysis",
                    "Num_sig_genes_in_cis_analysis", "Cis_enrichment_P","Num_genes_in_trans_analysis", "Num_sig_genes_in_trans_analysis","Trans_enrichment_P")

for (family in unique(manning_kinome_map$Family)){
  fam_genes = as.character(manning_kinome_map$KINASE[manning_kinome_map$Family==family])
# for (family in unique(manning_kinome_map$GroupName)){
#   fam_genes = as.character(manning_kinome_map$KINASE[manning_kinome_map$GroupName==family])
  numPathwayGene = length(fam_genes)
  inCisSet = sum(fam_genes %in% table_HUMAN_cis$KINASE)
  inTransSet = sum(fam_genes %in% c(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE))
  sigInCisSet = sum(fam_genes %in% table_HUMAN_cis_sig$KINASE)
  sigInTransSet = sum(fam_genes %in% c(table_HUMAN_trans_sig$KINASE,table_HUMAN_trans_sig$SUBSTRATE))

  cis_p = run_fisher(universe = unique(table_HUMAN_cis$KINASE),fam_genes,unique(table_HUMAN_cis_sig$KINASE))
  trans_p = run_fisher(universe = unique(c(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE)),fam_genes,unique(c(table_HUMAN_trans_sig$KINASE,table_HUMAN_trans_sig$SUBSTRATE)))

  a = c(family,numPathwayGene,inCisSet, sigInCisSet, cis_p, inTransSet, sigInTransSet, trans_p)
  stats=rbind(stats,a)
  row.names(stats)=NULL
}
stats = data.frame(stats)
stats$Cis_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Cis_enrichment_P)), method="BH")
stats$Trans_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Trans_enrichment_P)), method="BH")
stats = stats[order(stats$Trans_enrichment_P),]

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/sig_pairs_K_fam_enrichment.txt", sep="")
write.table(stats, file=tn, quote=F, sep = '\t', row.names = FALSE)

# kinase group
stats = matrix(,,ncol=8)
colnames(stats) = c("group_name","Num_group_genes","Num_group_in_cis_analysis",
                    "Num_sig_genes_in_cis_analysis", "Cis_enrichment_P","Num_genes_in_trans_analysis", "Num_sig_genes_in_trans_analysis","Trans_enrichment_P")

for (group in unique(manning_kinome_map$GroupName)){
  group_genes = as.character(manning_kinome_map$KINASE[manning_kinome_map$GroupName==group])
  # for (group in unique(manning_kinome_map$GroupName)){
  #   group_genes = as.character(manning_kinome_map$KINASE[manning_kinome_map$GroupName==group])
  numPathwayGene = length(group_genes)
  inCisSet = sum(group_genes %in% table_HUMAN_cis$KINASE)
  inTransSet = sum(group_genes %in% c(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE))
  sigInCisSet = sum(group_genes %in% table_HUMAN_cis_sig$KINASE)
  sigInTransSet = sum(group_genes %in% c(table_HUMAN_trans_sig$KINASE,table_HUMAN_trans_sig$SUBSTRATE))
  
  cis_p = run_fisher(universe = unique(table_HUMAN_cis$KINASE),group_genes,unique(table_HUMAN_cis_sig$KINASE))
  trans_p = run_fisher(universe = unique(c(table_HUMAN_trans$KINASE,table_HUMAN_trans$SUBSTRATE)),group_genes,unique(c(table_HUMAN_trans_sig$KINASE,table_HUMAN_trans_sig$SUBSTRATE)))
  
  a = c(group,numPathwayGene,inCisSet, sigInCisSet, cis_p, inTransSet, sigInTransSet, trans_p)
  stats=rbind(stats,a)
  row.names(stats)=NULL
}
stats = data.frame(stats)
stats$Cis_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Cis_enrichment_P)), method="BH")
stats$Trans_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Trans_enrichment_P)), method="BH")
stats = stats[order(stats$Trans_enrichment_P),]

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/sig_pairs_K_group_enrichment.txt", sep="")
write.table(stats, file=tn, quote=F, sep = '\t', row.names = FALSE)

##### auto-phosphorylation #####
table_HUMAN_cis_auto = table_HUMAN_cis[table_HUMAN_cis$coef_pro_kin>.9,]
table_HUMAN_cis_auto$FDR_pro_kin_plot = table_HUMAN_cis_auto$FDR_pro_kin
table_HUMAN_cis_auto$FDR_pro_kin_plot[table_HUMAN_cis_auto$FDR_pro_kin_plot < 10^(-10)] = 10^(-10)
p = ggplot(table_HUMAN_cis_auto,aes(y=KINASE, x=coef_pro_kin, color=-log10(FDR_pro_kin_plot)))
#p = p + geom_text_repel(aes(label=SUB_MOD_RSD), size=3)
p = p + geom_text(aes(label=SUB_MOD_RSD), size=3, alpha=0.7)
#p = p + geom_point(stroke = 0 )
p = p + theme_bw()
p = p + geom_vline(xintercept = 1, alpha=0.3) #+ xlim(-2,2)
p = p + theme(axis.title = element_text(size=14), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + theme(legend.position = "bottom")
p = p + labs(x = "Coefficient for kinase protein expression", y="Kinase")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/Cis_kinase_autophosphorylation.pdf',sep ="")
ggsave(file=fn, height=8, useDingbats=FALSE)

# ##### compare RTK vs. non RTK in cis effects (autophosphorylation) #####
# # in Broad run non-RTK showed more tail than RTK
# RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
# rtk = as.vector(t(RTK_file))
# 
# table_HUMAN_cis_sig$RTK = "Other"
# table_HUMAN_cis_sig$RTK[table_HUMAN_cis_sig$KINASE %in% rtk] = "RTK"
# table_HUMAN_cis_auto = table_HUMAN_cis[table_HUMAN_cis$coef_pro_kin>1,]
# 
# p = ggplot(table_HUMAN_cis_sig,aes(x=coef_pro_kin,fill=RTK),color=NA)
# #p = ggplot(table_cis_outlier_removed_m_auto,aes(x=coef_pro_kin,fill=RTK),color=NA)
# p = p + geom_density(alpha=0.8)
# p = p + theme_bw()
# p = p + geom_vline(xintercept =1,alpha=0.5)
# p = p + labs(x="regression coefficient")
# p
# fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/Cis_protein_substrate_density_trans_edited_vs_RTK.pdf',sep ="")
# ggsave(file=fn, height=3, width=6, useDingbats=FALSE)

# p = ggplot(table_cis_outlier_removed_m,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin), color=RTK))
# p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
# p = p + geom_point(alpha=0.2,  stroke = 0 )
# p = p + geom_text(aes(label= ifelse(-log10(FDR_pro_kin)>5 | (Cancer=="OV" & -log10(FDR_pro_kin)>2), as.character(pair), NA)),size=1.5,alpha=0.5)
# p = p + theme_bw() + xlim(1,2)
# p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
# p = p + labs(x = "Coefficient for kinase protein expression", y="-log10(FDR)")
# p
# fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/Cis_protein_substrate_volcano_trans_edited_vs_RTK_effect_gt1.pdf',sep ="")
# ggsave(file=fn, height=3, width=6, useDingbats=FALSE)


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