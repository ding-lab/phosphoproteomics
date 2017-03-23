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

# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory
load("/Users/khuang/bin/2015-08-01_Gene_Set.RData")
KEGG_signaling = KEGG[c(grep("signaling", names(KEGG)),grep("Cell cycle", names(KEGG)))]

# manning kinome allow family-wise analysis
manning_kinome = read.table(header=TRUE, quote = "", sep="\t",file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Gene_family/knome_render_peerj-01-126-s001_wgene.txt")
#manning_kinome = manning_kinome[,-c(12:14)]
manning_kinome_map = manning_kinome[,c(15,9,10)]
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
    f.test = fisher.test(test.table)
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
  k_s_table = read.delim(paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
}
k_s_table_RXN = k_s_table[,c("GENE","SUB_GENE","IN_VIVO_RXN","IN_VITRO_RXN")]
colnames(k_s_table_RXN) = c("KINASE","SUBSTRATE","IN_VIVO_RXN","IN_VITRO_RXN")
k_s_table_RXN = k_s_table_RXN[!duplicated(c(k_s_table_RXN$KINASE,k_s_table$SUBSTRATE)),] # should do an overlap instead, but mostly consistent

# input regression processed data -----------------------------------------
table_2can <- read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))
table_2can = merge(table_2can,manning_kinome_map,by="KINASE",all.x=T,all.y=F)

table_2can_new = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited_phosphonetwork.txt",sep = ""))
table_2can_new = table_2can_new[!(table_2can_new$pair %in% table_2can$pair),] # only keep the ones not in previous
table_2can_new = merge(table_2can_new,manning_kinome_map,by="KINASE",all.x=T)

# save results for supplementary tables
table_2can_cis_sig = table_2can[table_2can$self & table_2can$FDR_pro_kin < sig & table_2can$coef_pro_kin > 0,]
table_2can_trans_sig = table_2can[!table_2can$self & table_2can$FDR_pho_kin < sig & table_2can$coef_pho_kin > 0,]
table_2can_new_trans_sig = table_2can_new[!table_2can_new$self & table_2can_new$FDR_pho_kin < sig & table_2can_new$coef_pho_kin > 0,]

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited_cis_sig.txt", sep="")
write.table(table_2can_cis_sig, file=tn, quote=F, sep = '\t', row.names = FALSE)
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited_trans_sig.txt", sep="")
write.table(table_2can_trans_sig, file=tn, quote=F, sep = '\t', row.names = FALSE)
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited_phosphonetwork_trans_sig.txt",sep = "")
write.table(table_2can_new_trans_sig, file=tn, quote=F, sep = '\t', row.names = FALSE)

##### scan for pathway enrichment
cis_results = rbind(table_2can[table_2can$self,],table_2can_new[table_2can_new$self,])
trans_results = rbind(table_2can[!table_2can$self,],table_2can_new[!table_2can_new$self,])
cis_results_sig = cis_results[cis_results$FDR_pro_kin < sig & cis_results$coef_pro_kin > 0,]
trans_results_sig = trans_results[trans_results$FDR_pho_kin < sig & trans_results$coef_pho_kin > 0,]

trans_sig_site_c = length(unique(paste(trans_results_sig$SUBSTRATE,trans_results_sig$SUB_MOD_RSD)))
trans_all_site_c = length(unique(paste(trans_results$SUBSTRATE,trans_results$SUB_MOD_RSD)))

cat("Number of all genes in significantly associated cis pairs:",length(unique(cis_results$KINASE)),"\n")
cat("Number of all genes in significantly associated trans pairs:",length(unique(c(trans_results$KINASE,trans_results$SUBSTRATE))),"\n")

cat("Number of sig phosphosites in trans analysis:",trans_sig_site_c,"\n")
cat("Number of all phosphosites in trans analysis:",trans_all_site_c,"\n")
cat("Percentage of sig phosphosites in trans analysis:",trans_sig_site_c/trans_all_site_c,"\n")

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
#   inCisSet = sum(pathway_genes %in% cis_results$KINASE)
#   inTransSet = sum(pathway_genes %in% c(trans_results$KINASE,trans_results$SUBSTRATE))
#   sigInCisSet = sum(pathway_genes %in% cis_results_sig$KINASE)
#   sigInTransSet = sum(pathway_genes %in% c(trans_results_sig$KINASE,trans_results_sig$SUBSTRATE))
#   
#   cis_p = run_fisher(universe = unique(cis_results$KINASE),pathway_genes,unique(cis_results_sig$KINASE))
#   trans_p = run_fisher(universe = unique(c(trans_results$KINASE,trans_results$SUBSTRATE)),pathway_genes,unique(c(trans_results_sig$KINASE,trans_results_sig$SUBSTRATE)))
#   
#   a = c(hsa, pathway_name,numPathwayGene,inCisSet, sigInCisSet, cis_p, inTransSet, sigInTransSet, trans_p)
#   stats=rbind(stats,a)
#   row.names(stats)=NULL
# }
# stats = data.frame(stats)
# stats$Cis_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Cis_enrichment_P)), method="BH")
# stats$Trans_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Trans_enrichment_P)), method="BH")
# 
# tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/sig_pairs_KEGG_enrichment.txt", sep="")
# write.table(stats, file=tn, quote=F, sep = '\t', row.names = FALSE)

# kinase family
stats = matrix(,,ncol=8)
colnames(stats) = c("family_name","Num_Family_genes","Num_family_in_cis_analysis",
                    "Num_sig_genes_in_cis_analysis", "Cis_enrichment_P","Num_genes_in_trans_analysis", "Num_sig_genes_in_trans_analysis","Trans_enrichment_P")

for (family in unique(manning_kinome_map$Family)){
  fam_genes = as.character(manning_kinome_map$KINASE[manning_kinome_map$Family==family])
# for (family in unique(manning_kinome_map$GroupName)){
#   fam_genes = as.character(manning_kinome_map$KINASE[manning_kinome_map$GroupName==family])
  numPathwayGene = length(fam_genes)
  inCisSet = sum(fam_genes %in% cis_results$KINASE)
  inTransSet = sum(fam_genes %in% c(trans_results$KINASE,trans_results$SUBSTRATE))
  sigInCisSet = sum(fam_genes %in% cis_results_sig$KINASE)
  sigInTransSet = sum(fam_genes %in% c(trans_results_sig$KINASE,trans_results_sig$SUBSTRATE))

  cis_p = run_fisher(universe = unique(cis_results$KINASE),fam_genes,unique(cis_results_sig$KINASE))
  trans_p = run_fisher(universe = unique(c(trans_results$KINASE,trans_results$SUBSTRATE)),fam_genes,unique(c(trans_results_sig$KINASE,trans_results_sig$SUBSTRATE)))

  a = c(family,numPathwayGene,inCisSet, sigInCisSet, cis_p, inTransSet, sigInTransSet, trans_p)
  stats=rbind(stats,a)
  row.names(stats)=NULL
}
stats = data.frame(stats)
stats$Cis_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Cis_enrichment_P)), method="BH")
stats$Trans_enrichment_FDR = p.adjust(as.numeric(as.character(stats$Trans_enrichment_P)), method="BH")

tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/sig_pairs_K_fam_enrichment.txt", sep="")
write.table(stats, file=tn, quote=F, sep = '\t', row.names = FALSE)

#####

table_2can = merge(table_2can,k_s_table_RXN,by=c("KINASE","SUBSTRATE"),all.x=T,all.y=F)
table_2can$evidence="none"
table_2can$evidence[table_2can$IN_VITRO_RXN=="X"]="in vitro"
table_2can$evidence[table_2can$IN_VIVO_RXN=="X"] = "in vivo"
table_2can$evidence[table_2can$IN_VIVO_RXN=="X" & table_2can$IN_VITRO_RXN=="X"] = "both"

table_2can$coef_pho_kin_filtered = remove_outliers(table_2can$coef_pho_kin)
table_2can$coef_pro_kin_filtered = remove_outliers(table_2can$coef_pro_kin)

table_cis_outlier_removed_m = table_2can[table_2can$SELF=="cis" & !is.na(table_2can$coef_pro_kin_filtered),]
table_trans_outlier_removed_m = table_2can[table_2can$SELF=="trans" & !is.na(table_2can$coef_pho_kin_filtered),]

##### compare RTK vs. non RTK in cis effects (autophosphorylation) #####
RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
rtk = as.vector(t(RTK_file))

table_cis_outlier_removed_m$RTK = "Other"
table_cis_outlier_removed_m$RTK[table_cis_outlier_removed_m$KINASE %in% rtk] = "RTK"
table_cis_outlier_removed_m_auto = table_cis_outlier_removed_m[table_cis_outlier_removed_m$coef_pro_kin>1,]

p = ggplot(table_cis_outlier_removed_m,aes(x=coef_pro_kin,fill=RTK),color=NA)
#p = ggplot(table_cis_outlier_removed_m_auto,aes(x=coef_pro_kin,fill=RTK),color=NA)
p = p + geom_density(alpha=0.8)
p = p + theme_bw()
p = p + geom_vline(xintercept =1,alpha=0.5)
p = p + labs(x="regression coefficient")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/Cis_protein_substrate_density_trans_edited_vs_RTK.pdf',sep ="")
ggsave(file=fn, height=3, width=6, useDingbats=FALSE)

p = ggplot(table_cis_outlier_removed_m,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin), color=RTK))
p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.2,  stroke = 0 )
p = p + geom_text(aes(label= ifelse(-log10(FDR_pro_kin)>5 | (Cancer=="OV" & -log10(FDR_pro_kin)>2), as.character(pair), NA)),size=1.5,alpha=0.5)
p = p + theme_bw() + xlim(1,2)
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase protein expression", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/Cis_protein_substrate_volcano_trans_edited_vs_RTK_effect_gt1.pdf',sep ="")
ggsave(file=fn, height=3, width=6, useDingbats=FALSE)


##### compare in vitro vs. in vivo #####
# cis
cat("Cis evidence vs. validaton","\n")
table(table_2can[table_2can$SELF=="cis",]$FDR_pro_kin<0.05 & table_2can[table_2can$SELF=="cis",]$coef_pro_kin>0,table_2can[table_2can$SELF=="cis",]$evidence)
           
p = ggplot(table_cis_outlier_removed_m,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin), color=evidence))
p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.2,  stroke = 0 )
#p = p + geom_text(aes(label= ifelse(-log10(FDR_pro_kin)>5 | (Cancer=="OV" & -log10(FDR_pro_kin)>2), as.character(pair), NA)),size=1.5,alpha=0.5)
p = p + theme_bw()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase protein expression", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/Cis_protein_substrate_volcano_trans_edited_in_vivo_vitro.pdf',sep ="")
ggsave(file=fn, height=3, width=6, useDingbats=FALSE)

# trans
cat("Trans evidence vs. validaton","\n")
table(table_2can[table_2can$SELF=="trans",]$FDR_pho_kin<0.05 & table_2can[table_2can$SELF=="trans",]$coef_pho_kin>0,table_2can[table_2can$SELF=="trans",]$evidence)

# trans volcano plotting module -------------------------------------------------
p = ggplot(table_trans_outlier_removed_m,aes(x=coef_pho_kin, y=-log10(FDR_pho_kin), color=evidence))
p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.05 ,  stroke = 0 )
#p = p + geom_text(aes(label= ifelse(-log10(FDR_pho_kin)>plot_fdr_scale | (Cancer=="OV" & -log10(FDR_pho_kin)>2), as.character(pair), NA)),size=1.5,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase phosphorylation level", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/Trans_',protein,'_substrate_volcano_trans_edited_in_vivo_vitro.pdf',sep ="")
ggsave(file=fn, height=3, width=6, useDingbats=FALSE)


# # Size~signficance --------------------------------------------------------
# table_PDX <- read_delim("~/Box Sync/pan3can_shared_data/analysis_results/tables/kinase_PDX_substrate_regression_trans_edited.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# table_2can_pdx <- rbind(table_2can,table_PDX)
# # table_2can$sig[table_2can$self] <- table_2can$FDR_pro_kin[table_2can$self] < 0.05
# # table_2can$sig[!table_2can$self] <- table_2can$FDR_pho_kin[!table_2can$self] < 0.05
# table_2can_pdx$sig[table_2can_pdx$self] <- table_2can_pdx$FDR_pro_kin[table_2can_pdx$self] < 0.05
# table_2can_pdx$sig[!table_2can_pdx$self] <- table_2can_pdx$FDR_pho_kin[!table_2can_pdx$self] < 0.05
# 
# p <- ggplot(table_2can_pdx, aes(x = sig , y = Size))
# p = p + facet_grid(SELF~Cancer,scales = "free_y", space = "free_y", drop=T)#, , space = "free_y",scales = "free_y")#, space = "free", scales = "free")
# p <- p + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
# # p = p + geom_text(aes(label= ifelse( (Size > 60 & !sig ), as.character(pair), NA)),size=1.5,alpha=0.5)
# p = p + theme_bw()
# p = p + theme_nogrid()
# p
# fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/Size~sig_violin_FDR_',sig,'.pdf',sep ="")
# ggsave(file=fn, height=5, width=6)
# 
