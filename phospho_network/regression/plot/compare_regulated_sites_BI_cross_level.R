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
  cat("P value:",p,"\n")
  cat("OR:", OR,"\n")
  #return(p)
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

table_HUMAN_cis_rna = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_cis_rna.txt",sep = ""))
table_HUMAN_trans_rna = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_trans_rna.txt",sep = ""))

table_HUMAN_trans_pro = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_trans_pro.txt",sep = ""))

table_HUMAN_cis_sig = table_HUMAN_cis[table_HUMAN_cis$FDR_pro_kin < sig & table_HUMAN_cis$coef_pro_kin > 0,]
table_HUMAN_trans_sig = table_HUMAN_trans[table_HUMAN_trans$FDR_pho_kin < sig & table_HUMAN_trans$coef_pho_kin > 0,]
table_HUMAN_cis_rna_sig = table_HUMAN_cis_rna[table_HUMAN_cis_rna$FDR_pro_kin < sig & table_HUMAN_cis_rna$coef_pro_kin > 0,]
table_HUMAN_trans_rna_sig = table_HUMAN_trans_rna[table_HUMAN_trans_rna$FDR_pho_kin < sig & table_HUMAN_trans_rna$coef_pho_kin > 0,]
table_HUMAN_trans_pro_sig = table_HUMAN_trans_pro[table_HUMAN_trans_pro$FDR_pho_kin < sig & table_HUMAN_trans_pro$coef_pho_kin > 0,]

### cis comparison ###
# look into two levels
colnames(table_HUMAN_cis_rna)[4:16] = paste("RNA",colnames(table_HUMAN_cis_rna)[4:16],sep="_")
table_HUMAN_cis_compare = merge(table_HUMAN_cis,table_HUMAN_cis_rna,by=colnames(table_HUMAN_cis_rna)[c(1:3,17:18)],all.x=T)
p = ggplot(table_HUMAN_cis_compare,aes(x=-log10(FDR_pro_kin), y=-log10(RNA_FDR_pro_kin)))
#p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.2,shape=16)
p = p + geom_text_repel(aes(label= ifelse(-log10(FDR_pro_kin)+log10(RNA_FDR_pro_kin) > 7, as.character(pair), NA)),size=2,alpha=0.8)
p = p + theme_bw() + theme_nogrid()
p = p + geom_vline(xintercept = -log10(0.05), alpha=0.3) + geom_hline(yintercept = -log10(0.05), alpha=0.3)
p = p + labs(x = "PRO -log10(FDR)", y="RNA -log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/Cis_kinase_substrate_rna_vs_pro_moel.pdf',sep ="")
ggsave(file=fn, height=3, width=4.5, useDingbats=FALSE)

### trans comparison ###
# look into two levels
colnames(table_HUMAN_trans_rna)[4:16] = paste("RNA",colnames(table_HUMAN_trans_rna)[4:16],sep="_")
colnames(table_HUMAN_trans_pro)[4:16] = paste("PRO",colnames(table_HUMAN_trans_pro)[4:16],sep="_")
table_HUMAN_trans_compare = merge(table_HUMAN_trans,table_HUMAN_trans_rna,by=colnames(table_HUMAN_trans_rna)[c(1:3,17:18)],all.x=T)
table_HUMAN_trans_compare = merge(table_HUMAN_trans_compare,table_HUMAN_trans_pro,by=colnames(table_HUMAN_trans_pro)[c(1:3,17:18)],all.x=T)

p = ggplot(table_HUMAN_trans_compare,aes(x=-log10(PRO_FDR_pho_kin), y=-log10(RNA_FDR_pho_kin)))
#p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.2,shape=16)
p = p + geom_text_repel(aes(label= ifelse(-log10(PRO_FDR_pho_kin)+log10(RNA_FDR_pho_kin) > 3.5, as.character(pair), NA)),size=2,alpha=0.8)
p = p + theme_bw() + theme_nogrid()
p = p + geom_vline(xintercept = -log10(0.05), alpha=0.3) + geom_hline(yintercept = -log10(0.05), alpha=0.3)
p = p + labs(x = "PRO -log10(FDR)", y="RNA -log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/trans_kinase_substrate_rna_vs_pro_moel.pdf',sep ="")
ggsave(file=fn, height=3, width=4.5, useDingbats=FALSE)

#cor.test(x=-log10(table_HUMAN_trans_compare$FDR_pho_kin), y=-log10(table_HUMAN_trans_compare$PRO_FDR_pho_kin))
p = ggplot(table_HUMAN_trans_compare,aes(x=-log10(FDR_pho_kin), y=-log10(PRO_FDR_pho_kin)))
#p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.1,shape=16)
p = p + geom_text_repel(aes(label= ifelse(-log10(FDR_pho_kin)+log10(PRO_FDR_pho_kin) > 8, as.character(pair), NA)),size=2,alpha=0.8)
p = p + theme_bw() + theme_nogrid()
p = p + geom_vline(xintercept = -log10(0.05), alpha=0.3) + geom_hline(yintercept = -log10(0.05), alpha=0.3)
p = p + labs(x = "PHO -log10(FDR)", y="PRO -log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/trans_kinase_substrate_pro_vs_pho_moel.pdf',sep ="")
ggsave(file=fn, height=3, width=4.5, useDingbats=FALSE)

##### venn diagrams #####
library(VennDiagram)

compare = list(Cis_RNA = table_HUMAN_cis_rna_sig$pair, Cis_Pro = table_HUMAN_cis_sig$pair)
cat("CIS:Enrichment of quantitative-associated pairs against others in the universe of tested pairs:", length(intersect(table_HUMAN_cis_rna$pair,table_HUMAN_cis$pair)), "\n")
cat("Enrichment against RNA pairs\n")
run_fisher(universe = intersect(table_HUMAN_cis_rna$pair,table_HUMAN_cis$pair),pathway_genes = unique(table_HUMAN_cis_rna_sig$pair),sig_genes = unique(table_HUMAN_cis_sig$pair))

venn.plot = venn.diagram(
  x = compare, 
  filename = NULL,
  lwd = 4,
  fill = c("royalblue1", "orange1"),
  alpha = 0.2,
  label.col = "black",
  cex = 1.5,
  cat.cex = 1,
  cat.default.pos = c("text"), 
  margin=0.2
)
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/Cis_kinase_substrate_rna_vs_pro.pdf',sep ="")
pdf(file=fn)
grid.draw(venn.plot)
dev.off()

compare_trans = list(Trans_RNA = table_HUMAN_trans_rna_sig$pair, Trans_Pro = table_HUMAN_trans_pro_sig$pair, Trans_Pho = table_HUMAN_trans_sig$pair)
#all_pairs = intersect(intersect(table_HUMAN_trans_rna$pair, table_HUMAN_trans_pro$pair), table_HUMAN_trans$pair)
all_pairs = table_HUMAN_trans$pair
cat("TRANS:Enrichment of quantitative-associated pairs against others in the universe of tested pairs:", length(all_pairs), "\n")
cat("Enrichment against RNA pairs\n")
run_fisher(universe = all_pairs,pathway_genes = unique(table_HUMAN_trans_rna_sig$pair),sig_genes = unique(table_HUMAN_trans_sig$pair))
cat("Enrichment against PRO pairs\n")
run_fisher(universe = all_pairs,pathway_genes = unique(table_HUMAN_trans_pro_sig$pair),sig_genes = unique(table_HUMAN_trans_sig$pair))
cat("Enrichment of PRO against RNA pairs\n")
run_fisher(universe = all_pairs,pathway_genes = unique(table_HUMAN_trans_rna_sig$pair),sig_genes = unique(table_HUMAN_trans_pro_sig$pair))

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
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/Trans_kinase_substrate_rna_vs_pro_vs_pho.pdf',sep ="")

pdf(file=fn)
grid.draw(venn.plot)
dev.off()
