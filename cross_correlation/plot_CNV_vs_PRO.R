##### plot_CNV_vs_PRO.R #####
# Kuan-lin Huang @ WashU 2017 Feb
# Plotting and analysis for cross level correlation between CNV, RNA, and Proteome data

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/cross_correlation/")
source("/Users/khuang/bin/LIB_exp.R")

baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"
cancer_genes = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/cgenes_and_druggable.list', header=FALSE, stringsAsFactors = F)
k_s_table = read.table(header=T, stringsAsFactors = F, quote = "", fill=T, sep = "\t","/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt")
genes = unique(c(as.vector(t(cancer_genes)),k_s_table$GENE))

### input files ###

RTK_file = read.table("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/reference_files/RTKs_list.txt")
RTK = as.vector(t(RTK_file))

BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted_amino_acid.txt",sep=""))
BRCA_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_CNV_formatted_normalized.txt",sep=""))
BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted_normalized.txt",sep="")) # use not normalized, doesn't matter for Spearman
BRCA_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
# BRCA_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_RPPA_formatted.txt",sep=""))
# BRCA_RPPA_f = RPPA2geneName(BRCA_RPPA)
# 
# CRC_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_CNV_formatted_normalized.txt",sep=""))
# CRC_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_mRNA_formatted.txt",sep="")) # use not normalized, doesn't matter for Spearman
# CRC_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_PRO_formatted_normalized.txt",sep=""))
# CRC_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_RPPA_formatted.txt",sep=""))
# 
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_SOMATIC_formatted_amino_acid.txt",sep=""))
OV_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_CNV_formatted_normalized.txt",sep=""))
OV_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_mRNA_formatted_normalized.txt",sep="")) # use not normalized, doesn't matter for Spearman
OV_JHU_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_PRO_formatted_normalized.txt",sep=""))
OV_PNNL_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))
OV_PNNL_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))
# OV_JHU_Gly = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_GLY_formatted_normalized.txt",sep=""))
# OV_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_RPPA_formatted.txt",sep=""))

BRCA_CNV.RNA = read.table(header=TRUE, sep="\t", file="figures/2017-02-05/2017-02-05_KH_BRCA\ CNV\ vs.\ RNA_p_correlation.txt")
BRCA_RNA.Pro = read.table(header=TRUE, sep="\t", file="figures/2017-02-05/2017-02-05_KH_BRCA\ RNA\ vs.\ PRO_p_correlation.txt")
BRCA_Pro.Pho = read.table(header=TRUE, sep="\t", file="figures/2017-02-05/2017-02-05_KH_BRCA\ PRO\ vs.\ PHO_p_correlation.txt")
BRCA_CNV.Pro = read.table(header=TRUE, sep="\t", file="figures/2017-02-05/2017-02-05_KH_BRCA\ CNV\ vs.\ PRO_p_correlation.txt")

OV_CNV.RNA = read.table(header=TRUE, sep="\t", file="figures/2017-02-05/2017-02-05_KH_OV\ CNV\ vs.\ RNA_p_correlation.txt")
OV_RNA.Pro = read.table(header=TRUE, sep="\t", file="figures/2017-02-05/2017-02-05_KH_OV\ RNA\ vs.\ PRO_p_correlation.txt")
OV_Pro.Pho = read.table(header=TRUE, sep="\t", file="figures/2017-02-05/2017-02-05_KH_OV\ PRO\ vs.\ PHO_p_correlation.txt")
OV_CNV.Pro = read.table(header=TRUE, sep="\t", file="figures/2017-02-05/2017-02-05_KH_OV\ CNV\ vs.\ PRO_p_correlation.txt")

OV_PNNL_CNV = OV_CNV[,colnames(OV_CNV) %in% colnames(OV_PNNL_Pro)]
OV_PNNL_RNA = OV_RNA[,colnames(OV_RNA) %in% colnames(OV_PNNL_Pro)]

##### merge #####
colnames(BRCA_CNV.RNA)[2] = "CNV.RNA_correlation_coef"
colnames(BRCA_CNV.RNA)[3] = "CNV.RNA_P"
colnames(BRCA_RNA.Pro)[2] = "RNA.Pro_correlation_coef"
colnames(BRCA_RNA.Pro)[3] = "RNA.Pro_P"
colnames(BRCA_Pro.Pho)[2] = "Pro.Pho_correlation_coef"
colnames(BRCA_Pro.Pho)[3] = "Pro.Pho_P"
colnames(BRCA_CNV.Pro)[2] = "CNV.Pro_correlation_coef"
colnames(BRCA_CNV.Pro)[3] = "CNV.Pro_P"
BRCA_correlation_merge = merge(BRCA_CNV.RNA,BRCA_RNA.Pro,by="gene")
BRCA_correlation_merge = merge(BRCA_correlation_merge,BRCA_CNV.Pro,by="gene")
BRCA_correlation_merge = merge(BRCA_correlation_merge,BRCA_Pro.Pho,by="gene")
#BRCA_all_correlation = rbind(BRCA_CNV.RNA,BRCA_RNA.Pro,BRCA_Pro.Pho,BRCA_CNV.Pro)

colnames(OV_CNV.RNA)[2] = "CNV.RNA_correlation_coef"
colnames(OV_CNV.RNA)[3] = "CNV.RNA_P"
colnames(OV_RNA.Pro)[2] = "RNA.Pro_correlation_coef"
colnames(OV_RNA.Pro)[3] = "RNA.Pro_P"
colnames(OV_Pro.Pho)[2] = "Pro.Pho_correlation_coef"
colnames(OV_Pro.Pho)[3] = "Pro.Pho_P"
colnames(OV_CNV.Pro)[2] = "CNV.Pro_correlation_coef"
colnames(OV_CNV.Pro)[3] = "CNV.Pro_P"
OV_correlation_merge = merge(OV_CNV.RNA,OV_RNA.Pro,by="gene")
OV_correlation_merge = merge(OV_correlation_merge,OV_CNV.Pro,by="gene")
OV_correlation_merge = merge(OV_correlation_merge,OV_Pro.Pho,by="gene")

# define CNV events: which genes have variable CNVs
# currently use CNV threshold of 0.5 (CNA) and -0.5 (CND)
CNV_count = rowSums(BRCA_CNV > 0.5 | BRCA_CNV < -0.5, na.rm=T)
CNA_value = rowSums(BRCA_CNV > 0.5, na.rm=T)
CND_value = rowSums(BRCA_CNV < -0.5, na.rm=T)
CNV_category = CNA_value - CND_value
BRCA_CNV_counts = cbind(row.names(BRCA_CNV),data.frame(CNV_count),data.frame(CNV_category))
colnames(BRCA_CNV_counts)[1] = "gene"

CNV_count = rowSums(OV_PNNL_CNV > 0.5 | OV_PNNL_CNV < -0.5, na.rm=T)
CNA_value = rowSums(OV_PNNL_CNV > 0.5, na.rm=T)
CND_value = rowSums(OV_PNNL_CNV < -0.5, na.rm=T)
CNV_category = CNA_value - CND_value
OV_CNV_counts = cbind(row.names(OV_PNNL_CNV),data.frame(CNV_count),data.frame(CNV_category))
colnames(OV_CNV_counts)[1] = "gene"

BRCA_fn = paste(pd,"BRCA_cross_lvl_corr_CNV.tsv",sep="_")
OV_fn = paste(pd,"OV_cross_lvl_corr_CNV.tsv",sep="_")
write.table(BRCA_correlation_merge, row.names = F, quote=F, sep = '\t', file=BRCA_fn)
write.table(OV_correlation_merge, row.names = F, quote=F, sep = '\t', file=OV_fn)
##### plot to look into correlation vs. SV etc #####
# CNV count vs. correlation

BRCA_correlation_merge_g = BRCA_correlation_merge[BRCA_correlation_merge$gene %in% genes,]
OV_correlation_merge_g = OV_correlation_merge[OV_correlation_merge$gene %in% genes,]
BRCA_correlation_merge_g_sd = merge(BRCA_correlation_merge_g,BRCA_CNV_counts,by="gene")
OV_correlation_merge_g_sd = merge(OV_correlation_merge_g,OV_CNV_counts,by="gene")
BRCA_correlation_merge_g_sd$cancer = "Breast cancer (n=77)"
OV_correlation_merge_g_sd$cancer = "Ovarian cancer (n=84)"
cancers_correlation_merge_g_sd = rbind(BRCA_correlation_merge_g_sd,OV_correlation_merge_g_sd)
cancers_correlation_merge_g_sd$CNV_category_plot = cancers_correlation_merge_g_sd$CNV_category
cancers_correlation_merge_g_sd$CNV_category_plot[cancers_correlation_merge_g_sd$CNV_category_plot > 30] = 30
cancers_correlation_merge_g_sd$CNV_category_plot[cancers_correlation_merge_g_sd$CNV_category_plot < -30] = -30

p = ggplot(cancers_correlation_merge_g_sd,aes(y=CNV.Pro_correlation_coef,x=CNV_count, color=CNV_category_plot))#, color=CNVg))
p = p + facet_grid(.~cancer,scale="free_x")
p = p + geom_point(size=1, alpha=0.8)
#p = p + geom_text(aes(label=ifelse(gene %in% c("PIK3CA", "EGFR", "FOXA1", "ERBB2", "KMT2C", "PTEN", "RB1", "MAP2K4"),as.character(gene),NA)))
p = p + geom_text(aes(label=ifelse(CNV.Pro_correlation_coef> 0.4 & ((CNV_count>14 & cancer =="Breast cancer (n=77)") | (CNV_count>24 & cancer =="Ovarian cancer (n=84)")) ,as.character(gene),NA)),size=2,vjust=-0.5)
p = p + scale_color_gradientn(colours=RdBu1024, limit=c(-30,30))
p = p + theme_bw()
p = p + ylim(-0.5,1)# + xlim(0,1)
p = p + labs(x="# of samples with CNV events",y="Correlation coefficient between CNV & Protein")
p = p + expand_limits(x=0)
p
fn = paste(pd, "CNV.PRO_correlation_vs_CNVcount_genes.pdf",sep="_")
ggsave(file=fn, height=5, width=10, useDingbats=FALSE)


##### plotting #####
#@ a gene list
genes = c("PTK2", "PARP1", "MAP2K4", "ERBB2", "PRKCI", "PTEN", "RB1", "MAP2K2")
get_genes_input = function(data,genelist,dataname){
  data_genes = data.frame(t(data[genelist,]))
  data_genes.m = melt(as.matrix(data_genes))
  colnames(data_genes.m) = c("Sample","Gene",dataname)
  return(data_genes.m)
}
BRCA_CNV_genes = get_genes_input(BRCA_CNV,genes,"CNV")
BRCA_Pro_genes = get_genes_input(BRCA_Pro,genes,"Protein")
BRCA_Pho_genes = get_genes_input(BRCA_Pho,genes,"Phosphoprotein")
BRCA_mut_genes = get_genes_input(BRCA_mut,genes,"Mutation")

OV_CNV_genes = get_genes_input(OV_CNV,genes,"CNV")
OV_Pro_genes = get_genes_input(OV_PNNL_Pro,genes,"Protein")
OV_Pho_genes = get_genes_input(OV_PNNL_Pho,genes,"Phosphoprotein")
OV_mut_genes = get_genes_input(OV_mut,genes,"Mutation")

# merge by sample
BRCA_Pro_merge = merge(BRCA_CNV_genes, BRCA_Pro_genes, by = c("Sample","Gene"), all=T)
BRCA_Pro_Pho_merge = merge(BRCA_Pro_merge, BRCA_Pho_genes, by = c("Sample","Gene"), all=T)
BRCA_mut_Pho_merge = merge(BRCA_Pro_Pho_merge, BRCA_mut_genes, by = c("Sample","Gene"), all=T)
BRCA_mut_Pho_merge$cancer = "BRCA"

OV_Pro_merge = merge(OV_CNV_genes, OV_Pro_genes, by = c("Sample","Gene"), all=T)
OV_Pro_Pho_merge = merge(OV_Pro_merge, OV_Pho_genes, by = c("Sample","Gene"), all=T)
OV_mut_Pho_merge = merge(OV_Pro_Pho_merge, OV_mut_genes, by = c("Sample","Gene"), all=T)
OV_mut_Pho_merge$cancer = "OV"

mut_Pho_merge = rbind(BRCA_mut_Pho_merge,OV_mut_Pho_merge)
# mut_Pho_merge$Phosphoprotein_plot = mut_Pho_merge$Phosphoprotein # for some reason this wouldn't work
# mut_Pho_merge$Phosphoprotein_plot[!is.na(mut_Pho_merge$Phosphoprotein_plot) & mut_Pho_merge$Phosphoprotein_plot>3] = 3
# mut_Pho_merge$Phosphoprotein_plot[!is.na(mut_Pho_merge$Phosphoprotein_plot) & (mut_Pho_merge$Phosphoprotein_plot<-3)] = -3

mut_Pho_merge_plot = mut_Pho_merge[!is.na(mut_Pho_merge$CNV) & !is.na(mut_Pho_merge$Protein),]
p = ggplot(data=mut_Pho_merge_plot)
p = p + facet_grid(Gene~cancer,scale="free",space="free",drop=T)
p = p + geom_point(aes(x=CNV, y=Protein, color=Phosphoprotein), size=1, alpha=0.5) #+ geom_point(aes(x=Status, y=value)) 
p = p + geom_text(aes(x=CNV, y=Protein, label = ifelse(!is.na(Mutation) & !(Mutation %in% c("silent","wt")),as.character(Mutation),NA)),size=2,nudge_y =0.12,alpha=0.2)
p = p + scale_colour_gradientn(na.value="grey", colours=RdBu1024, limits=c(-4.2,4.2))
p = p + theme_nogrid() #+ guides(fill=FALSE) 
p = p + geom_vline(xintercept=0, alpha=0.5)
p = p + labs(x="log2(Ploidity/2)") + xlim(-3,3)
p
fn = paste(pd, "genelist_CNV_vs_pro_pho.pdf", sep="_")
ggsave(file=fn, height=10, useDingbats=FALSE)

# use ERBB2 to find phosphoprotein threshold
mut_Pho_merge_plot_ERBB2 = mut_Pho_merge_plot[mut_Pho_merge_plot$Gene =="ERBB2" & mut_Pho_merge_plot$cancer == "BRCA",]
IQR = quantile(mut_Pho_merge_plot_ERBB2$Phosphoprotein, probs=0.75, na.rm=T) - quantile(mut_Pho_merge_plot_ERBB2$Phosphoprotein, probs=0.25, na.rm=T) 
mut_Pho_merge_plot_ERBB2$Phosphoprotein_mad = (mut_Pho_merge_plot_ERBB2$Phosphoprotein - median(mut_Pho_merge_plot_ERBB2$Phosphoprotein))/IQR

IQR = quantile(mut_Pho_merge_plot_ERBB2$Protein, probs=0.75, na.rm=T) - quantile(mut_Pho_merge_plot_ERBB2$Protein, probs=0.25, na.rm=T) 
mut_Pho_merge_plot_ERBB2$Protein_mad = (mut_Pho_merge_plot_ERBB2$Protein - median(mut_Pho_merge_plot_ERBB2$Protein))/IQR

getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
p = ggplot(data=mut_Pho_merge_plot_ERBB2)
p = p + facet_grid(Gene~cancer,scale="free",space="free",drop=T)
p = p + geom_point(aes(x=CNV, y=Phosphoprotein_mad, color=Protein_mad),alpha=0.8) #+ geom_point(aes(x=Status, y=value)) 
p = p + geom_text(aes(x=CNV, y=Phosphoprotein_mad, label = ifelse(!is.na(Mutation) & !(Mutation %in% c("silent","wt")),as.character(Mutation),NA)),size=2,nudge_y =0.12,alpha=0.2)
p = p + scale_colour_gradientn(na.value="grey", colours=getPalette(100))#, limits=c(-4.2,4.2))
p = p + theme_nogrid() #+ guides(fill=FALSE) 
p = p + geom_vline(xintercept=0, alpha=0.5)
p = p + geom_vline(xintercept=0.5, alpha=0.2)
p = p + geom_hline(yintercept=1, alpha=0.2)
#p = p + geom_vline(xintercept=0, alpha=0.5)
p = p + labs(x="log2(Ploidity/2)", y = "Phosphorylation score") + xlim(-3,3)
p
fn = paste(pd, "genelist_CNV_vs_pho_pro_ERBB2.pdf", sep="_")
ggsave(file=fn, useDingbats=FALSE)

# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
# #getPalette = colorRampPalette(c("#fed976","#e31a1c","#800026"))
# p = ggplot(data=mut_Pho_merge)
# p = p + facet_grid(cancer~Gene,scale="free",space="free",drop=T)
# p = p + geom_point(aes(x=Protein, y=Phosphoprotein, colour=log2(CNV)), size=1) #+ geom_point(aes(x=Status, y=value)) 
# p = p + geom_text(aes(x=Protein, y=Phosphoprotein, label = ifelse(!(Mutation %in% c("silent","wt")),as.character(Mutation),NA)),size=2,nudge_y =0.12,alpha=0.2)
# p = p + scale_colour_gradientn(na.value=NA, colours=getPalette(100))
# p = p + theme_nogrid() + guides(fill=FALSE) 
# #p = p + geom_abline(slope=1, alpha=0.5)
# #p = p + labs(x="log2(CNV)")
# p
# fn = paste(pd, "genelist_pro_vs_pho_wCNV.pdf", sep="_")
# ggsave(file=fn, width=20, useDingbats=FALSE)