##### phosphosite_summary.R #####
# Kuan-lin Huang @ WashU 2016 Feb
# cross-correlation between different phosphosites

baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/"

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/phospho_summary/")
source("/Users/khuang/bin/LIB_exp.R")

brca_file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
ov_file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
brca_pho = read.table(row.names=1, header=TRUE, sep="\t", file= brca_file)
ov_pho = read.table(row.names=1, header=TRUE, sep="\t", file= ov_file)

num_sites = length(unique(c(row.names(brca_pho),row.names(ov_pho))))
brca_genes = gsub("(.*?)\\.(.).*","\\1",row.names(brca_pho))
ov_genes = gsub("(.*?)\\.(.).*","\\1",row.names(ov_pho))
num_genes = length(unique(c(brca_genes,ov_genes)))
cat("number of unique genes characterized:", num_genes)
cat("number of unique sites characterized:", num_sites)

# BRCA
# non_unique = gsub("\\.[0-9]+","",row.names(brca_pho))
# brca_residue = gsub("(.*)\\.(.).*","\\2",non_unique)
non_unique = gsub(".*:","",row.names(brca_pho))
brca_residue = substring(non_unique, 1, 1)
BRCA = table(brca_residue)
BRCA

brca_pho_10NA = brca_pho[rowSums(!is.na(brca_pho))>=10,]
# non_unique = gsub("\\.[0-9]+","",row.names(brca_pho_10NA))
# brca_residue_10NA = gsub("(.*)\\.(.).*","\\2",non_unique)
non_unique = gsub(".*:","",row.names(brca_pho_10NA))
brca_residue_10NA = substring(non_unique, 1, 1)
BRCA10NA = table(brca_residue_10NA)
BRCA10NA

# OV
# non_unique = gsub("\\.[0-9]+","",row.names(ov_pho))
# OV_residue = gsub("(.*)\\.(.).*","\\2",non_unique)
non_unique = gsub(".*:","",row.names(ov_pho))
OV_residue = substring(non_unique, 1, 1)
OV = table(OV_residue)
OV

OV_pho_10NA = ov_pho[rowSums(!is.na(ov_pho))>=10,]
# non_unique = gsub("\\.[0-9]+","",row.names(OV_pho_10NA))
# OV_residue_10NA = gsub("(.*)\\.(.).*","\\2",non_unique)
non_unique = gsub(".*:","",row.names(OV_pho_10NA))
OV_residue_10NA = substring(non_unique, 1, 1)
OV10NA = table(OV_residue_10NA)
OV10NA

# both
#residue_table = rbind(BRCA,BRCA10NA,OV,OV10NA)
residue_table = rbind(BRCA10NA,OV10NA)
residue_table_m = melt(residue_table)
colnames(residue_table_m) = c("Cancer","Residue","Number_of_Phosphosites")

fn = paste(pd, 'phosphosite_brca.ov_sty_summary.pdf',sep ="_")
p = ggplot(residue_table_m,aes(x=Cancer, y=Number_of_Phosphosites, fill=Residue))
p = p + geom_bar(stat="identity") + theme_bw() + theme_nogrid()
#p = p + labs(x = "Data type", y="# of phosphosites")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
ggsave(file=fn, height = 3, width=4, useDingbats=FALSE)

##### investigate expression vs. variance #####
RTK_file = read.table("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/reference_files/RTKs_list.txt")
RTK = as.vector(t(RTK_file))

BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))
gene = row.names(BRCA_RNA)
gene_exp = rowMeans(BRCA_RNA,na.rm=T)
RNA_gene_exp = cbind(gene, gene_exp)

site = row.names(brca_pho_10NA)
gene = gsub(":.*","",site)
site_SD = rowSds(as.matrix(brca_pho_10NA),na.rm=T)
pho_site_gene_SD = cbind(gene, site, site_SD)

RNA_gene_exp_pho_site_gene_SD = merge(pho_site_gene_SD, RNA_gene_exp, by="gene", all.x=T)

RNA_gene_exp_pho_site_gene_SD$gene_exp = as.numeric(as.character(RNA_gene_exp_pho_site_gene_SD$gene_exp))
RNA_gene_exp_pho_site_gene_SD$site_SD = as.numeric(as.character(RNA_gene_exp_pho_site_gene_SD$site_SD))

RNA_gene_exp_pho_site_gene_SD_RTK = RNA_gene_exp_pho_site_gene_SD[RNA_gene_exp_pho_site_gene_SD$gene %in% RTK,] 

p = ggplot(RNA_gene_exp_pho_site_gene_SD_RTK,aes(x=gene_exp, y=site_SD))
#p = p + facet_grid(self~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
#p = p + geom_point(alpha=0.5)
#p = p + geom_text(aes(label= ifelse(gene_exp > 10 & site_SD > 2.6, as.character(site), NA)),size=2,alpha=0.5)
p = p + geom_text(aes(label= ifelse(gene %in% RTK, as.character(site), NA), color = gene),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
#p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + theme(legend.position="none")
p = p + labs(x = "Average gene expression (log2)", y="Phosphosite expression SD")
p = p + expand_limits(x = 0, y = 0) + xlim(0,15)
p
fn = paste(pd, 'phosphositeSD_vs_geneExp_RTK.pdf',sep ="_")
ggsave(file=fn, height=5, width=5, useDingbats=FALSE)

## violin plot of expression
BRCA_Clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file=paste(baseD,"BRCA/BRCA_clinical_summary.txt",sep=""))
brca_pho_10NA_RTK = brca_pho_10NA[gene %in% RTK,]
brca_pho_10NA_RTK_t = t(brca_pho_10NA_RTK)
brca_pho_10NA_RTK_t_m = melt(brca_pho_10NA_RTK_t)
colnames(brca_pho_10NA_RTK_t_m) = c("Sample","Site","Expression")
BRCA_Clin_t = as.data.frame(t(BRCA_Clin))
BRCA_Clin_t$Sample = row.names(BRCA_Clin_t)
brca_pho_10NA_RTK_clin = merge(brca_pho_10NA_RTK_t_m, BRCA_Clin_t, by="Sample")
brca_pho_10NA_RTK_clin$Gene = gsub(":.*","",brca_pho_10NA_RTK_clin$Site)

# plotting
get.clinical.scale = function() {
  # Set1 colors
  colors = c(NA, "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  # use Perou's intrinsic subtype colors instead
  colors = c(NA, "#101010", NA, "#636363", "red", "pink", "#ffeda0", "#deebf7", "#3182bd", "#FFFF33", "#A65628", "#F781BF", "#999999") #positive is dark grey       
  
  color.names = c("wt","mut","negative", "positive", "Basal", "Her2", "CLDN low", "LumB", "LumA") # add lumA and check color
  #color.names = c("negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB")
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="Intrinsic Subtype", values=colors)
  
  return(clinical.color.scale)
}

p = ggplot(data=brca_pho_10NA_RTK_clin)
#p = p + facet_wrap(Gene~Site, nrow=5)
p = p + facet_grid(.~Gene,scale="free", space="free")
#p = p + geom_boxplot(aes(x=Site, y=Expression, fill=NULL),alpha=0.1, outlier.shape = NA) 
p = p + geom_jitter(aes(x=Site, y=Expression, color=pam50), size=1, alpha=0.4) #+ geom_point(aes(x=Status, y=value)) 
#p = p + geom_text(aes(x=Species, y=Pho, label = ifelse(outlier,as.character(Sample),NA)),size=2)
p = p + theme_nogrid() + guides(fill=FALSE) 
p = p + labs(x = "", y = "Phosphosite expression")
p = p + get.clinical.scale()
p = p + theme(text = element_text(colour="black", size=16), axis.ticks.x = element_blank(),
              axis.text.x = element_text(colour="black", size=5,angle=90,vjust=0.5),#element_blank(), axis.title.x = element_blank(),  
              axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
p = p + theme(legend.position="bottom") + ylim(-7.2,7.2)
p
fn = paste(pd, "BRCA_phosite_RTK_by_subtype.pdf", sep="_")
ggsave(file=fn, width=12, height=6,useDingbats=FALSE)