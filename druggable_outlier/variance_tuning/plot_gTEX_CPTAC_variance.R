##### plot_variance.R #####
# Kuan-lin Huang @ WashU 2016 Jan
# use gTEX expression data as normal background distribution
# plot relative variance in samples across genes
# command run: Rscript calc_CPTAC_rna_gene_variance.R > CPTAC_rpkm_tissue_gene_stat.txt

##### dependencies #####
setwd("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier/variance_tuning")
source("/Users/khuang/bin/LIB_exp.R")

drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv_list.txt_list', header=FALSE, stringsAsFactors = F)
druggable = as.vector(t(drugList))

gTEX_RNA = read.table(quote="", header=TRUE, sep="\t", file="gTEX_rpkm_tissue_gene_stat_log2.txt")
CPTAC_RNA = read.table(quote="", header=TRUE, sep="\t", file="CPTAC_rpkm_tissue_gene_stat_log2.txt")

gTEX_RNA$cohort = "gTEX"
CPTAC_RNA$cohort = "CPTAC"

all_RNA = rbind(gTEX_RNA, CPTAC_RNA)

# plot gTEX distribution across tissue type
for (gene in druggable){
  #all_RNA_gene = all_RNA[all_RNA$gene %in% gene,]
  fn = paste(pd, gene, 'gTEX_mean_SD.pdf',sep ="_")
  gTEX_RNA_a = gTEX_RNA[gTEX_RNA$gene %in% gene,]
  p = ggplot(data = gTEX_RNA_a, aes(x=tissue, y=mean, fill=tissue)) + 
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
  p = p + 
    theme(legend.position = "none",axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=10, angle = 90), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())
  p = p + theme_nogrid()
  p
  ggsave(file=fn, height=2.5, width=8, useDingbats=FALSE)
}


gTEX_breast = gTEX_RNA[gTEX_RNA$tissue=="Breast",]
CPTAC_breast = CPTAC_RNA[CPTAC_RNA$tissue=="CPTAC_BRCA",]
breast_m = merge(gTEX_breast,CPTAC_breast, by = "gene")
breast_m$Expression = "low"
breast_m[breast_m$median.x + breast_m$median.y >= 9,]$Expression = "high"
breast_m$label = breast_m$gene

fn = paste(pd, "breast_SD_gTEX-CPTAC.pdf",sep ="_")
p = ggplot(data = breast_m , aes(x=SD.x, y=SD.y, colour=Expression))
p = p + geom_point(alpha=0.1, size=1.5) + xlab("gTEX RNA SD") + ylab("CPTAC RNA SD") 
p = p + theme_bw() 
p = p + geom_text(aes(label=ifelse(as.character(gene) %in% druggable,as.character(gene),'')),hjust=-0.05,just=0,size=3,alpha=0.7)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn, height=5, width=5, useDingbats=FALSE)

fn = paste(pd, "breast_IQR_gTEX-CPTAC.pdf",sep ="_")
p = ggplot(data = breast_m , aes(x=IQR.x, y=IQR.y, colour=Expression)) 
p = p + geom_point(alpha=0.1, size=1.5) + xlab("gTEX RNA IQR") + ylab("CPTAC RNA IQR") 
p = p + theme_bw() 
p = p + geom_text(aes(label=ifelse(as.character(gene) %in% druggable,as.character(gene),'')),hjust=-0.05,just=0,size=3,alpha=0.7)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn, height=5, width=5, useDingbats=FALSE)