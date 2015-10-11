##### find_plot_outlier.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run outlier analysis for 3 cancer types and plot the result

##### dependencies #####
setwd("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier")
source("/Users/khuang/bin/LIB_exp.R")
source("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier/outlier.R")
system("mkdir logs")
logFile = paste("logs/", date, "_outlier_analysis.log", sep="")
sink(file=logFile)
#sink(file=NULL)
kinaseList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/2014-12-05_CPTAC_Kinase.MATRIX.v3b5_sheet1_genes.list', header=FALSE, stringsAsFactors = F)
drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv_list.txt_list', header=FALSE, stringsAsFactors = F)
kinome = as.vector(t(kinaseList))
druggable = as.vector(t(drugList))

# function to format CDAP proteome data processed by Kuan
format_pro = function(Pro.m){
  colnames(Pro.m) = sub(".Unshared.Log.Ratio","",colnames(Pro.m)) 
  Pro.m.d = Pro.m[Pro.m$Gene %in% druggable,]
  row.names(Pro.m.d) = Pro.m.d$Gene
  Pro.m.dn = Pro.m.d[,-1]
  return(Pro.m.dn)
}

format_crc = function(Pro.m){
  row.names(Pro.m) = make.names(Pro.m$Gene, uniq=T)
  Gene = Pro.m$Gene
  Pro.m = as.matrix(Pro.m[,-1])
  colnames(Pro.m) = sub(".Unshared.Spectral.Counts","",colnames(Pro.m)) 
  
  # quantile normalization using function from limma and log2 transformation: Nature CRC proteogenomics 2014
  Pro.mn = normalizeQuantiles(Pro.m,ties=T)
  Pro.mnl = log2(Pro.mn)
  Pro.m.d = Pro.mnl[Gene %in% druggable,]
  
  return(Pro.m.d)
}

# function to normalize CNV based on a log10 scale
normalize_CNV = function(CNV.m){
  colnames(CNV.m) = paste(colnames(CNV.m),".01A", sep="")
  row.names(CNV.m) = sub(" ","", row.names(CNV.m))
  CNV.n.m = as.matrix(CNV.m)
  for (i in 1:nrow(CNV.n.m)){
    CNV.n.m[i,]=log(CNV.n.m[i,]/mean(CNV.n.m[i,], na.rm=T), base=2)
  } 
  CNV.n.md = CNV.n.m[row.names(CNV.n.m) %in% druggable,]
  return(CNV.n.md)
}

# function to normalize RSEM
format_RSEM = function(RSEM.m){ # should be normalized using the 75% quantile method already
  RSEM.m.d = RSEM.m[RSEM.m$Hybridization.REF %in% druggable,]
  row.names(RSEM.m.d) = RSEM.m.d$Hybridization.REF
  RSEM.m.d = RSEM.m.d[,-1]
  RSEM.m.d.n = as.matrix(RSEM.m.d)
  for (i in 1:nrow(RSEM.m.d.n)){
    RSEM.m.d.n[i,]=log(RSEM.m.d.n[i,]/mean(RSEM.m.d.n[i,], na.rm=T), base=2)
  } 
  return(RSEM.m.d.n)
}

##### BRCA #####
### CNV ###
BRCA_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_CNV/BRCA_105_CNV.txt")
# "NA " fail, be careful next time about the new space character
BRCA_CNV.d = normalize_CNV(BRCA_CNV)
BRCA_CNV_druggable = find_outlier(BRCA_CNV.d, name = "BRCA druggable CNV")

### RNA ###
BRCA_RNA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RNA/TCGA_Breast_BI_RSEM.tsv.parsed_hugoified")
BRCA_RNA.d = format_RSEM(BRCA_RNA)
BRCA_RNA_druggable = find_outlier(BRCA_RNA.d, name = "BRCA druggable RNA")

### Proteome ###
BRCA_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Proteome_CDAP.r2/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pro.d = format_pro(BRCA_Pro)
BRCA_Pro_druggable = find_outlier(BRCA_Pro.d, name = "BRCA druggable proteome")

### Phosphoproteome ###
BRCA_Pho = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Phosphoproteome_CDAP.r2/TCGA_Breast_BI_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
BRCA_Pho.d = format_pro(BRCA_Pho)
BRCA_Pho_druggable = find_outlier(BRCA_Pho.d, name = "BRCA druggable phosphoproteome")

### all levels ###
# do this later, may be more benefitial to do the all outlier table and summarize that instead

##### OV #####
### CNV ###
OV_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_CNV/OV_173_CNV.txt")
OV_CNV.d = normalize_CNV(OV_CNV)
OV_CNV_druggable = find_outlier(OV_CNV.d, name = "OV druggable CNV")

### RNA ###
OV_RNA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RNA/TCGA_OV_RSEM.tsv.parsed_hugoified")
OV_RNA.d = format_RSEM(OV_RNA)
OV_RNA_druggable = find_outlier(OV_RNA.d, name = "OV druggable RNA")

### Proteome ###
OV_JHU_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/JHU_Proteome_CDAP.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_PNNL_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/PNNL_Proteome_CDAP.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")

OV_JHU_Pro.d = format_pro(OV_JHU_Pro)
OV_PNNL_Pro.d = format_pro(OV_PNNL_Pro)

OV_JHU_Pro_druggable = find_outlier(OV_JHU_Pro.d, name = "OV JHU druggable proteome")
OV_PNNL_Pro_druggable = find_outlier(OV_PNNL_Pro.d, name = "OV PNNL druggable proteome")

### Phosphoproteome ###
OV_Pho = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/PNNL_Phosphoproteome_CDAP.r2/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_Pho.d = format_pro(OV_Pho)
OV_Pho_druggable = find_outlier(OV_Pho.d, name = "OV PNNL druggable phosphoproteome")

### Glycoproteome ###
OV_Gly = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/JHU_Glycoproteome_CDAP.r2/TCGA_Ovarian_JHU_Glycoproteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
OV_Gly.d = format_pro(OV_Gly)
OV_Gly_druggable = find_outlier(OV_Gly.d, name = "OV JHU druggable glycoproteome")

### merging the two proteome? ###
### all levels ###

##### CRC #####
### CNV ###
CRC_CNV = read.table(na.strings="NA ",row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_CNV/CRC_88_CNV.txt")
CRC_CNV.d = normalize_CNV(CRC_CNV)
CRC_CNV_druggable = find_outlier(CRC_CNV.d, name = "CRC druggable CNV")

### RNA ###
CRC_RNA = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_RNA/TCGA_COADREAD_RSEM_combined.tsv.parsed_hugoified")
CRC_RNA.d = format_RSEM(CRC_RNA)
CRC_RNA_druggable = find_outlier(CRC_RNA.d, name = "CRC druggable RNA")

### Proteome ###
CRC_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/CRC/VU_Proteome_CDAP.r2/TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv_hugoified.unshared_spectal_counts.txt")
# spectral count: look into the paper to see how to normalize the count data
CRC_Pro.d = format_crc(CRC_Pro)
CRC_Pro_druggable = find_outlier(CRC_Pro.d, name = "CRC druggable proteome")

### all levels ###

sink(file=NULL)

##### previous codes#####
if (FALSE){

##### OUTLIER IN ITRAQ PROTEOME #####
ITRAQ = read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')
ITRAQ.d = ITRAQ[row.names(ITRAQ) %in% druggable,]
ITRAQ_druggable = find_outlier(ITRAQ.d, name = "ITRAQ druggable proteome")

##### OUTLIER IN LFQ#####
LFQ=read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned_collapsed.txt')
LFQ.d = LFQ[row.names(LFQ) %in% druggable,]
LFQ_druggable = find_outlier(LFQ.d, "LFQ druggable proteome")

### overlap between ITRAQ and LFQ
LFQ_top_outlier2=LFQ_druggable$top_outlier[,1:5]
LFQ_top_outlier_zscore2=LFQ_druggable$top_outlier_zscore[,1:5]
LFQ_top_outlier_boolean2=LFQ_druggable$top_outlier_boolean[,1:5]
top_outlier2=ITRAQ_druggable$top_outlier[,1:5]
top_outlier_zscore2=ITRAQ_druggable$top_outlier_zscore[,1:5]
top_outlier_boolean2=ITRAQ_druggable$top_outlier_boolean[,1:5]
top_overlap = merge(LFQ_top_outlier2, top_outlier2, by = "row.names", all.x=T)
colnames(top_overlap)[1] = "WHIM"
top_zscore_overlap = merge(LFQ_top_outlier_zscore2, top_outlier_zscore2, by = "row.names", all.x=T)
colnames(top_zscore_overlap)[1] = "WHIM"
top_boolean_overlap = merge(LFQ_top_outlier_boolean2, top_outlier_boolean2, by = "row.names", all.x=T)
colnames(top_boolean_overlap)[1] = "WHIM"

# plot the overlap 
top_overlap.m <- melt(top_overlap, id.var = "WHIM")
top_zscore_overlap.m <- melt(top_zscore_overlap, id.var = "WHIM")
top_boolean_overlap.m <- melt(top_boolean_overlap, id.var = "WHIM")
colnames(top_boolean_overlap.m)[3]="outlier"
colnames(top_zscore_overlap.m)[3]="outlier_score"
top_overlap.m$variable = sub("druggable proteome ","",top_overlap.m$variable)
top_zscore_overlap.m$variable = sub("druggable proteome ","",top_zscore_overlap.m$variable)
top_boolean_overlap.m$variable = sub("druggable proteome ","",top_boolean_overlap.m$variable)

pdf(paste(pd,'ITRAQ_LFQ_top5_outlier_score_proteome.pdf',sep ="_"), height=6, width=9)
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette = colorRampPalette(YlGnBu)
# YlOrRd = brewer.pal(9, "YlOrRd") 
# getPalette = colorRampPalette(YlOrRd)
outlier.colors=c("NA", "#000000")

p = ggplot()
p = p + geom_tile(data=top_zscore_overlap.m, aes(y=as.factor(WHIM), x=variable, fill=outlier_score), linetype="blank") + scale_fill_gradientn(colours=getPalette(100))
p = p + geom_tile(data=top_boolean_overlap.m, aes(y=as.factor(WHIM), x=variable, color=outlier), fill=NA, size=0.5) + scale_colour_manual(values = outlier.colors)
p = p + geom_text(data=top_overlap.m,aes(y=as.factor(WHIM), x=variable, label = value), color="red", size=3)
p = p + ylab("Sample") + xlab("Top druggable protein") + theme_bw() + 
  theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=12), axis.text.y = element_text(,colour="black", size=12))
p
dev.off()

##### OUTLIER IN ITRAQ PHOSPHO ##### 
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
# 56651 phosphosites
ITRAQpho=ITRAQpho[,-c(1,2)]
cat("Original number of ITRAQ phosphosites: 56651\n")
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQpho.na = ITRAQpho[,-c(17,18,20)]
genes = sub("-NP.*", "", row.names(ITRAQpho.na))
row.names(ITRAQpho.na) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("\\. _.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("\\._.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub(" _.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("_.*","",row.names(ITRAQpho.na)), unique=T)
ITRAQpho.na.d = ITRAQpho.na[genes %in% druggable, ] # 1167 phosphosites
cat("Druggable list filtered ITRAQ phosphosites: 1167\n")

# outlier analysis
ITRAQpho_druggable = find_outlier(ITRAQpho.na.d, "ITRAQ druggable phosphoproteome", h=12)


##### OUTLIER IN LFQ phospho #####
LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)
LFQpho$sites=sub(".*\\(","",LFQpho$phospho_site)
LFQpho$sites=sub("\\)","",LFQpho$site)
LFQpho$sites=paste(LFQpho$gene_name,LFQpho$sites, sep=".")
row.names(LFQpho) = make.names(LFQpho$sites, unique=T)
colnames(LFQpho) = sub("_P", "", colnames(LFQpho))
colnames(LFQpho) = sub("Intensity.W", "WHIM", colnames(LFQpho))
cat("Original number of LFQ phosphosites: 18229\n")
# 18229 phosphosites
LFQpho.na = LFQpho[,-c(19:25)]
genes = sub("\\..*", "", row.names(LFQpho.na))
LFQpho.na.d = LFQpho.na[genes %in% druggable, ] #331 phosphosites
cat("Druggable list filtered LFQ phosphosites: 331\n")

LFQpho_druggable = find_outlier(LFQpho.na.d, "LFQ druggable phosphoproteome", h=12)

### find extreme phosphos that are not indicated by pro 
# phospho (ITRAQpho_outlier_zscore) and pro (outlier_zscore) in iTRAQ of the same sample, see how they correlate in z-score
ITRAQpho_outlier_zscore.dt = as.data.frame(ITRAQpho_druggable$outlier_zscore)
ITRAQpho_outlier_zscore.dt$Gene = sub("\\..*", "", row.names(ITRAQpho_outlier_zscore.dt))
ITRAQpho_outlier_zscore.dt$Site = row.names(ITRAQpho_outlier_zscore.dt)
outlier_zscore.dt = as.data.frame(ITRAQ_druggable$outlier_zscore)
outlier_zscore.dt$Gene = sub("\\..*", "", row.names(outlier_zscore.dt))

ITRAQpho_outlier_zscore.m = melt(ITRAQpho_outlier_zscore.dt, id.var = c("Gene","Site"))
outlier_zscore.dt.m = melt(outlier_zscore.dt, id.var = "Gene")
pro_pho_overlap = merge(ITRAQpho_outlier_zscore.m, outlier_zscore.dt.m, by = c("Gene","variable"), all.x=T)

corPP = cor(pro_pho_overlap$value.x, pro_pho_overlap$value.y, use = 'pairwise.complete.obs', method = "pearson")
cat(paste("\n","##### Pearson correlation for protein and mRNA outlier score:",corPP, "#####\n", sep=" "))

#pro_pho_overlap$gene = rep("Other", nrow(pro_pho_overlap))
#pro_pho_overlap[pro_pho_overlap$Gene == "ERBB2",]$gene="ERBB2"
#pro_pho_overlap[pro_pho_overlap$Gene == "PAK1",]$gene="PAK1"

pdf(paste(pd,'ITRAQ_phospho_vs_pro_score.pdf', sep="_"),height=7, width =7)

p = ggplot(data = pro_pho_overlap, aes(x=value.y, y=value.x, colour=variable))#, shape=gene)) 
p = p + geom_point(alpha=0.3, size=1.5) + xlab("Protein expression outlier score") + ylab("Phosphosite expression outlier score") 
p = p + scale_shape_manual(values=c(0,16,2)) + guides(colour=FALSE)
p = p + theme_bw() + coord_fixed() #+ theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))
p = p + xlim(-5,5) + ylim(-5,5) + geom_vline(xintercept=1.5,alpha=0.2) + geom_hline(yintercept=1.5,alpha=0.2)
#p = p + geom_text(aes(label=ifelse(value.x > 3 & value.x > value.y + 2,Site,'')),hjust=-0.05,just=0,size=3,alpha=0.7)
#p = p + scale_x_continuous(expand=c(0.02,0)) + scale_y_continuous(expand=c(0.02,0))
p = p + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
#p
#dev.off()

p2 = ggplot(data = pro_pho_overlap,aes(x=value.y, colour = variable, fill=variable)) + 
  geom_density(alpha=0.1) + guides(colour=FALSE) + #geom_vline(xintercept=1.5, alpha=0.2)  + 
  xlim(-5,5) + 
  theme_bw() +
  theme0(plot.margin = unit(c(1,0,-1.52,3.5),"lines"))
#theme0(plot.margin = unit(c(1,0,-0.48,2.2),"lines"))

p3 = ggplot(data = pro_pho_overlap,aes(x=value.x, colour = variable, fill=variable)) + 
  geom_density(alpha=0.1) + guides(colour=FALSE) + #geom_vline(xintercept=1.5, alpha=0.2) + 
  xlim(-5,5) + 
  theme_bw() + 
  coord_flip()  +
  theme0(plot.margin = unit(c(0.7,0,2.8,-0.85),"lines"))
#theme0(plot.margin = unit(c(0,1,1.2,-0.48),"lines"))

grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3))

dev.off()
#http://stackoverflow.com/questions/17370460/scatterplot-with-alpha-transparent-histograms-in-r

##### RNA-Seq #####
RNA = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")
row.names(RNA)=make.names(RNA$gene, unique=T)
RNA=RNA[,-1] #16209 genes
RNA = RNA[,colnames(RNA) %in% colnames(ITRAQ)]
RNA.d = RNA[row.names(RNA) %in% druggable,]

RNA_druggable = find_outlier(RNA.d, "RNA druggable genes")#, w=12)

##### find outlier proteome expression not identified by transcriptome #####
outlier_zscore.dt = as.data.frame(unlist(ITRAQ_druggable$outlier_zscore))
outlier_zscore.dt$Gene = row.names(outlier_zscore.dt)
RNA_outlier_zscore.dt = as.data.frame(unlist(RNA_druggable$outlier_zscore))
RNA_outlier_zscore.dt$Gene = row.names(RNA_outlier_zscore.dt)

outlier_zscore.dt.m = melt(outlier_zscore.dt, id.var = "Gene")
RNA_outlier_zscore.dt.m = melt(RNA_outlier_zscore.dt, id.var = "Gene")
pro_rna_overlap = merge(RNA_outlier_zscore.dt.m, outlier_zscore.dt.m, by = c("Gene","variable"))

corPR=cor(pro_rna_overlap$value.x, pro_rna_overlap$value.y, use = 'pairwise.complete.obs', method="pearson")
cat(paste("\n","##### Pearson correlation for protein and mRNA outlier score:",corPR, "#####\n", sep=" "))

#pro_rna_overlap$gene = rep("Other", nrow(pro_rna_overlap))
#pro_rna_overlap[pro_rna_overlap$Gene == "ERBB2",]$gene="ERBB2"
#pro_rna_overlap[pro_rna_overlap$Gene == "AKT2",]$gene="AKT2"
#pro_rna_overlap[pro_rna_overlap$Gene == "AURKA",]$gene="AURKA" #found to be amplified in basal breast tumor

pdf(paste(pd,'rna_vs_pro_outlier_score.pdf', sep="_"),height=7, width =7)
p = ggplot(data = pro_rna_overlap, aes(x=value.x, y=value.y, colour=variable))#, shape=gene)) 
p = p + geom_point(alpha=0.3, size=1.5) + xlab("mRNA expression outlier score") + ylab("Protein expression outlier score")
p = p + scale_shape_manual(values=c(0,2,16)) + guides(colour=FALSE)
p = p + theme_bw() + coord_fixed() #+ theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))
p = p + xlim(-5,5) + ylim(-5,5) + geom_vline(xintercept=1.5,alpha=0.2) + geom_hline(yintercept=1.5,alpha=0.2)
p = p + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))

p2 = ggplot(data = pro_rna_overlap,aes(x=value.x, colour = variable, fill=variable)) + 
  geom_density(alpha=0.1) + guides(colour=FALSE) + #geom_vline(xintercept=1.5, alpha=0.2) + 
  xlim(-5,5) + 
  theme_bw() +
  theme0(plot.margin = unit(c(1,0,-1.52,3.5),"lines"))
  #theme0(plot.margin = unit(c(0,0,0,0),"lines"))

p3 = ggplot(data = pro_rna_overlap,aes(x=value.y, colour = variable, fill=variable)) + 
  geom_density(alpha=0.1) + guides(colour=FALSE) + #geom_vline(xintercept=1.5, alpha=0.2) +  
  xlim(-5,5) + 
  theme_bw() + 
  coord_flip()  +
  theme0(plot.margin = unit(c(0.7,0,2.8,-0.85),"lines"))
  #theme0(plot.margin = unit(c(0,0,0,0),"lines"))

grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3))
dev.off()

# remove_axis =
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank())
# p = p + remove_axis
# p2 = p2 + remove_axis
# p3= p3 + remove_axis
# p4 = p4 + remove_axis
# 
# xyplot.gb = ggplot_build(p)
# densityY.gb = ggplot_build(p2)
# densityX.gb = ggplot_build(p3)
# 
# # densityY.gb$panel$ranges[[1]]$y.range = xyplot.gb$panel$ranges[[1]]$y.range   
# # densityX.gb$panel$ranges[[1]]$x.range = xyplot.gb$panel$ranges[[1]]$x.range
# 
# xyplot.gt = ggplot_gtable(xyplot.gb)
# densityY.gt = ggplot_gtable(densityY.gb)
# densityX.gt = ggplot_gtable(densityX.gb)
# 
# #xyplot.gt_panel = subset(xyplot.gt$layout, name == "panel")
# densityY.gt$heights = xyplot.gt$heights
# densityX.gt$widths = xyplot.gt$widths
# 
# main.grob = arrangeGrob(densityY.gt, p4, xyplot.gt, densityX.gt, ncol=2, nrow=2) #widths=c(0.3,0.7), heights=c(0.7,0.3), 
# grid.draw(main.grob)
# #pdf(file=paste(pd,'rna_vs_pro_outlier_score.pdf', sep="_"), useDingbats=FALSE)
# dev.off()
# 
# pdf(paste(pd,'rna_vs_pro_outlier_score2.pdf', sep="_"))
# grid.arrange(arrangeGrob(densityY.gt,ncol=2,widths=c(3,1)),
#              arrangeGrob(xyplot.gt,densityX.gt,ncol=2,widths=c(3,1)),
#              heights=c(1,3))
##### CNV #####
CNV = read.table(row.names=1,header=TRUE, sep="\t", file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified_normalized.tsv")
CNV.d = CNV[row.names(CNV) %in% druggable,]
# find outliers
CNV_druggable = find_outlier(CNV.d, "CNV druggable genes", w=12)

### all levels: CNV, RNA, proteome, phosphoproteome ### 
# keep the 24 WHIMs in ITRAQ proteome
# all_levels=list("LFQ_druggable"=LFQ_druggable,"RNA_druggable"=RNA_druggable,"CNV_druggable"=CNV_druggable)
# top_overlap_all = merge(ITRAQpho_druggable$top_outlier[,c(1:5)], ITRAQ_druggable$top_outlier[,c(1:5)], by = "row.names", all.x=T)
# for (i in all_levels){
#   rownames(top_overlap_all)=top_overlap_all[,1]
#   top_overlap_all=top_overlap_all[,-1]
#   top_overlap_all = merge(top_overlap_all, i, by = "row.names", all.x=T)
# }
# colnames(top_overlap_all)[1] = "sample"
# 
# all_levels2=list("top_outlier_zscore"=top_outlier_zscore,"LFQ_top_outlier_zscore"=LFQ_top_outlier_zscore,"RNA_top_outlier_zscore"=RNA_top_outlier_zscore,"CNV_top_outlier_zscore"=CNV_top_outlier_zscore)
# ITRAQpho_top_outlier_zscore.n = ITRAQpho_top_outlier_zscore/mean(ITRAQpho_top_outlier_zscore)
# LFQpho_top_outlier_zscore.n = LFQpho_top_outlier_zscore/mean(LFQpho_top_outlier_zscore)
# top_overlap_all_zscore = merge(ITRAQpho_top_outlier_zscore.n, LFQpho_top_outlier_zscore.n, by = "row.names", all.x=T)
# for (i in all_levels2){
#   rownames(top_overlap_all_zscore)=top_overlap_all_zscore[,1]
#   top_overlap_all_zscore=top_overlap_all_zscore[,-1]
#   i.n=i/mean(i)
#   top_overlap_all_zscore = merge(top_overlap_all_zscore, i.n, by = "row.names", all.x=T)
# }
# colnames(top_overlap_all_zscore)[1] = "sample"
# 
# all_levels3=list("top_outlier_boolean"=top_outlier_boolean,"LFQ_top_outlier_boolean"=LFQ_top_outlier_boolean,"RNA_top_outlier_boolean"=RNA_top_outlier_boolean,"CNV_top_outlier_boolean"=CNV_top_outlier_boolean)
# top_overlap_all_boolean = merge(ITRAQpho_top_outlier_boolean, LFQpho_top_outlier_boolean, by = "row.names", all.x=T)
# for (i in all_levels3){
#   rownames(top_overlap_all_boolean)=top_overlap_all_boolean[,1]
#   top_overlap_all_boolean=top_overlap_all_boolean[,-1]
#   top_overlap_all_boolean = merge(top_overlap_all_boolean, i, by = "row.names", all.x=T)
# }
# colnames(top_overlap_all_boolean)[1] = "sample"
# 
# top_overlap_all.m <- melt(top_overlap_all, id.var = "sample")
# top_overlap_all_zscore.m <- melt(top_overlap_all_zscore, id.var = "sample")
# top_overlap_all_boolean.m <- melt(top_overlap_all_boolean, id.var = "sample")
# 
# fn=paste(pd,'all_level_top5_mzscore.pdf',sep ="_")
# YlGnBu = brewer.pal(9, "YlGnBu") 
# getPalette = colorRampPalette(YlGnBu)
# outlier.colors=c("NA", "#000000")
# 
# p = ggplot()
# p = p + geom_tile(data=top_overlap_all_zscore.m, aes(x=as.factor(sample), y=variable, fill=value), linetype="blank") + scale_fill_gradientn(colours=getPalette(100), na.value="white", name="relative modified z-score")# values=rescale(seq(0,6,by=12/99)))
# p = p + geom_tile(data=top_overlap_all_boolean.m, aes(x=as.factor(sample), y=variable, color=value), fill=NA, size=0.5) + scale_colour_manual(values = outlier.colors, name="KDE outlier")
# p = p + geom_text(data=top_overlap_all.m ,aes(x=as.factor(sample), y=variable, label = value), color="red", size=1.5, angle=90)
# p = p + xlab("sample") + ylab("top druggable targets") + theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=12), axis.text.y = element_text(colour="black", size=12))
# p
# ggsave(file=fn, width = 210, height = 400, units = "mm")


##### BRCA77 human data ##### 
###proteome###
BRCA77 = read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt')
BRCA77.d = BRCA77[row.names(BRCA77) %in% druggable,]
BRCA77_druggable = find_outlier(BRCA77.d, "BRCA77 druggable proteome", h=6, w=24)

###phosphoproteome###
BRCA77pho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
# 33239 phosphosites
row.names(BRCA77pho) = make.names(BRCA77pho$Gene.site,unique=T)
BRCA77pho.na = BRCA77pho[,-c(1,2)]
cat("Original number of BRCA77 ITRAQ phosphosites: 33239\n")
genes = sub(".NP.*", "", row.names(BRCA77pho.na))
BRCA77pho.na.d = BRCA77pho.na[genes %in% druggable, ] #651 phosphosites
row.names(BRCA77pho.na.d) = make.names(sub(".NP_\\d+_"," ",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("\\. _.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("\\._.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub(" _.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("_.*","",row.names(BRCA77pho.na.d)), unique=T)
cat("Druggable list filtered BRCA77 ITRAQ phosphosites: 651\n")

# find outliers
BRCA77pho_druggable = find_outlier(BRCA77pho.na.d, "BRCA77 druggable phosphoproteome", h=10, w=24)

# ### find extreme phosphos that are not indicated by pro 
# BRCA77pho_outlier_zscore.dt = as.data.frame(BRCA77pho_outlier_zscore)
# BRCA77pho_outlier_zscore.dt$Gene = sub("\\..*", "", row.names(BRCA77pho_outlier_zscore.dt))
# BRCA77pho_outlier_zscore.dt$Site = row.names(BRCA77pho_outlier_zscore.dt)
# BRCA77_outlier_zscore.dt = as.data.frame(BRCA77_outlier_zscore)
# BRCA77_outlier_zscore.dt$Gene = sub("\\..*", "", row.names(BRCA77_outlier_zscore.dt))
# 
# BRCA77pho_outlier_zscore.m = melt(BRCA77pho_outlier_zscore.dt, id.var = c("Gene","Site"))
# BRCA77_outlier_zscore.dt.m = melt(BRCA77_outlier_zscore.dt, id.var = "Gene")
# BRCA77_pro_pho_overlap = merge(BRCA77pho_outlier_zscore.m, BRCA77_outlier_zscore.dt.m, by = c("Gene","variable"), all.x=T)
# 
# cor(BRCA77_pro_pho_overlap$value.x, BRCA77_pro_pho_overlap$value.y, use = 'pairwise.complete.obs')
# 
# fn = paste(pd,'BRCA77_phospho_vs_pro_mzscore.pdf', sep="_")
# p = ggplot(data = BRCA77_pro_pho_overlap, aes(x=value.y, y=value.x, colour=variable, label=Site)) 
# p = p + geom_point(alpha=0.3) + xlab("druggable proteome z-score") + ylab("druggable phosphosites z-score") + theme_bw()
# p = p + xlim(-5,7.5) + ylim(-5,7.5) + geom_vline(xintercept=2,alpha=0.2) + geom_hline(yintercept=2,alpha=0.2)
# p = p + geom_text(aes(label=ifelse(value.x > 2.5 & value.x > value.y + 2 & Gene=="ERBB2",paste(variable,Site),'')),hjust=-0.05,just=0,size=3,alpha=0.7)
# p
# ggsave(file=fn, height=10, width=10)


##### test for enrichment in overlapping outlier #####
# limit ITRAQ outlier to the set that overlapped with LFQ
ITRAQ_L = ITRAQ_druggable$top_outlier_boolean[row.names(ITRAQ_druggable$top_outlier_boolean) %in% row.names(LFQ_druggable$top_outlier_boolean),]
ITRAQ_o = sum(ITRAQ_L, na.rm=T)
ITRAQ_no = dim(ITRAQ_L)[1]*dim(ITRAQ_L)[2]
LFQ_o = sum(LFQ_druggable$top_outlier_boolean, na.rm=T)
LFQ_no = dim(LFQ_druggable$top_outlier_boolean)[1]*dim(LFQ_druggable$top_outlier_boolean)[2]
overlap = 5

a = matrix(c(5,ITRAQ_o-5,LFQ_o-5,LFQ_no-(ITRAQ_o-5)-(LFQ_o-5)-5),nrow=2,ncol=2)
cat("\n##### Enrichment of overlapping iTRAQ and LFQ outliers#####\n")
fisher.test(a)

##### find outlier using joint TCGA-WHIM cohort #####
merged_proteome=merge(ITRAQ.d, BRCA77.d, by = "row.names") #63 genes
row.names(merged_proteome) = merged_proteome[,1]
merged_proteome = merged_proteome[,-1]
merged_druggable = find_outlier(merged_proteome, name="TCGA WHIM merged druggable proteome", whim_only=T)
}
