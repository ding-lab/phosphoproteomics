##### find_plot_pQTL.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run pQTL analysis for 3 cancer types and plot the result using RPPA, compare to proteomic results

##### dependencies #####
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

setwd(paste(baseD,"pan3can_analysis/mutation_impact/", sep=""))
source("/Users/khuang/bin/LIB_exp.R")
source("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/mutation_impact/mutation_impact.R")

pd = paste(pd,"all", sep="_")
# system("mkdir logs")
# logFile = paste("logs/", date, "_mutation_impact_analysis.log", sep="")
# sink(file=logFile)
#sink(file=NULL)
kinaseList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/2014-12-05_CPTAC_Kinase.MATRIX.v3b5_sheet1_genes.list', header=FALSE, stringsAsFactors = F)
drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv_list.txt_list', header=FALSE, stringsAsFactors = F)
cancer_genes127 = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/Kandoth_127SMG_list_brca_refseq_ncbi.txt.list', header=FALSE, stringsAsFactors = F)

kinome = as.vector(t(kinaseList))
druggable = as.vector(t(drugList))
druggable = as.vector(t(cancer_genes127))
SMGs = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C","TP53", "NF1", "KRAS", "BRCA1", "BRCA2", "CDK12",
         "TP53","KRAS","APC","PIK3CA","SMAD4")

##### BRCA #####
### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep=""))
brcaGenes = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C")
BRCA_mut_g = BRCA_mut[row.names(BRCA_mut) %in% brcaGenes,]

### RPPA ###
BRCA_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_RPPA_formatted.txt",sep=""))
BRCA_RPPA_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_RPPA,name="BRCA RPPA")

##### OV #####
### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_SOMATIC_formatted.txt",sep=""))
ovGenes = c("TP53", "NF1", "KRAS", "BRCA1", "BRCA2", "CDK12")
OV_mut_g = OV_mut[row.names(OV_mut) %in% ovGenes,]

### RPPA ###
OV_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_RPPA_formatted.txt",sep=""))
OV_RPPA_diff_exp = find_diff_exp(OV_mut_g,OV_RPPA,name="OV RPPA")

##### CRC #####
### Mutation matrix ###
CRC_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_SOMATIC_formatted.txt",sep=""))
crcGenes = c("TP53","KRAS","APC","PIK3CA","SMAD4")
CRC_mut_g = CRC_mut[row.names(CRC_mut) %in% crcGenes,]
### RPPA ###
CRC_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_RPPA_formatted.txt",sep=""))
CRC_RPPA_diff_exp = find_diff_exp(CRC_mut_g,CRC_RPPA,name="CRC RPPA")


##### merge everything and plot heatmap #####

### fold change
BRCA_RPPA_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_RPPA_diff_exp$fold_change)))
colnames(BRCA_RPPA_diff_exp_fc_m)[3] = "BRCA_RPPA"
OV_RPPA_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_RPPA_diff_exp$fold_change)))
colnames(OV_RPPA_diff_exp_fc_m)[3] = "OV_RPPA"
CRC_RPPA_diff_exp_fc_m = melt(as.matrix(do.call(cbind,CRC_RPPA_diff_exp$fold_change)))
colnames(CRC_RPPA_diff_exp_fc_m)[3] = "CRC_RPPA"

fc = merge(BRCA_RPPA_diff_exp_fc_m, OV_RPPA_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
fc = merge(fc, CRC_RPPA_diff_exp_fc_m, by = c("Var1","Var2"), all=T)

fc_m = melt(fc, id.var=c("Var1","Var2"))
colnames(fc_m)[4] = "FC"

## fdr
BRCA_RPPA_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_RPPA_diff_exp$t_fdr)))
colnames(BRCA_RPPA_diff_exp_fdr_m)[3] = "BRCA_RPPA"
OV_RPPA_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_RPPA_diff_exp$t_fdr)))
colnames(OV_RPPA_diff_exp_fdr_m)[3] = "OV_RPPA"
CRC_RPPA_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,CRC_RPPA_diff_exp$t_fdr)))
colnames(CRC_RPPA_diff_exp_fdr_m)[3] = "CRC_RPPA"

fdr = merge(BRCA_RPPA_diff_exp_fdr_m, OV_RPPA_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
fdr = merge(fdr, CRC_RPPA_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)

fdr_m = melt(fdr, id.var=c("Var1","Var2"))
colnames(fdr_m)[4] = "FDR"

fc_fdr = merge(fc_m, fdr_m, by=c("Var1","Var2","variable"))

## plot

fdr.colors=c("NA", "#000000")
min_d = min(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
max_d = max(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
bound = max(c(max_d, -min_d))
fc_fdr$sig = as.numeric(as.character(fc_fdr$FDR)) <= 0.1
fc_fdr$FC_2 = as.numeric(fc_fdr$FC)
fc_fdr[!is.na(fc_fdr$FC_2) & fc_fdr$FC_2>=2,]$FC_2=2
fc_fdr[!is.na(fc_fdr$FC_2) & fc_fdr$FC_2<=-2,]$FC_2=-2
fc_fdr.v = fc_fdr[rowSums(is.na(fc_fdr))<4,]

fn = paste(pd, 'merged_diff_exp_RPPA.pdf',sep ="_")
p = ggplot(data=fc_fdr.v)
#p = p + facet_grid(.~variable, scales = "fixed", space = "free", drop=TRUE)
p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=16))
p = p + coord_equal()
p
ggsave(file=fn, width=10.5, height=24, useDingbats=FALSE)


### plot selected markers ###

for (i in 3:5){
  fdr[,i] = as.numeric(as.character(fdr[,i]))
}
markers = unique(fdr[rowSums(fdr[,c(3:5)]<=0.1, na.rm=T) >=1,]$Var1)
fc_fdr_s = fc_fdr.v[fc_fdr.v$Var1 %in% as.character(markers),]
#fc_fdr_s2 = fc_fdr_s[fc_fdr_s$FDR<=0.2,]

fn = paste(pd, 'merged_diff_exp_RPPA_fdr0.1in1.pdf',sep ="_")
p = ggplot(data=fc_fdr_s)
#p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
p = p + coord_fixed()
p
ggsave(file=fn, width=10.5, height=10, useDingbats=FALSE)

##### plot collective violin plot #####
### make violin plots for selected mutated gene-marker combo###

plot_diff_exp_violin_c = function (Mgene, gene){
  BRCA_mut_gene = data.frame(t(BRCA_mut[Mgene,]))
  BRCA_mut_gene$carrier = as.character(BRCA_mut_gene[,1]) != "wt" & as.character(BRCA_mut_gene[,1]) != "silent"
  Bgene = sub("_.*","",row.names(BRCA_RPPA))
  BRCA_RPPA_gene = t(BRCA_RPPA[which(Bgene==gene),])
  # merge by sample
  BRCA_RPPA_merge = merge(BRCA_mut_gene, BRCA_RPPA_gene, by = "row.names")
  BRCA_RPPA_merge$data = "BRCA_RPPA"
  
  CRC_mut_gene = data.frame(t(CRC_mut[Mgene,]))
  CRC_mut_gene$carrier = as.character(CRC_mut_gene[,1]) != "wt" & as.character(CRC_mut_gene[,1]) != "silent"
  Cgene = sub("_.*","",row.names(CRC_RPPA))
  CRC_RPPA_gene = t(CRC_RPPA[which(Cgene==gene),])
  CRC_RPPA_merge = merge(CRC_mut_gene, CRC_RPPA_gene, by = "row.names")
  CRC_RPPA_merge$data = "CRC_RPPA"
  
  OV_mut_gene = data.frame(t(OV_mut[Mgene,]))
  OV_mut_gene$carrier = as.character(OV_mut_gene[,1]) != "wt" & as.character(OV_mut_gene[,1]) != "silent"
  Ogene = sub("_.*","",row.names(OV_RPPA))
  OV_RPPA_gene = t(OV_RPPA[which(Ogene==gene),])
  OV_RPPA_merge = merge(OV_mut_gene, OV_RPPA_gene, by = "row.names")
  OV_RPPA_merge$data = "OV_RPPA"
  
  colnames(BRCA_RPPA_merge)[4]=gene
  colnames(CRC_RPPA_merge)[4]=gene
  colnames(OV_RPPA_merge)[4]=gene
  
  gene_all_lvl = rbind(BRCA_RPPA_merge,CRC_RPPA_merge,OV_RPPA_merge)
  if (ncol(gene_all_lvl)>5){gene_all_lvl = gene_all_lvl[,c(1,2,3,4,ncol(gene_all_lvl))]}
  colnames(gene_all_lvl) = c("Sample","Mutation_Type","Mutation_Status","Expression","Dataset")
  
  # plot violin plots faceted by marker genes
  fn = paste(pd, Mgene, gene, "RPPA_mutational_impact_violin.pdf", sep="_")
  p = ggplot(data=gene_all_lvl)
  p = p + facet_grid(.~Dataset)
  p = p + geom_violin(aes(x=Mutation_Status, y=Expression, fill=Mutation_Status),alpha=0.5) + guides(fill=FALSE) 
  p = p + geom_jitter(aes(x=Mutation_Status, y=Expression, color=Mutation_Type)) #+ geom_point(aes(x=Status, y=value)) 
  p = p + labs(x = paste(Mgene,"Mutation Status"), y = paste(gene, "Expression")) + theme_bw()
  p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
  p
  ggsave(file=fn, width=14.5, limitsize=FALSE, useDingbats=FALSE)
}

# plot any pair of interest
if (FALSE){
  plot_diff_exp_violin_c("TP53","TP53")
  plot_diff_exp_violin_c("TP53","ESR1")
  plot_diff_exp_violin_c("TP53","IGF1R")
  plot_diff_exp_violin_c("TP53","GATA3")
  plot_diff_exp_violin_c("TP53","EGFR")
  plot_diff_exp_violin_c("TP53","CDH1")
  plot_diff_exp_violin_c("GATA3","EGFR")
  plot_diff_exp_violin_c("TP53","CHEK2")
  plot_diff_exp_violin_c("NF1","TP53")
  plot_diff_exp_violin_c("CDH1","CDH1")
  plot_diff_exp_violin_c("CDH1","PDGFRB")
  plot_diff_exp_violin_c("KRAS","MAP2K1")
}