##### find_plot_pQTL.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run pQTL analysis for 3 cancer types and plot the result

####### positional effect: mutated aa vs. phospho aa ###
##### dependencies #####
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

setwd(paste(baseD,"pan3can_analysis/mutation_impact/", sep=""))
source("/Users/khuang/bin/LIB_exp.R")
source("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/mutation_impact/mutation_impact.R")

#pd = paste(pd,"druggable",sep="_")
# system("mkdir logs")
# logFile = paste("logs/", date, "_mutation_impact_analysis.log", sep="")
# sink(file=logFile)
#sink(file=NULL)
kinaseList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/2014-12-05_CPTAC_Kinase.MATRIX.v3b5_sheet1_genes.list', header=FALSE, stringsAsFactors = F)
drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv_list.txt_list', header=FALSE, stringsAsFactors = F)
#cancer_genes = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/Kandoth_127SMG_list_brca_refseq_ncbi.txt.list', header=FALSE, stringsAsFactors = F)
cancer_genes = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/cgenes_and_druggable.list', header=FALSE, stringsAsFactors = F)

kinome = as.vector(t(kinaseList))
#druggable = as.vector(t(drugList))
druggable = as.vector(t(cancer_genes))
SMGs = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C","TP53", "NF1", "KRAS", "BRCA1", "BRCA2", "CDK12",
         "TP53","KRAS","APC","PIK3CA","SMAD4")

##### BRCA #####
### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep=""))
brcaGenes = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C")
BRCA_mut_g = BRCA_mut[row.names(BRCA_mut) %in% brcaGenes,]

### Phosphosites ###
BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv")
BRCA_Pho_genes = gsub("\\:.*","",row.names(BRCA_Pho))
BRCA_Pho_c = BRCA_Pho[BRCA_Pho_genes %in% druggable,]

#BRCA_Pho_c_10NA = BRCA_Pho_c[rowSums(!is.na(BRCA_Pho_c))>=10,]
BRCA_Pho_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_Pho_c,name="BRCA Phosphosites")

### all levels ###
# do this later, may be more benefitial to do the all outlier table and summarize that instead

##### OV #####
### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_SOMATIC_formatted.txt",sep=""))
ovGenes = c("TP53", "NF1", "KRAS", "BRCA1", "BRCA2", "CDK12")
OV_mut_g = OV_mut[row.names(OV_mut) %in% ovGenes,]

### Proteome ###
OV_Pho = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA.tsv")
OV_Pho_genes = gsub("\\:.*","",row.names(OV_Pho))
OV_Pho_c = OV_Pho[OV_Pho_genes %in% druggable,]
# OV_Pho_c_10NA = OV_Pho_c[rowSums(!is.na(OV_Pho_c))>=10,]
OV_Pho_diff_exp = find_diff_exp(OV_mut_g,OV_Pho_c,name="OV Phosphosites")

### merging the two proteome? ###
### all levels ###


##### merge everything and plot #####
# plot heatmap gene by gene #

### fold change
BRCA_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$fold_change)))
colnames(BRCA_Pho_diff_exp_fc_m)[3] = "BRCA_PHO"
OV_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_Pho_diff_exp$fold_change)))
colnames(OV_Pho_diff_exp_fc_m)[3] = "OV_PHO"

fc = merge(BRCA_Pho_diff_exp_fc_m, OV_Pho_diff_exp_fc_m, by = c("Var1","Var2"), all=T)

fc_m = melt(fc, id.var=c("Var1","Var2"))
colnames(fc_m)[4] = "FC"

## fdr
BRCA_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$t_fdr)))
colnames(BRCA_Pho_diff_exp_fdr_m)[3] = "BRCA_PHO"
OV_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_Pho_diff_exp$t_fdr)))
colnames(OV_Pho_diff_exp_fdr_m)[3] = "OV_PHO"

fdr = merge(BRCA_Pho_diff_exp_fdr_m, OV_Pho_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)

fdr_m = melt(fdr, id.var=c("Var1","Var2"))
colnames(fdr_m)[4] = "FDR"

fc_fdr = merge(fc_m, fdr_m, by=c("Var1","Var2","variable"))

## plot

fdr.colors=c("NA", "#000000")
min_d = min(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
max_d = max(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
bound = max(c(max_d, -min_d))
fc_fdr$sig = as.numeric(as.character(fc_fdr$FDR)) < 0.05
fc_fdr$FC_2 = as.numeric(fc_fdr$FC)
fc_fdr[!is.na(fc_fdr$FC_2) & fc_fdr$FC_2>=2,]$FC_2=2
fc_fdr[!is.na(fc_fdr$FC_2) & fc_fdr$FC_2<=-2,]$FC_2=-2
fc_fdr.v = fc_fdr[rowSums(is.na(fc_fdr))<4,]
## dropping lvls, doesn't work yet
# fc_fdr.v = droplevels(fc_fdr.v)
# fc_fdr.v$Var1 = factor(fc_fdr.v$Var1)
# fc_fdr.v$Var2 = factor(fc_fdr.v$Var2)
# fc_fdr.v[fc_fdr.v$variable == "BRCA_PRO",]$Var2 = droplevels(fc_fdr.v[fc_fdr.v$variable == "BRCA_PRO",]$Var2) 

fn = paste(pd, 'phosphosites_merged_diff_exp_cgenes_druggable.pdf',sep ="_")
p = ggplot(data=fc_fdr.v)
#p = p + facet_grid(.~variable, scales = "fixed", space = "free", drop=TRUE)
p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR < 0.05"),values = fdr.colors)
p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=6),axis.ticks = element_blank())#element_text(colour="black", size=16))
p = p + coord_equal()
p
ggsave(file=fn, width=10.5, height=32, useDingbats=FALSE)


### plot selected markers ###
for (i in 3:4){
  fdr[,i] = as.numeric(as.character(fdr[,i]))
}
markers = unique(fdr[rowSums(fdr[,c(3:4)]<0.05, na.rm=T) >=1,]$Var1)
fc_fdr_s = fc_fdr.v[fc_fdr.v$Var1 %in% as.character(markers),]
#fc_fdr_s2 = fc_fdr_s[fc_fdr_s$FDR<=0.2,]

fn = paste(pd, 'phosphosites_merged_diff_exp_cgenes_druggable_fdr0.05in1.pdf',sep ="_")
p = ggplot(data=fc_fdr_s)
#p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR < 0.05"),values = fdr.colors)
p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
p = p + coord_fixed()
p
ggsave(file=fn, width=8, height=13, useDingbats=FALSE)

### for grant ### 
if (FALSE){
  markers = c("AR","BRAF","CDH1","CHEK2","EGFR","ERBB2","ESR1","GATA3","PIK3CA","PDGFRA","TP53")
  fc_fdr_s = fc_fdr.v[fc_fdr.v$Var1 %in% as.character(markers),]
  fc_fdr_s = fc_fdr_s[fc_fdr_s$variable %in% c("BRCA_PRO","BRCA_PHO","OV_PNNL_PRO","OV_PNNL_PHO","CRC_PRO"),]
  fc_fdr_s$variable = fc_fdr_s$variable[drop=T]
  levels(fc_fdr_s$variable) = c("BRCA_PRO","BRCA_PHO","OV_PRO","OV_PHO","CRC_PRO")
  
  fn = paste(pd, 'merged_diff_exp_cgenes_druggable_U24_label.pdf',sep ="_")
  p = ggplot(data=fc_fdr_s)
  #p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
  p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
  p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
  p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
  p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
    theme(axis.title = element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=6), 
          axis.text.y = element_text(colour="black", size=6), axis.ticks = element_blank(),
          strip.text.x = element_text(size = 6))#element_text(colour="black", size=16))
  p = p + coord_equal()
  #p = p + theme(legend.position="none")
  p
  ggsave(file=fn, width=10, height=5, useDingbats=FALSE)
}

##### plot collective violin plot #####
### make violin plots for selected mutated gene-marker combo###

#   PERHAPS INCLUDE VENKATA STYLE HEATMAP AT SOME POINT # 
plot_diff_exp_violin_c = function (Mgene, gene){
  BRCA_mut_gene = data.frame(t(BRCA_mut[Mgene,]))
  BRCA_mut_gene$carrier = as.character(BRCA_mut_gene[,1]) != "wt" & as.character(BRCA_mut_gene[,1]) != "silent"
  BRCA_Pho_gene = t(BRCA_Pho[BRCA_Pho_genes == gene,])

  # merge by sample
  BRCA_Pho_merge = merge(BRCA_mut_gene, BRCA_Pho_gene, by = "row.names")
  BRCA_Pho_merge$data = "BRCA_PHO"

  OV_mut_gene = data.frame(t(OV_mut[Mgene,]))
  OV_mut_gene$carrier = as.character(OV_mut_gene[,1]) != "wt" & as.character(OV_mut_gene[,1]) != "silent"
  OV_Pho_gene = t(OV_Pho[OV_Pho_genes == gene,])
  
  # merge by sample
  OV_Pho_merge = merge(OV_mut_gene, OV_Pho_gene, by = "row.names")
  OV_Pho_merge$data = "OV_PHO"
  
  BRCA_Pho_merge_m = melt(BRCA_Pho_merge, id =c("Row.names","carrier",Mgene,"data"))
  OV_Pho_merge_m = melt(OV_Pho_merge, id =c("Row.names","carrier",Mgene,"data"))
  
  if (ncol(BRCA_Pho_merge_m) == ncol(OV_Pho_merge_m)){
    gene_all_lvl = rbind(BRCA_Pho_merge_m,OV_Pho_merge_m)
  } else {gene_all_lvl = BRCA_Pho_merge_m}
  
  colnames(gene_all_lvl) = c("Sample","Mutation_Status","Mutation_Type","Dataset","Phosphosite","Expression")
  
  max_bound = max(gene_all_lvl$Expression,na.rm=T) + 0.5
  min_bound = min(gene_all_lvl$Expression,na.rm=T) - 0.5
  
  num_sites = length(unique(gene_all_lvl$Phosphosite))
  # plot violin plots faceted by marker genes
  fn = paste(pd, Mgene, gene, "mutational_impact_phosphosite_violin.pdf", sep="_")
  p = ggplot(data=gene_all_lvl)
  p = p + facet_grid(Phosphosite~Dataset)
  p = p + geom_violin(aes(x=Mutation_Status, y=Expression, fill=Mutation_Status),alpha=0.5) + guides(fill=FALSE) 
  #p = p + geom_point(aes(x=Mutation_Status, y=Expression, color=Mutation_Type, position = position_jitter(w = 0.1, h = 0)))
  p = p + geom_jitter(aes(x=Mutation_Status, y=Expression, color=Mutation_Type)) #+ geom_point(aes(x=Status, y=value)) 
  #p = p + geom_dotplot(aes(binaxis = "y",x=Mutation_Status, y=as.numeric(Expression), color=Mutation_Type, 
  #                         stackdir = "center",binwidth = 5, stackgroups = TRUE,method = "histodot")) #+ geom_point(aes(x=Status, y=value)) 
  p = p + labs(x = paste(Mgene,"Mutation Status"), y = paste(gene, "Expression"))
  p = p + ylim(c(min_bound,max_bound))
  p = p + theme_bw() + theme_nogrid()
  p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
  p = p + scale_color_manual(values = set1[1:length(unique(gene_all_lvl$Mutation_Type))])
  p
  ggsave(file=fn, height = num_sites*1.5, limitsize=FALSE, useDingbats=FALSE)
}

# plot any pair of interest
if (FALSE){
  
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
  plot_diff_exp_violin_c("PIK3CA","AKT1")
  plot_diff_exp_violin_c("PIK3CA","AKT2")
  plot_diff_exp_violin_c("PIK3CA","AKT3")
  plot_diff_exp_violin_c("PIK3CA","MTOR")
  
  # cis
  plot_diff_exp_violin_c("TP53","TP53")
  plot_diff_exp_violin_c("GATA3","GATA3")
  plot_diff_exp_violin_c("EGFR","EGFR")
  plot_diff_exp_violin_c("PIK3CA","PIK3CA")
  
  
  plot_diff_exp_violin_c("IGF1R","IGF1R")
  plot_diff_exp_violin_c("MET","MET")
}

gene="MDM2"
