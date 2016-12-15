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
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted_amino_acid.txt",sep=""))
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
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_SOMATIC_formatted_amino_acid.txt",sep=""))
ovGenes = c("TP53", "NF1", "KRAS", "BRCA1", "BRCA2", "CDK12")
OV_mut_g = OV_mut[row.names(OV_mut) %in% ovGenes,]

### Proteome ###
OV_Pho = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA.tsv")
OV_Pho_genes = gsub("\\:.*","",row.names(OV_Pho))
OV_Pho_c = OV_Pho[OV_Pho_genes %in% druggable,]
#OV_Pho_c_10NA = OV_Pho_c[rowSums(!is.na(OV_Pho_c))>=10,]
# may want to limit analysis to the most widely observed phosphosite for each protein
OV_Pho_diff_exp = find_diff_exp(OV_mut_g,OV_Pho_c,name="OV Phosphosites")

# draw violin plot for each significant correlation

##### plot collective violin plot #####
### make violin plots for selected mutated gene-marker combo###

#   PERHAPS INCLUDE VENKATA STYLE HEATMAP AT SOME POINT # 
plot_mut_heatmap = function (Mgene, gene){
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
  
  colnames(gene_all_lvl) = c("Sample","Mutation_Status","Amino_acid_change","Dataset","Phosphosite","Expression")
  
  # this line is optional: whether to remove samples with all NA or not
  gene_all_lvl = gene_all_lvl[!is.na(gene_all_lvl$Expression),]
  
  # so we get all the samples with the same aa change
  gene_all_lvl$Amino_acid_change_u = make.names(gene_all_lvl$Amino_acid_change, unique=T)
  # reorder mutation and phosphorylation
  #pos1 = sub("p.[A-Z]","",gene_all_lvl$Amino_acid_change)
  #gene_all_lvl$mut_pos = sub("[A-Z]$","",pos1)
   
  temp = gsub("[[:alpha:]]", "", gene_all_lvl$Amino_acid_change)
  gene_all_lvl$mut_pos = gsub("\\.","",temp)
  
  
  #gene_all_lvl$pho_pos = sub(".*\\.[a-z]","",gene_all_lvl$Phosphosite)
  gene_all_lvl$pho_pos = sub(".*[a-z]([0-9]+)","\\1",gene_all_lvl$Phosphosite)
  
  gene_all_lvl$Amino_acid_change_u = factor(gene_all_lvl$Amino_acid_change_u, levels = gene_all_lvl$Amino_acid_change_u[order(gene_all_lvl$mut_pos)])
  gene_all_lvl$Phosphosite = factor(gene_all_lvl$Phosphosite, levels = gene_all_lvl$Phosphosite[order(gene_all_lvl$pho_pos)])
  
  num_mut = length(unique(gene_all_lvl$Amino_acid_change_u))
  num_sites = length(unique(gene_all_lvl$Phosphosite))
  bound = max(max(gene_all_lvl$Expression,na.rm=T), -min(gene_all_lvl$Expression,na.rm=T)) 
  #min_bound = min(gene_all_lvl$Expression,na.rm=T) 
  
  fn = paste(pd, Mgene, gene, "mutational_impact_pos.pdf", sep="_")
  p = ggplot(data=gene_all_lvl)
  p = p + facet_grid(Dataset ~ ., drop=T, scales = "free_x", space = "free")
  #p = p + facet_wrap(~Dataset, ncol=1, scales = "free_x", drop=TRUE)
  p = p + geom_tile(aes(x=Amino_acid_change_u, y=Phosphosite, fill=as.numeric(Expression)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-bound,bound))
  #p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
  p = p + labs(x = paste(Mgene,"Mutation"), y = "")#paste(gene, "Expression"))
  p = p + theme_bw() + theme_nogrid() +
    theme(axis.title = element_text(size=14), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=8), 
          axis.text.y = element_text(colour="black", size=12), axis.ticks = element_blank())#element_text(colour="black", size=16))
  p = p + coord_equal()
  #p = p + theme(legend.position="none")
  p
  ggsave(file=fn, height=num_sites*1, width = num_mut*0.2, useDingbats=FALSE)
  
  ## plot a version 
  # plot violin plots faceted by marker genes
#   fn = paste(pd, Mgene, gene, "mutational_impact_phosphosite_violin.pdf", sep="_")
#   p = ggplot(data=gene_all_lvl)
#   p = p + facet_grid(Phosphosite~Dataset)
#   p = p + geom_violin(aes(x=Mutation_Status, y=Expression, fill=Mutation_Status),alpha=0.5) + guides(fill=FALSE) 
#   #p = p + geom_point(aes(x=Mutation_Status, y=Expression, color=Mutation_Type, position = position_jitter(w = 0.1, h = 0)))
#   p = p + geom_jitter(aes(x=Mutation_Status, y=Expression, color=Mutation_Type)) #+ geom_point(aes(x=Status, y=value)) 
#   #p = p + geom_dotplot(aes(binaxis = "y",x=Mutation_Status, y=as.numeric(Expression), color=Mutation_Type, 
#   #                         stackdir = "center",binwidth = 5, stackgroups = TRUE,method = "histodot")) #+ geom_point(aes(x=Status, y=value)) 
#   p = p + labs(x = paste(Mgene,"Mutation Status"), y = paste(gene, "Expression"))
#   p = p + ylim(c(min_bound,max_bound))
#   p = p + theme_bw() + theme_nogrid()
#   p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
#                 axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
#   p = p + scale_color_manual(values = set1[1:length(unique(gene_all_lvl$Mutation_Type))])
#   p
#   ggsave(file=fn, height = num_sites*1.5, limitsize=FALSE, useDingbats=FALSE)
  
}

# plot any pair of interest
if (FALSE){
  
  plot_mut_heatmap("TP53","ESR1")
  plot_mut_heatmap("TP53","IGF1R")
  plot_mut_heatmap("TP53","GATA3")
  plot_mut_heatmap("TP53","EGFR")
  plot_mut_heatmap("TP53","CDH1")
  plot_mut_heatmap("GATA3","EGFR")
  plot_mut_heatmap("TP53","CHEK2")
  plot_mut_heatmap("NF1","TP53")
  plot_mut_heatmap("CDH1","CDH1")
  plot_mut_heatmap("CDH1","PDGFRB")
  plot_mut_heatmap("KRAS","MAP2K1")
  plot_mut_heatmap("PIK3CA","PIK3CA")
  plot_mut_heatmap("PIK3CA","AKT1")
  plot_mut_heatmap("PIK3CA","AKT2")
  plot_mut_heatmap("PIK3CA","AKT3")
  plot_mut_heatmap("PIK3CA","MTOR")
  
  # cis
  plot_mut_heatmap("TP53","TP53")
  plot_mut_heatmap("GATA3","GATA3")
  plot_mut_heatmap("EGFR","EGFR")
  plot_mut_heatmap("PIK3CA","PIK3CA")
  
  
  plot_mut_heatmap("IGF1R","IGF1R")
}

gene = "TP53"
Mgene = "TP53"
