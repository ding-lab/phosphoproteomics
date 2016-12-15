##### find_plot_pQTL.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run pQTL analysis for 3 cancer types and plot the result

####### positional effect: mutated aa vs. phospho aa ###
##### dependencies #####
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

setwd(paste(baseD,"pan3can_analysis/mutation_impact/", sep=""))
source("/Users/khuang/bin/LIB_exp.R")
source("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/mutation_impact/mutation_impact.R")


##### BRCA #####
### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted_amino_acid.txt",sep=""))
#brcaGenes = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C")
#BRCA_mut_g = BRCA_mut[row.names(BRCA_mut) %in% brcaGenes,]

### Phosphosites ###
BRCA_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
# BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv",sep=""))
#BRCA_Pho_genes = gsub("\\:.*","",row.names(BRCA_Pho))

##### OV #####
### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_SOMATIC_formatted_amino_acid.txt",sep=""))
#ovGenes = c("TP53", "NF1", "KRAS", "BRCA1", "BRCA2", "CDK12")
#OV_mut_g = OV_mut[row.names(OV_mut) %in% ovGenes,]

### Proteome ###
# OV_JHU_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_PRO_formatted_normalized.txt",sep=""))

OV_PNNL_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))
OV_PNNL_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))
# OV_PNNL_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA.tsv",sep=""))

# ##### CRC #####
# ### Mutation matrix ###
# CRC_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_SOMATIC_formatted_amino_acid.txt",sep=""))
# 
# ### Proteome ###
# CRC_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_PRO_formatted_normalized.txt",sep=""))
# CRC_Pro_c = CRC_Pro[row.names(CRC_Pro) %in% druggable,]
# # CRC_Pro_diff_exp = find_diff_exp(CRC_mut_g,CRC_Pro_c,name="CRC Proteome")

##### plot collective violin plot #####
### make violin plots for selected mutated gene-marker combo###

plot_diff_exp_violin_c = function (Mgene, gene){
  BRCA_mut_gene = data.frame(t(BRCA_mut[Mgene,]))
  #BRCA_mut_gene$carrier = as.character(BRCA_mut_gene[,1]) != "wt" & as.character(BRCA_mut_gene[,1]) != "silent"
  BRCA_Pro_gene = t(BRCA_Pro[gene,])
  BRCA_Pho_gene = t(BRCA_Pho[gene,])
  # merge by sample
  BRCA_Pro_merge = merge(BRCA_mut_gene, BRCA_Pro_gene, by = "row.names")
  BRCA_Pro_merge$data = "BRCA_PRO"
  BRCA_Pho_merge = merge(BRCA_mut_gene, BRCA_Pho_gene, by = "row.names")
  BRCA_Pho_merge$data = "BRCA_PHO"
  
#   CRC_mut_gene = data.frame(t(CRC_mut[Mgene,]))
#   #CRC_mut_gene$carrier = as.character(CRC_mut_gene[,1]) != "wt" & as.character(CRC_mut_gene[,1]) != "silent"
#   CRC_Pro_gene = t(CRC_Pro[gene,])
#   CRC_Pro_merge = merge(CRC_mut_gene, CRC_Pro_gene, by = "row.names")
#   CRC_Pro_merge$data = "CRC_PRO"
  
  OV_mut_gene = data.frame(t(OV_mut[Mgene,]))
  #OV_mut_gene$carrier = as.character(OV_mut_gene[,1]) != "wt" & as.character(OV_mut_gene[,1]) != "silent"
#   OV_JHU_Pro_gene = t(OV_JHU_Pro[gene,])
#   OV_JHU_Pro_merge = merge(OV_mut_gene, OV_JHU_Pro_gene, by = "row.names")
#   OV_JHU_Pro_merge$data = "OV_JHU_PRO"
  OV_PNNL_Pro_gene = t(OV_PNNL_Pro[gene,])
  OV_PNNL_Pro_merge = merge(OV_mut_gene, OV_PNNL_Pro_gene, by = "row.names")
  OV_PNNL_Pro_merge$data = "OV_PRO"
  OV_PNNL_Pho_gene = t(OV_PNNL_Pho[gene,])
  OV_PNNL_Pho_merge = merge(OV_mut_gene, OV_PNNL_Pho_gene, by = "row.names")
  OV_PNNL_Pho_merge$data = "OV_PHO"
  
  colnames(BRCA_Pro_merge)[4]=gene
  colnames(BRCA_Pho_merge)[4]=gene
#   colnames(CRC_Pro_merge)[4]=gene
#   colnames(OV_JHU_Pro_merge)[4]=gene
  colnames(OV_PNNL_Pro_merge)[4]=gene
  colnames(OV_PNNL_Pho_merge)[4]=gene
  
  #gene_all_lvl = rbind(BRCA_Pro_merge,BRCA_Pho_merge,CRC_Pro_merge,OV_JHU_Pro_merge,OV_PNNL_Pro_merge,OV_PNNL_Pho_merge)
  gene_all_lvl = rbind(BRCA_Pro_merge,BRCA_Pho_merge,OV_PNNL_Pro_merge,OV_PNNL_Pho_merge)
  colnames(gene_all_lvl) = c("Sample","Mutation","Expression","Dataset")
  gene_all_lvl$Mutation = as.character(gene_all_lvl$Mutation)
  temp = substr(gene_all_lvl$Mutation,1,nchar(gene_all_lvl$Mutation)-1)
  
  gene_all_lvl$mis_pos = gsub("p\\.[A-Z]","",temp)
  gene_all_lvl$mis_pos = as.numeric(gene_all_lvl$mis_pos)
  
  gene_all_lvl$type = "Truncation"
  gene_all_lvl[!is.na(gene_all_lvl$mis_pos) & !grepl("\\*",gene_all_lvl$Mutation),]$type = "Missense"

  temp = gsub("[[:alpha:]]", "", gene_all_lvl$Mutation)
  temp = gsub("\\*", "", temp)
  gene_all_lvl$mut_pos = gsub("\\.","",temp)
  gene_all_lvl$mut_pos = as.numeric(gene_all_lvl$mut_pos)
  
  gene_all_lvl$mut_pos[gene_all_lvl$Mutation=="wt"] = 0
  gene_all_lvl$type[gene_all_lvl$Mutation=="wt"] = "Wild type"
  # draw wildtype average#
  
  # plot violin plots faceted by marker genes
  fn = paste(pd, Mgene, gene, "mutational_impact_violin.pdf", sep="_")
  p = ggplot(data=gene_all_lvl)
  p = p + facet_grid(Dataset~., scale="free_y")
  p = p + geom_point(aes(x=mut_pos, y=Expression, color=type),alpha=0.1) + guides(fill=FALSE) 
  p = p + labs(x = paste(Mgene,"Mutation Position"), y = paste(gene, "Expression")) + theme_bw()
  p = p + geom_text(aes(x=mut_pos, y=Expression, label = Mutation, color=type, stringsAsFactors=FALSE), size=2,alpha=0.8)
  p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
  p
  ggsave(file=fn, width=10, useDingbats=FALSE)
  
}

plot_diff_exp_violin_brca = function (Mgene, gene){
  BRCA_mut_gene = data.frame(t(BRCA_mut[Mgene,]))
  #BRCA_mut_gene$carrier = as.character(BRCA_mut_gene[,1]) != "wt" & as.character(BRCA_mut_gene[,1]) != "silent"
  BRCA_Pro_gene = t(BRCA_Pro[gene,])
  BRCA_Pho_gene = t(BRCA_Pho[gene,])
  # merge by sample
  BRCA_Pro_merge = merge(BRCA_mut_gene, BRCA_Pro_gene, by = "row.names")
  BRCA_Pro_merge$data = "BRCA_PRO"
  BRCA_Pho_merge = merge(BRCA_mut_gene, BRCA_Pho_gene, by = "row.names")
  BRCA_Pho_merge$data = "BRCA_PHO"
  
  colnames(BRCA_Pro_merge)[4]=gene
  colnames(BRCA_Pho_merge)[4]=gene

  gene_all_lvl = rbind(BRCA_Pro_merge,BRCA_Pho_merge)
  colnames(gene_all_lvl) = c("Sample","Mutation","Expression","Dataset")
  gene_all_lvl$Mutation = as.character(gene_all_lvl$Mutation)
  temp = substr(gene_all_lvl$Mutation,1,nchar(gene_all_lvl$Mutation)-1)
  
  gene_all_lvl$mis_pos = gsub("p\\.[A-Z]","",temp)
  gene_all_lvl$mis_pos = as.numeric(gene_all_lvl$mis_pos)
  
  gene_all_lvl$type = "Truncation"
  if (sum(!is.na(gene_all_lvl$mis_pos) & !grepl("\\*",gene_all_lvl$Mutation))>0){
    gene_all_lvl[!is.na(gene_all_lvl$mis_pos) & !grepl("\\*",gene_all_lvl$Mutation),]$type = "Missense"
  }
  
  temp = gsub("[[:alpha:]]", "", gene_all_lvl$Mutation)
  temp = gsub("\\*", "", temp)
  gene_all_lvl$mut_pos = gsub("\\.","",temp)
  gene_all_lvl$mut_pos = as.numeric(gene_all_lvl$mut_pos)
  
  gene_all_lvl$mut_pos[gene_all_lvl$Mutation=="wt"] = 0
  gene_all_lvl$type[gene_all_lvl$Mutation=="wt"] = "Wild type"
  # draw wildtype average#
  
  # plot violin plots faceted by marker genes
  fn = paste(pd, Mgene, gene, "mutational_impact_violin.pdf", sep="_")
  p = ggplot(data=gene_all_lvl)
  p = p + facet_grid(Dataset~., scale="free_y")
  p = p + geom_point(aes(x=mut_pos, y=Expression, color=type),alpha=0.1) + guides(fill=FALSE) 
  p = p + labs(x = paste(Mgene,"Mutation Position"), y = paste(gene, "Expression")) + theme_bw()
  p = p + geom_text(aes(x=mut_pos, y=Expression, label = Mutation, color=type, stringsAsFactors=FALSE), size=2,alpha=0.8)
  p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
  p
  ggsave(file=fn, width=10, useDingbats=FALSE)
  
}

# plot_diff_exp_violin_pho_c = function (Mgene, site){
#   gene = NA
#   BRCA_mut_gene = data.frame(t(BRCA_mut[Mgene,]))
#   #BRCA_mut_gene$carrier = as.character(BRCA_mut_gene[,1]) != "wt" & as.character(BRCA_mut_gene[,1]) != "silent"
# #   BRCA_Pro_gene = t(BRCA_Pro[gene,])
#   BRCA_Pho_gene = t(BRCA_Pho[site,])
#   # merge by sample
# #   BRCA_Pro_merge = merge(BRCA_mut_gene, BRCA_Pro_gene, by = "row.names")
# #   BRCA_Pro_merge$data = "BRCA_PRO"
#   BRCA_Pho_merge = merge(BRCA_mut_gene, BRCA_Pho_gene, by = "row.names")
#   BRCA_Pho_merge$data = "BRCA_PHO"
#   
# #   CRC_mut_gene = data.frame(t(CRC_mut[Mgene,]))
# #   #CRC_mut_gene$carrier = as.character(CRC_mut_gene[,1]) != "wt" & as.character(CRC_mut_gene[,1]) != "silent"
# #   CRC_Pro_gene = t(CRC_Pro[gene,])
# #   CRC_Pro_merge = merge(CRC_mut_gene, CRC_Pro_gene, by = "row.names")
# #   CRC_Pro_merge$data = "CRC_PRO"
#   
#   OV_mut_gene = data.frame(t(OV_mut[Mgene,]))
#   #OV_mut_gene$carrier = as.character(OV_mut_gene[,1]) != "wt" & as.character(OV_mut_gene[,1]) != "silent"
# #   OV_JHU_Pro_gene = t(OV_JHU_Pro[gene,])
# #   OV_JHU_Pro_merge = merge(OV_mut_gene, OV_JHU_Pro_gene, by = "row.names")
# #   OV_JHU_Pro_merge$data = "OV_JHU_PRO"
# #   OV_PNNL_Pro_gene = t(OV_PNNL_Pro[gene,])
# #   OV_PNNL_Pro_merge = merge(OV_mut_gene, OV_PNNL_Pro_gene, by = "row.names")
# #   OV_PNNL_Pro_merge$data = "OV_PNNL_PRO"
#   OV_PNNL_Pho_gene = t(OV_PNNL_Pho[site,])
#   OV_PNNL_Pho_merge = merge(OV_mut_gene, OV_PNNL_Pho_gene, by = "row.names")
#   OV_PNNL_Pho_merge$data = "OV_PNNL_PHO"
#   
# #   colnames(BRCA_Pro_merge)[4]=gene
#   colnames(BRCA_Pho_merge) = c("Sample","Mutation","Expression","Dataset")
# #   colnames(CRC_Pro_merge)[4]=gene
# #   colnames(OV_JHU_Pro_merge)[4]=gene
# #   colnames(OV_PNNL_Pro_merge)[4]=gene
#   colnames(OV_PNNL_Pho_merge)= c("Sample","Mutation","Expression","Dataset")
#   
#   #gene_all_lvl = rbind(BRCA_Pro_merge,BRCA_Pho_merge,CRC_Pro_merge,OV_JHU_Pro_merge,OV_PNNL_Pro_merge,OV_PNNL_Pho_merge)
#   gene_all_lvl = rbind(BRCA_Pho_merge,OV_PNNL_Pho_merge)
#   colnames(gene_all_lvl) = c("Sample","Mutation","Expression","Dataset")
#   gene_all_lvl$Mutation = as.character(gene_all_lvl$Mutation)
#   temp = substr(gene_all_lvl$Mutation,1,nchar(gene_all_lvl$Mutation)-1)
#   
#   gene_all_lvl$mis_pos = gsub("p\\.[A-Z]","",temp)
#   gene_all_lvl$mis_pos = as.numeric(gene_all_lvl$mis_pos)
#   
#   gene_all_lvl$type = "Truncation"
#   gene_all_lvl[!is.na(gene_all_lvl$mis_pos) & !grepl("\\*",gene_all_lvl$Mutation),]$type = "Missense"
#   
#   temp = gsub("[[:alpha:]]", "", gene_all_lvl$Mutation)
#   temp = gsub("\\*", "", temp)
#   gene_all_lvl$mut_pos = gsub("\\.","",temp)
#   gene_all_lvl$mut_pos = as.numeric(gene_all_lvl$mut_pos)
#   
#   gene_all_lvl$mut_pos[gene_all_lvl$Mutation=="wt"] = 0
#   gene_all_lvl$type[gene_all_lvl$Mutation=="wt"] = "Wild type"
#   # draw wildtype average#
#   
#   gene_all_lvl = gene_all_lvl[!is.na(gene_all_lvl$Expression) & !is.na(gene_all_lvl$mut_pos),]
#   
#   # plot violin plots faceted by marker genes
#   fn = paste(pd, Mgene, site, "mutational_impact_violin.pdf", sep="_")
#   p = ggplot(data=gene_all_lvl)
#   p = p + facet_grid(Dataset~., scale="free_y")
#   p = p + geom_point(aes(x=mut_pos, y=Expression, color=type),alpha=0.1) + guides(fill=FALSE) 
#   p = p + labs(x = paste(Mgene,"Mutation Position"), y = paste(site, "Expression")) + theme_bw()
#   p = p + geom_text(aes(x=mut_pos, y=Expression, label = Mutation, color=type, stringsAsFactors=FALSE), size=2,alpha=0.8)
#   p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
#                 axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
#   p
#   ggsave(file=fn, width=10, useDingbats=FALSE)
#   
# }

# plot any pair of interest
if (FALSE){
  plot_diff_exp_violin_c("TP53","TP53")
  plot_diff_exp_violin_c("GATA3","GATA3")
  plot_diff_exp_violin_c("CDH1","PDGFRB")
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
  plot_diff_exp_violin_pho_c("TP53","TP53:NP_000537.3:s314")
}

plot_diff_exp_violin_c("TP53","TP53")
plot_diff_exp_violin_c("TP53","CHEK2")
plot_diff_exp_violin_brca("GATA3","GATA3")
plot_diff_exp_violin_brca("CDH1","PDGFRB")
plot_diff_exp_violin_brca("CDH1","PDGFRA")
plot_diff_exp_violin_brca("TP53","PLK1")
