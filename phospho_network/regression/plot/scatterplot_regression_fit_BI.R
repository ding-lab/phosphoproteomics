### scatterplot_regression_fit.R ### 
# Kuan Huang @ 2017 March

library(ggplot2)
library(grid)
require(plyr)

lm_eqn = function(df){
  m = glm(pho_sub_norm ~ pro_kinase, data = df);
  eq <- substitute(italic(pho_sub_norm) == a + b %.% italic(pro_kinase),
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2)));
  as.character(as.expression(eq));
}

get.clinical.scale = function() {
  colors = c(NA, "#101010", NA, "#636363", "#CE2427","#EF5591","#FFFF33","#8FBCE5","#423996","#50A547") #positive is dark grey       
  color.names = c("wt","Normal","negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB","LumA","Lymphoma")
  
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="mutation", values=colors)
  
  return(clinical.color.scale)
}

draw_cis_correlation = function(gene){
  # merge and preprocess data
  pro_data_g = data.frame(t(pro_data[pro_data$X==gene,-1]))
  colnames(pro_data_g) = "kinase"
  pro_data_g$sample = row.names(pro_data_g)
  
  pho_data_g = data.frame(t(pho_data[pho_data$Gene==gene,-c(1,2,3)]))
  colnames(pho_data_g) = gsub("\\.",":",colnames(pho_data_g))
  
  # pick the top sites
  table_HUMAN_cis_g = table_HUMAN_cis[table_HUMAN_cis$KINASE==gene & table_HUMAN_cis$FDR_pro_kin<sig,]
  table_HUMAN_cis_g_top3 = gsub(":.*:",":",table_HUMAN_cis_g$pair[order(table_HUMAN_cis_g$FDR_pro_kin)][1:3])
  pho_data_g_sele = pho_data_g[,colnames(pho_data_g) %in% table_HUMAN_cis_g_top3]
  pho_data_g_sele$sample = row.names(pho_data_g)
  
  data = merge(pro_data_g,pho_data_g_sele,by="sample")
  data$sample = gsub("[0-9][0-9]TCGA","01A",data$sample)
  data$subtype = "Normal"
  for (cohort in c("Her2","LumA","LumB","Basal")) {
    subtype_sample = na.omit(colnames(clinical)[clinical[1,]==cohort])
    data$subtype[data$sample %in% subtype_sample] = cohort
  }
  data$subtype[data$subtype=="Her2"]="HER2-E"
  
  data_m = melt(data, id.vars=c("kinase","sample","subtype"))
  data_m$gene = gene
  p = ggplot(data_m,aes(x=kinase, y=value, color = subtype))
  p = p + facet_grid(variable~gene)
  p = p + geom_point(stroke=0, alpha=0.8) 
  p = p + theme_bw()
  p = p + get.clinical.scale()
  p = p + labs(x = "", y="")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(legend.position="bottom")
  p
  fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/KSpair/',gene,'_cis_site_regressions.pdf',sep ="")
  ggsave(file=fn, height = 4, width = 4, useDingbats=FALSE)
}

draw_trans_correlation = function(gene){
  # merge and preprocess data
  pho_gdata_g = data.frame(t(pho_gdata[pho_gdata$X==gene,-1]))
  colnames(pho_gdata_g) = "kinase"
  pho_gdata_g$sample = row.names(pho_gdata_g)
  
  # pick the top sites
  table_HUMAN_trans_g = table_HUMAN_trans[table_HUMAN_trans$KINASE==gene & table_HUMAN_trans$FDR_pho_kin<sig,]
  table_HUMAN_trans_g_top3 = paste(table_HUMAN_trans_g$SUBSTRATE[order(table_HUMAN_trans_g$FDR_pho_kin)][1:3],
                                   table_HUMAN_trans_g$SUB_MOD_RSD[order(table_HUMAN_trans_g$FDR_pho_kin)][1:3],sep=":")
  pho_data_g_sele = data.frame(t(pho_data[rownames(pho_data) %in% table_HUMAN_trans_g_top3,-c(1:3)]))
  pho_data_g_sele$sample = row.names(pho_data_g)
  
  data = merge(pho_gdata_g,pho_data_g_sele,by="sample")
  data$sample = gsub("[0-9][0-9]TCGA","01A",data$sample)
  data$subtype = "Normal"
  for (cohort in c("Her2","LumA","LumB","Basal")) {
    subtype_sample = na.omit(colnames(clinical)[clinical[1,]==cohort])
    data$subtype[data$sample %in% subtype_sample] = cohort
  }
  data$subtype[data$subtype=="Her2"]="HER2-E"
  
  data_m = melt(data, id.vars=c("kinase","sample","subtype"))
  data_m$gene = gene
  p = ggplot(data_m,aes(x=kinase, y=value, color = subtype))
  p = p + facet_grid(variable~gene)
  p = p + geom_point(stroke=0, alpha=0.8) 
  p = p + theme_bw()
  p = p + get.clinical.scale()
  p = p + labs(x = "", y="")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(legend.position="bottom")
  p
  fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/KSpair/',gene,'_trans_site_regressions.pdf',sep ="")
  ggsave(file=fn, height = 4, width = 4, useDingbats=FALSE)
}

draw_trans_correlation_pair = function(kinase,substrate,sites, h){
  # merge and preprocess data
  pho_gdata_g = data.frame(t(pho_gdata[pho_gdata$X==kinase,-1]))
  colnames(pho_gdata_g) = "kinase"
  pho_gdata_g$sample = row.names(pho_gdata_g)
  
  pho_data_g = data.frame(t(pho_data[pho_data$Gene==substrate,-c(1,2,3)]))
  #colnames(pho_data_g) = gsub("\\.",":",colnames(pho_data_g))
  colnames(pho_data_g) = gsub(".*\\.","",colnames(pho_data_g))
  pho_data_g_sele = pho_data_g[,sites]
  pho_data_g_sele$sample = row.names(pho_data_g_sele)
  
  data = merge(pho_gdata_g,pho_data_g_sele,by="sample")
  data$sample = gsub("[0-9][0-9]TCGA","01A",data$sample)
  data$subtype = "Normal"
  for (cohort in c("Her2","LumA","LumB","Basal")) {
    subtype_sample = na.omit(colnames(clinical)[clinical[1,]==cohort])
    data$subtype[data$sample %in% subtype_sample] = cohort
  }
  data$subtype[data$subtype=="Her2"]="HER2-E"
  
  data_m = melt(data, id.vars=c("kinase","sample","subtype"))
  data_m$gene = kinase
  p = ggplot(data_m,aes(x=kinase, y=value, color = subtype))
  p = p + facet_grid(variable~gene)
  p = p + geom_point(stroke=0, alpha=0.8) 
  p = p + theme_bw()
  p = p + get.clinical.scale()
  p = p + labs(x = "", y="")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(legend.position="bottom")
  p
  fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/KSpair/',kinase,'_vs_',substrate,'_trans_site_regressions.pdf',sep ="")
  ggsave(file=fn, width = 4, height = h, useDingbats=FALSE)
}

draw_cis_correlation_all = function(gene){
  # merge and preprocess data
  pro_data_g = data.frame(t(pro_data[pro_data$X==gene,-1]))
  colnames(pro_data_g) = "kinase"
  pro_data_g$sample = row.names(pro_data_g)
  
  pho_data_g = data.frame(t(pho_data[pho_data$Gene==gene,-c(1,2,3)]))
  #colnames(pho_data_g) = gsub("\\.",":",colnames(pho_data_g))
  colnames(pho_data_g) = gsub(".*\\.","",colnames(pho_data_g))
  pho_data_g$sample = row.names(pho_data_g)
  
  data = merge(pro_data_g,pho_data_g,by="sample")
  data$sample = gsub("[0-9][0-9]TCGA","01A",data$sample)
  data$subtype = "Normal"
  for (cohort in c("Her2","LumA","LumB","Basal")) {
    subtype_sample = na.omit(colnames(clinical)[clinical[1,]==cohort])
    data$subtype[data$sample %in% subtype_sample] = cohort
  }
  data$subtype[data$subtype=="Her2"]="HER2-E"
  
  data_m = melt(data, id.vars=c("kinase","sample","subtype"))
  data_m$gene = gene
  p = ggplot(data_m,aes(x=kinase, y=value, color = subtype))
  p = p + facet_grid(variable~gene)
  p = p + geom_point(stroke=0, alpha=0.8) 
  p = p + theme_bw()
  p = p + get.clinical.scale()
  p = p + labs(x = "", y="")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(legend.position="bottom")
  p
  fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/KSpair/',gene,'_cis_site_regressions_all.pdf',sep ="")
  ggsave(file=fn, width = 4, useDingbats=FALSE)
}

draw_cis_correlation_sele = function(gene,sites,h=5){
  # merge and preprocess data
  pro_data_g = data.frame(t(pro_data[pro_data$X==gene,-1]))
  colnames(pro_data_g) = "kinase"
  pro_data_g$sample = row.names(pro_data_g)
  
  pho_data_g = data.frame(t(pho_data[pho_data$Gene==gene,-c(1,2,3)]))
  #colnames(pho_data_g) = gsub("\\.",":",colnames(pho_data_g))
  colnames(pho_data_g) = gsub(".*\\.","",colnames(pho_data_g))
  pho_data_g_sele = pho_data_g[,sites]
  pho_data_g_sele$sample = row.names(pho_data_g_sele)
  
  data = merge(pro_data_g,pho_data_g_sele,by="sample")
  data$sample = gsub("[0-9][0-9]TCGA","01A",data$sample)
  data$subtype = "Normal"
  for (cohort in c("Her2","LumA","LumB","Basal")) {
    subtype_sample = na.omit(colnames(clinical)[clinical[1,]==cohort])
    data$subtype[data$sample %in% subtype_sample] = cohort
  }
  data$subtype[data$subtype=="Her2"]="HER2-E"
  
  data_m = melt(data, id.vars=c("kinase","sample","subtype"))
  data_m$gene = gene
  p = ggplot(data_m,aes(x=kinase, y=value, color = subtype))
  p = p + facet_grid(variable~gene)
  p = p + geom_point(stroke=0, alpha=0.8) 
  p = p + theme_bw()
  p = p + get.clinical.scale()
  p = p + labs(x = "", y="")
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(legend.position="bottom")
  p
  fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/KSpair/',gene,'_cis_site_regressions_sele.pdf',sep ="")
  ggsave(file=fn, width = 4, height=h,useDingbats=FALSE)
}
##### data files #####

table_HUMAN_cis = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_cis_sig_fam.txt",sep = ""))
table_HUMAN_trans = read.delim(stringsAsFactors = F,paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_trans_sig_fam.txt",sep = ""))
clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))

# input -------------------------------------------------------------------

# calculate mean phosphorylation level (gene) in subtypes ------------------------
HUMAN_pho_g = paste(baseD,"pan3can_shared_data/BRCA/BRCA77_unimodal_phosphoproteome-ratio-norm_collapsed.txt",sep="")
pho_gdata = read.delim(HUMAN_pho_g)

# calculate mean phosphorylation level (site) in subtypes ------------------------
pho_data = read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA77_unimodal_phosphoproteome-ratio-norm_wGpos_cleaned.txt",sep=""))
row.names(pho_data) = gsub("\\.",":",make.names(gsub(":NP.*:",":",pho_data$Gene.site),unique=T))

# calculate mean protein expression level in subtypes ------------------------
pro_data <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt",sep=""))

### main ###
# # draw regression for the top genes
# cis_genes = names(sort(table(table_HUMAN_cis$KINASE), decreasing = T))[1:15]
# for (gene in cis_genes){ draw_cis_correlation(gene)}
# 
# trans_genes = names(sort(table(table_HUMAN_trans$KINASE), decreasing = T))[1:15]
# for (gene in trans_genes){ draw_trans_correlation(gene)}

draw_cis_correlation_all("AKT1")
draw_cis_correlation_all("MAPK3")
draw_cis_correlation_all("MAP2K6")


RAF1_sites = c("S29","S43","T260")
draw_cis_correlation_sele("RAF1",RAF1_sites,4)
draw_trans_correlation_pair("MAPK3","RAF1",RAF1_sites,4)
draw_trans_correlation_pair("PAK1","RAF1",RAF1_sites,4)

AKT1_sites = c("S129","T308","T450")
draw_cis_correlation_sele("AKT1",AKT1_sites,4)

ERBB2_sites = c("S1151","S998","Y1248")
draw_cis_correlation_sele("ERBB2",ERBB2_sites,4)

draw_trans_correlation_pair("MAP3K5","MAP2K6")
