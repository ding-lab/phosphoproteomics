##### cross_lvl_corr.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# cross level correlation between CNV, RNA, and Proteome data

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/cross_correlation/")
source("/Users/khuang/bin/LIB_exp.R")

baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

calc_corr = function (lvl1, lvl2, name="cross-data"){
  corr_m = matrix(,nrow=nrow(lvl1),ncol=3)
  s_size = length(intersect(colnames(lvl1),colnames(lvl2))) # maximum intersect possible
  # run for each gene
  for (i in 1:nrow(lvl1)){
    gene = row.names(lvl1[i,])
    corr_g = NA
    corr_p = NA
    
    lvl1_gene = t(lvl1[gene,])
    
    ## correlation test
    if (gene %in% row.names(lvl2)){
      lvl2_gene = t(lvl2[gene,])
      merged_g = merge(lvl1_gene,lvl2_gene, by="row.names")
      if (nrow(merged_g) >= s_size*0.3){ # require the genes to be observed in at least 30% of the overlapping samples
        corr_stat = try(cor.test(merged_g[,2],merged_g[,3], method = "pearson"), silent=T)
        if (!is(corr_stat,"try-error")){ # if the correlation test is carried out successfully
          corr_g = corr_stat$estimate 
          corr_p = corr_stat$p.value
        }
      }  
    }
    
    corr_m[i,] = c(gene, corr_g, corr_p)
  }
  colnames(corr_m)[3] = "P"
  colnames(corr_m)[2] = "correlation_coef"
  colnames(corr_m)[1] = "gene"
  write.table(corr_m, row.names = F, quote=F, sep = '\t', file=paste(pd,name,"p_correlation.txt",sep="_"))
  
  corr_m_df = as.data.frame(corr_m)
  corr_m_df$correlation_coef = as.numeric(as.character(corr_m_df$correlation_coef))
  medi = round(median(corr_m_df$correlation_coef, na.rm=T),digit=3)
  per_pos = round(sum(corr_m_df$correlation_coef > 0, na.rm=T)/length(corr_m_df$correlation_coef[!is.na(corr_m_df$correlation_coef)]),digit=3)
  text_sum = as.data.frame(matrix(c(-0.75, -0.75, 500, 300,medi,per_pos), ncol=3))
  colnames(text_sum) = c("X","Y","Lab")
  
  fn = paste(pd, name, 'p_correlation.pdf',sep ="_")
  p = ggplot()
  p = p + geom_histogram(data = corr_m_df,aes(x=correlation_coef), binwidth = 0.05) + theme_bw()
  p = p + geom_text(data = text_sum, aes(x=X, y=Y, label=Lab))
  p = p + xlim(-1,1) + geom_vline(xintercept = medi, color="red")
  p
  ggsave(file=fn, useDingbats=FALSE)
  return(corr_m_df)
}

calc_R_corr = function (RPPA, lvl2, name="cross-data"){
  corr_m = matrix(,nrow=nrow(RPPA),ncol=2)
  s_size = length(intersect(colnames(RPPA),colnames(lvl2))) # maximum intersect possible
  for (i in 1:nrow(RPPA)){
    marker = row.names(RPPA[i,])
    corr_g = NA
    corr_p = NA
    
    RPPA_marker = t(RPPA[marker,])
    
    gene = sub("_.*","",marker)
    gene = sub(" .*","",gene)
    
    if (gene %in% row.names(lvl2)){
      lvl2_gene = t(lvl2[gene,])
      merged_g = merge(RPPA_marker,lvl2_gene, by="row.names")
      if (nrow(merged_g) >= s_size*0.3){ # require the genes to be observed in at least 30% of the overlapping samples
        corr_stat = try(cor.test(merged_g[,2],merged_g[,3], method = "pearson"), silent=T)
        if (!is(corr_stat,"try-error")){ # if the correlation test is carried out successfully
          corr_g = corr_stat$estimate 
          corr_p = corr_stat$p.value
        }
      }
    }     
    
    corr_m[i,] = c(row.names(RPPA)[i], corr_g, corr_p)
  }
  
  colnames(corr_m)[2] = "correlation"
  colnames(corr_m)[1] = "marker"
  write.table(corr_m, row.names = F, quote=F, sep = '\t', file=paste(pd,name,"p_correlation.txt",sep="_"))
  
  corr_m_df = as.data.frame(corr_m)
  corr_m_df$correlation = as.numeric(as.character(corr_m_df$correlation))
  medi = round(median(corr_m_df$correlation, na.rm=T),digit=3)
  per_pos = round(sum(corr_m_df$correlation > 0, na.rm=T)/length(corr_m_df$correlation[!is.na(corr_m_df$correlation)]),digit=3)
  text_sum = as.data.frame(matrix(c(-0.75, -0.75, 10, 5,medi,per_pos), ncol=3))
  colnames(text_sum) = c("X","Y","Lab")
  
  fn = paste(pd, name, 'p_correlation.pdf',sep ="_")
  p = ggplot()
  p = p + geom_histogram(data = corr_m_df,aes(x=correlation), binwidth = 0.05) + theme_bw()
  p = p + geom_text(data = text_sum, aes(x=X, y=Y, label=Lab))
  p = p + xlim(-1,1) + geom_vline(xintercept = medi, color="red")
  p
  ggsave(file=fn, useDingbats=FALSE)
  return(corr_m_df)
}

RPPA2geneName = function(m){
  m$gene = sub("_.*","",row.names(m))
  m$gene = sub(" .*","",m$gene)
  return(m)
}

##### MAIN #####
# focus on protein and phospho correlations
# load data
BRCA_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
OV_PNNL_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))
OV_PNNL_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))

# calculate correlation
BRCA_Pro.Pho = calc_corr(BRCA_Pro,BRCA_Pho, name="BRCA PRO vs. PHO")
OV_Pro.Pho = calc_corr(OV_PNNL_Pro,OV_PNNL_Pho, name="OV PRO vs. PHO")

### PRO.PHO ###
Pro.Pho = merge(BRCA_Pro.Pho,OV_Pro.Pho, by="gene")
colnames(Pro.Pho) = c("gene","BRCA_corr","BRCA_P","OV_corr","OV_P")

# compare
Pro.Pho$diff = Pro.Pho$BRCA_corr - Pro.Pho$OV_corr

p = ggplot(Pro.Pho,aes(x=BRCA_corr, y = OV_corr))
p = p + geom_point(alpha=0.3) + theme_nogrid()
p = p + xlim(-1,1) + ylim(-1,1)
p = p + labs( x = "BRCA Pro vs. Pho Correlation" , y = "OV Pro vs. Pho Correlation")
p = p + geom_abline(intercept = 0, slope = 1)
p = p + geom_text(aes(label = ifelse(abs(diff) >= 1 & (as.numeric(as.character(BRCA_P)) < 0.05) & (as.numeric(as.character(OV_P)) < 0.05) ,as.character(gene),NA)))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
fn = paste(pd, 'OV_vs_BRCA_PRO.PHO_correlation.pdf',sep ="_")
ggsave(file=fn, height=5, width=5, useDingbats=FALSE)




### previous correlations ###
# BRCA_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_CNV_formatted_normalized.txt",sep=""))
# BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted.txt",sep="")) # use not normalized, doesn't matter for Spearman
# BRCA_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
# BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
# BRCA_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_RPPA_formatted.txt",sep=""))
# BRCA_RPPA_f = RPPA2geneName(BRCA_RPPA)
# 
# CRC_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_CNV_formatted_normalized.txt",sep=""))
# CRC_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_mRNA_formatted.txt",sep="")) # use not normalized, doesn't matter for Spearman
# CRC_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_PRO_formatted_normalized.txt",sep=""))
# CRC_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_RPPA_formatted.txt",sep=""))
# 
# OV_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_CNV_formatted_normalized.txt",sep=""))
# OV_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_mRNA_formatted.txt",sep="")) # use not normalized, doesn't matter for Spearman
# OV_JHU_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_PRO_formatted_normalized.txt",sep=""))
# OV_PNNL_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))
# OV_PNNL_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))
# OV_JHU_Gly = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_GLY_formatted_normalized.txt",sep=""))
# OV_RPPA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_RPPA_formatted.txt",sep=""))

##### ANALYSIS #####
# BRCA_CNV.RNA = calc_corr(BRCA_CNV,BRCA_RNA, name="BRCA CNV vs. RNA")
# BRCA_RNA.Pro = calc_corr(BRCA_RNA,BRCA_Pro, name="BRCA RNA vs. PRO")
# BRCA_Pro.Pho = calc_corr(BRCA_Pro,BRCA_Pho, name="BRCA PRO vs. PHO")
# BRCA_CNV.Pro = calc_corr(BRCA_CNV,BRCA_Pro, name="BRCA CNV vs. PRO")
# 
# OV_CNV.RNA = calc_corr(OV_CNV,OV_RNA, name="OV CNV vs. RNA")
# OV_RNA.JHU_Pro = calc_corr(OV_RNA,OV_JHU_Pro, name="OV RNA vs. JHU PRO")
# OV_RNA.PNNL_Pro = calc_corr(OV_RNA,OV_PNNL_Pro, name="OV RNA vs. PNNL PRO")
# OV_Pro.PNNL_Pho = calc_corr(OV_PNNL_Pro,OV_PNNL_Pho, name="OV PNNL PRO vs. PNNL PHO")
# #OV_Pro.JHU_Gly = calc_corr(OV_JHU_Pro,OV_JHU_Gly, name="OV JHU PRO vs. JHU GLY")
# OV_CNV.JHU_Pro = calc_corr(OV_CNV,OV_JHU_Pro, name="OV CNV vs. JHU PRO")
# OV_CNV.PNNL_Pro = calc_corr(OV_CNV,OV_PNNL_Pro, name="OV CNV vs. PNNL PRO")
# 
# CRC_CNV.RNA = calc_corr(CRC_CNV,CRC_RNA, name="CRC CNV vs. RNA")
# CRC_RNA.Pro = calc_corr(CRC_RNA,CRC_Pro, name="CRC RNA vs. PRO")
# CRC_CNV.Pro = calc_corr(CRC_CNV,CRC_Pro, name="CRC CNV vs. PRO")
# 
# ##### RPPA ANALYSIS #####
# BRCA_Pro.RPPA = calc_R_corr(BRCA_RPPA,BRCA_Pro, name="BRCA RPPA vs. PRO")
# BRCA_Pho.RPPA = calc_R_corr(BRCA_RPPA,BRCA_Pho, name="BRCA RPPA vs. PHO")
# 
# OV_JHU_Pro.RPPA = calc_R_corr(OV_RPPA,OV_JHU_Pro, name="OV RPPA vs. JHU PRO")
# OV_PNNL_Pro.RPPA = calc_R_corr(OV_RPPA,OV_PNNL_Pro, name="OV RPPA vs. PNNL PRO")
# OV_Pho.RPPA = calc_R_corr(OV_RPPA,OV_PNNL_Pho, name="OV RPPA vs. PNNL PHO")
# 
# CRC_Pro.RPPA = calc_R_corr(CRC_RPPA,CRC_Pro, name="CRC RPPA vs. PRO")
# 
# ##### merge #####
# 
# 
# ### for U24 grant ###
# if (FALSE){
# 
#   ### CNV.RNA ###
#   BRCA_CNV.RNA$dataset = "BRCA"
#   OV_CNV.RNA$dataset = "OV"
#   CRC_CNV.RNA$dataset = "CRC"
#   
#   all_CNV.RNA = rbind(BRCA_CNV.RNA,OV_CNV.RNA,CRC_CNV.RNA)
#   all_CNV.RNA$dataset = as.factor(all_CNV.RNA$dataset)
#   
#   fn = paste(pd, 'pan3can_CNV.RNA_correlation.pdf',sep ="_")
#   p = ggplot(all_CNV.RNA,aes(x=correlation, fill=dataset))
#   p = p + geom_histogram(binwidth = 0.05, alpha=0.5, position="identity") + theme_bw()
#   p = p + xlim(-1,1) #+ geom_vline(xintercept = medi, color="red")
#   p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
#   p
#   ggsave(file=fn, height=5, width=8, useDingbats=FALSE)
#   
#   ### RNA.PRO ###
#   BRCA_RNA.Pro$dataset = "BRCA"
#   OV_RNA.PNNL_Pro$dataset = "OV"
#   CRC_RNA.Pro$dataset = "CRC"
#   
#   all_RNA.PRO = rbind(BRCA_RNA.Pro,OV_RNA.PNNL_Pro,CRC_RNA.Pro)
#   all_RNA.PRO$dataset = as.factor(all_RNA.PRO$dataset)
#   
#   fn = paste(pd, 'pan3can_RNA.PRO_correlation.pdf',sep ="_")
#   p = ggplot(all_RNA.PRO,aes(x=correlation, fill=dataset))
#   p = p + geom_histogram(binwidth = 0.05, alpha=0.5, position="identity") + theme_bw()
#   p = p + xlim(-1,1) #+ geom_vline(xintercept = medi, color="red")
#   p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
#   p
#   ggsave(file=fn, height=5, width=8, useDingbats=FALSE)
#   
#   ### CNV.Pro ###
#   BRCA_CNV.Pro$dataset = "BRCA"
#   OV_CNV.PNNL_Pro$dataset = "OV"
#   CRC_CNV.Pro$dataset = "CRC"
#   
#   all_CNV.PRO = rbind(BRCA_CNV.Pro,OV_CNV.PNNL_Pro,CRC_CNV.Pro)
#   all_CNV.PRO$dataset = as.factor(all_CNV.PRO$dataset)
#   
#   fn = paste(pd, 'pan3can_CNV.PRO_correlation.pdf',sep ="_")
#   p = ggplot(all_CNV.PRO,aes(x=correlation, fill=dataset))
#   p = p + geom_histogram(binwidth = 0.05, alpha=0.5, position="identity") + theme_bw()
#   p = p + xlim(-1,1) #+ geom_vline(xintercept = medi, color="red")
#   p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
#   p
#   ggsave(file=fn, height=5, width=8, useDingbats=FALSE)
#   
#   ### PRO.PHO ###
#   BRCA_Pro.Pho$dataset = "BRCA"
#   OV_Pro.PNNL_Pho$dataset = "OV"
#   
#   all_Pro.Pho = rbind(BRCA_Pro.Pho,OV_Pro.PNNL_Pho)
#   all_Pro.Pho$dataset = as.factor(all_Pro.Pho$dataset)
#   
#   fn = paste(pd, 'pan3can_PRO.PHO_correlation.pdf',sep ="_")
#   p = ggplot(all_Pro.Pho,aes(x=correlation, fill=dataset))
#   p = p + geom_histogram(binwidth = 0.05, alpha=0.5, position="identity") + theme_bw()
#   p = p + xlim(-1,1) #+ geom_vline(xintercept = medi, color="red")
#   p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
#   p
#   ggsave(file=fn, height=5, width=8, useDingbats=FALSE)
#   
#   ##### faceted graph that contains every level #####
#   all_CNV.RNA$Comparison = "CNV vs. RNA"
#   all_RNA.PRO$Comparison = "RNA vs. Protein"
#   all_CNV.PRO$Comparison = "CNV vs. Protein"
#   all_Pro.Pho$Comparison = "Protein vs. Phosphoprotein"
#   #colnames(all_RPPA.PRO)[1] = "gene"
#   all_groups = rbind(all_CNV.RNA,all_RNA.PRO,all_CNV.PRO,all_Pro.Pho)
#   all_groups$Comparison = factor(all_groups$Comparison, levels = c("CNV vs. RNA","RNA vs. Protein","CNV vs. Protein","Protein vs. Phosphoprotein"))
#   
#   fn = paste(pd, 'pan3cancer_crossLvl_p_correlation_histogram_U24_wo_legend.pdf',sep ="_")
#   p = ggplot(data = all_groups)
#   p = p + facet_grid(.~Comparison)
#   p = p + geom_density(aes(x=correlation, color=dataset, fill=dataset), alpha=0.2) + theme_bw()
#   p = p + labs(x="Pearson correlation coefficient", y="Density") + scale_x_continuous(breaks=c(-1,0,1), limits=c(-1,1))
#   p = p + theme(legend.position="none")
#   p = p + theme(text = element_text(colour="black", size=8), axis.text.x = element_text(colour="black", size=6), axis.text.y = element_text(colour="black", size=6))
#   p
#   ggsave(file=fn, useDingbats=FALSE, width = 5, height=2)
#   
# }
# 
# 
# ### CNV.RNA ###
# BRCA_CNV.RNA$dataset = "BRCA"
# OV_CNV.RNA$dataset = "OV"
# CRC_CNV.RNA$dataset = "CRC"
# 
# all_CNV.RNA = rbind(BRCA_CNV.RNA,OV_CNV.RNA,CRC_CNV.RNA)
# all_CNV.RNA$dataset = as.factor(all_CNV.RNA$dataset)
# 
# fn = paste(pd, 'pan3can_CNV.RNA_correlation.pdf',sep ="_")
# p = ggplot(all_CNV.RNA,aes(x=correlation, fill=dataset))
# p = p + geom_histogram(binwidth = 0.05, alpha=0.5, position="identity") + theme_bw()
# p = p + xlim(-1,1) #+ geom_vline(xintercept = medi, color="red")
# p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=5, width=8, useDingbats=FALSE)
# 
# ### RNA.PRO ###
# BRCA_RNA.Pro$dataset = "BRCA"
# OV_RNA.JHU_Pro$dataset = "OV_JHU"
# OV_RNA.PNNL_Pro$dataset = "OV_PNNL"
# CRC_RNA.Pro$dataset = "CRC"
# 
# all_RNA.PRO = rbind(BRCA_RNA.Pro,OV_RNA.JHU_Pro,OV_RNA.PNNL_Pro,CRC_RNA.Pro)
# all_RNA.PRO$dataset = as.factor(all_RNA.PRO$dataset)
# 
# fn = paste(pd, 'pan3can_RNA.PRO_correlation.pdf',sep ="_")
# p = ggplot(all_RNA.PRO,aes(x=correlation, fill=dataset))
# p = p + geom_histogram(binwidth = 0.05, alpha=0.5, position="identity") + theme_bw()
# p = p + xlim(-1,1) #+ geom_vline(xintercept = medi, color="red")
# p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=5, width=8, useDingbats=FALSE)
# 
# ### PRO.PHO ###
# BRCA_Pro.Pho$dataset = "BRCA"
# OV_Pro.PNNL_Pho$dataset = "OV_PNNL"
# 
# all_Pro.Pho = rbind(BRCA_Pro.Pho,OV_Pro.PNNL_Pho)
# all_Pro.Pho$dataset = as.factor(all_Pro.Pho$dataset)
# 
# fn = paste(pd, 'pan3can_PRO.PHO_correlation.pdf',sep ="_")
# p = ggplot(all_Pro.Pho,aes(x=correlation, fill=dataset))
# p = p + geom_histogram(binwidth = 0.05, alpha=0.5, position="identity") + theme_bw()
# p = p + xlim(-1,1) #+ geom_vline(xintercept = medi, color="red")
# p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=5, width=8, useDingbats=FALSE)
# 
# ### PRO.RPPA ###
# BRCA_Pro.RPPA$dataset = "BRCA"
# OV_JHU_Pro.RPPA$dataset = "OV_JHU"
# OV_PNNL_Pro.RPPA$dataset = "OV_PNNL"
# CRC_Pro.RPPA$dataset = "CRC"
# 
# all_RPPA.PRO = rbind(BRCA_Pro.RPPA,OV_JHU_Pro.RPPA,OV_PNNL_Pro.RPPA,CRC_Pro.RPPA)
# all_RPPA.PRO$dataset = as.factor(all_RPPA.PRO$dataset)
# 
# fn = paste(pd, 'pan3can_RPPA.PRO_correlation.pdf',sep ="_")
# p = ggplot(all_RPPA.PRO,aes(x=correlation, fill=dataset))
# p = p + geom_histogram(binwidth = 0.05, alpha=0.5, position="identity") + theme_bw()
# p = p + xlim(-1,1) #+ geom_vline(xintercept = medi, color="red")
# p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=5, width=8, useDingbats=FALSE)
# 
# ##### RPPA analysis to see if better correlations are seen in the validated markers #####
# 
# RPPA_info = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/collaborations/TCPA_2015-10-30/TCGA-marker-protein-validation.txt")
# RPPA_info$marker = paste(RPPA_info$Gene.s.,RPPA_info$Protein.Marker.ID,sep="_")
# 
# RPPA_info_merge.brca = merge(RPPA_info,BRCA_Pro.RPPA,by="row.names")
# RPPA_info_merge.ov.jhu = merge(RPPA_info,OV_JHU_Pro.RPPA,by="row.names")
# RPPA_info_merge.ov.pnnl = merge(RPPA_info,OV_PNNL_Pro.RPPA,by="row.names")
# RPPA_info_merge.crc = merge(RPPA_info,CRC_Pro.RPPA,by="row.names")
# 
# RPPA_info_merge.brca$set = "BRCA"
# RPPA_info_merge.ov.jhu$set = "OV JHU"
# RPPA_info_merge.ov.pnnl$set = "OV PNNL"
# RPPA_info_merge.crc$set = "CRC"
# 
# RPPA_info_merge_all = rbind(RPPA_info_merge.brca,RPPA_info_merge.crc,RPPA_info_merge.ov.jhu,RPPA_info_merge.ov.pnnl)
# 
# fn = paste(pd, "pan3can_RPPA.PRO_correlation_antibody_status.pdf", sep="_")
# p = ggplot(data=RPPA_info_merge_all)
# p = p + facet_grid(.~set)
# p = p + geom_violin(alpha=0.5, aes(x=Validation.Status, y=correlation, fill=Validation.Status)) + geom_jitter(aes(x=Validation.Status, y=correlation))  + guides(fill=FALSE)
# p = p + labs( x = "RPPA antibody status", y = "RPPA vs. Pro correlation") + theme_bw()
# p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14))
# p
# ggsave(file=fn)
# 
# ##### faceted graph that contains every level #####
# all_CNV.RNA$Comparison = "CNV vs. RNA"
# all_RNA.PRO$Comparison = "RNA vs. MS Protein"
# all_RPPA.PRO$Comparison = "RPPA vs. MS Protein"
# colnames(all_RPPA.PRO)[1] = "gene"
# all_groups = rbind(all_CNV.RNA,all_RNA.PRO,all_RPPA.PRO)
# 
# fn = paste(pd, 'pan3cancer_crossLvl_p_correlation_histogram.pdf',sep ="_")
# p = ggplot(data = all_groups)
# p = p + facet_grid(.~Comparison)
# p = p + geom_density(aes(x=correlation, color=dataset, fill=dataset), alpha=0.2) + theme_bw()
# p = p + xlim(-1,1) + labs(x="Pearson correlation coefficient")
# p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))
# p
# ggsave(file=fn, useDingbats=FALSE, width = 15)