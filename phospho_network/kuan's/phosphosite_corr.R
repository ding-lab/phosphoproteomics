##### phosphosite_corr.R #####
# Kuan-lin Huang @ WashU 2016 Feb
# cross-correlation between different phosphosites

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/phospho_network/")
source("/Users/khuang/bin/LIB_exp.R")

### function ###
calc_corr = function (x, name="data"){
  
  corr_m = vector("list")
  
  k = 1
  for (i in 1:nrow(x)){
    for (j in (i+1):nrow(x)){
      corr_p = try(cor(unlist(x[i,]),unlist(x[j,]), method = "pearson",use="pairwise.complete.obs"), silent=T)
      if (is(corr_p,"try-error")){ corr_p = NA} 
      corr_m[[k]] = c(row.names(x)[i], row.names(x)[j], corr_p)
      k = k + 1 
    } 
  }
  
  corr_table = do.call(rbind, corr_m)
  
  colnames(corr_table)[1] = "Phosphosite1"
  colnames(corr_table)[2] = "Phosphosite2"
  colnames(corr_table)[3] = "Pearson_R"
  write.table(corr_table, row.names = F, quote=F, sep = '\t', file=paste(pd,name,"phosites_p_correlation.txt",sep="_"))
  
  corr_m_df = as.data.frame(corr_table)
  corr_m_df$correlation = as.numeric(as.character(corr_m_df$correlation))
  medi = round(median(corr_m_df$correlation, na.rm=T),digit=3)
  per_pos = round(sum(corr_m_df$correlation > 0, na.rm=T)/length(corr_m_df$correlation[!is.na(corr_m_df$correlation)]),digit=3)
  text_sum = as.data.frame(matrix(c(-0.75, -0.75, 500, 300,medi,per_pos), ncol=3))
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

### files
brca_file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
ov_file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"

brca_pho = read.table(row.names=1, header=TRUE, sep="\t", file= brca_file)
ov_pho = read.table(row.names=1, header=TRUE, sep="\t", file= ov_file)


brca_pho_10NA = brca_pho[rowSums(!is.na(brca_pho))>=10,]
BRCA_PHO = calc_corr(brca_pho_10NA)

ov_pho_10NA = ov_pho[rowSums(!is.na(ov_pho))>=10,]
OV_PHO = calc_corr(ov_pho_10NA)

### CNV.RNA ###
BRCA_PHO$dataset = "BRCA"
OV_PHO$dataset = "OV"

all_PHO = rbind(BRCA_PHO,OV_PHO)
all_PHO$dataset = as.factor(all_PHO$dataset)

fn = paste(pd, 'brca.ov_PHO_correlation.pdf',sep ="_")
p = ggplot(all_PHO,aes(x=correlation, fill=dataset))
p = p + geom_histogram(binwidth = 0.05, alpha=0.5, position="identity") + theme_bw()
p = p + xlim(-1,1) #+ geom_vline(xintercept = medi, color="red")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
ggsave(file=fn, height=5, width=8, useDingbats=FALSE)