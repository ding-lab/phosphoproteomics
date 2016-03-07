##### plot_outlier_summary.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run outlier analysis for 3 cancer types and plot the result

##### dependencies #####
setwd("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier")
source("/Users/khuang/bin/LIB_exp.R")

summary = read.table(header=T, sep="\t","/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier/results/2016-03-07_outlier_summary_filtered_mRNA_500.txt")


#YlOrRd = brewer.pal("#FFFFFF","#fed976","#bd0026") 
getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
outlier.colors=c("NA", "#000000")

summary$Level = factor(summary$Level, levels = c("MUT","CNV", "RNA", "RPPA", "PRO", "PHO"))

#levels(summary$Level) = c("CNV", "RNA", "PRO", "PHO", "SMT")
summary$rounded_outlier_percentage = round(as.numeric(summary$Outlier_precentage), digits=1)
summary$Gene = factor(summary$Gene, levels = summary$Gene[order(summary$rounded_outlier_percentage)])

summary$truncated_outlier_percentage = summary$rounded_outlier_percentage
summary[summary$truncated_outlier_percentage>=30,]$truncated_outlier_percentage = 30

fn = paste(pd, 'pan3can_druggable_outliers_mRNA10_summary.pdf',sep ="_")
p = ggplot(data=summary)
p = p + facet_grid(.~Cancer,drop=T,scales = "free", space = "free")
#p = p + coord_equal()
p = p + geom_tile(aes(x=Level, y=Gene, fill=truncated_outlier_percentage), linetype="blank") + scale_fill_gradientn(name= "Percentage", colours=getPalette(100), na.value=NA, limit =c(0,30))
p = p + geom_text(aes(x=Level, y=Gene, label = rounded_outlier_percentage, stringsAsFactors=FALSE), color="black", size=3)
p = p  + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
ggsave(file=fn, height=10, width=15, useDingbats=FALSE)

##### short summary ##### 
if (FALSE){
# note this is a bit problematic at this point as didn't print any values lower than 1.5

YlOrRd = brewer.pal(9, "YlOrRd") 
getPalette = colorRampPalette(YlOrRd)
outlier.colors=c("NA", "#000000")

# showing less genes
summary2 = summary
for (gene in unique(summary2$Gene)){
  if (sum(summary2[summary$Gene == gene,]$Outlier_precentage, na.rm=T) < 10){
    summary2 = summary2[summary2$Gene != gene,]
    print(gene)
  }
}

#summary2$Level = factor(summary$Level, levels = c("CNV", "RNA", "PRO", "PHO", "GLY"))

#levels(summary$Level) = c("CNV", "RNA", "PRO", "PHO", "GLY")
#summary$rounded_outlier_percentage = round(as.numeric(summary$Outlier_precentage), digits=1)
#summary2$Gene = factor(summary2$Gene, levels = summary2$Gene[order(summary2$rounded_outlier_percentage)])
fn = paste(pd, 'pan3can_druggable_outliers_summary_short.pdf',sep ="_")
p = ggplot(data=summary2)
p = p + facet_grid(.~Cancer)
p = p + geom_tile(aes(x=Level, y=Gene, fill=rounded_outlier_percentage), linetype="blank") + scale_fill_gradientn(name= "Percentage", colours=getPalette(100), na.value=NA)#, limit =c(0,22))
p = p + geom_text(aes(x=Level, y=Gene, label = rounded_outlier_percentage, stringsAsFactors=FALSE), color="black", size=3)
p = p  + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=5),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
ggsave(file=fn, height=10, width=15, useDingbats=FALSE)
}

### for grant ###
if (FALSE){
  #markers = c("ERBB2", "IGF1R", "EGFR", "MAPK3", "KRAS", "CDK4", "ARAF","MDM2", "GNAS", "FGF4","AKT3")
  markers = c("ERBB2", "IGF1R", "EGFR", "KRAS", "CDK6", "ALK","MAPK3","PDGFRB","BRAF")
  summary_m = summary[summary$Gene %in% as.character(markers),]
  summary_m = summary_m[summary_m$Cancer %in% c("BRCA","CRC","OV PNNL"),]
  fn = paste(pd, 'pan3can_druggable_outliers_summary_pgdac_selected.pdf',sep ="_")
  p = ggplot(data=summary_m)
  p = p + facet_grid(.~Cancer,drop=T,scales = "free", space = "free")
  #p = p + coord_equal()
  p = p + geom_tile(aes(x=Level, y=Gene, fill=truncated_outlier_percentage), linetype="blank") + scale_fill_gradientn(name= "Percentage", colours=getPalette(100), na.value=NA)#, limit =c(0,30))
  p = p + geom_text(aes(x=Level, y=Gene, label = rounded_outlier_percentage, stringsAsFactors=FALSE), color="black", size=3)
  p = p  + theme_bw() + theme_nogrid() +
    theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=14))
  p
  ggsave(file=fn, height=3, width=9, useDingbats=FALSE)
}

