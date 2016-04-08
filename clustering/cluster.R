# cluster.R by Kuan Huang @ WashU 201506
# cluster WHIMs based on (1) iTRAQ (2) LFQ, (1) proteome (2) phosphoproteome 
# cluster WHIMs based on RNA-Seq, CNV

# dependencies
library(ggplot2)
library(reshape)
library(RColorBrewer)
library("gplots")
library(matrixStats)

# others
setwd("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/clustering/")
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/"
source("/Users/khuang/bin/LIB_exp.R")

## function

get.clinical.scale = function() {
  # Set1 colors
  #colors = c(NA, "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999") # positive is red
  # colors = c(NA, "#636363", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999") #positive is dark grey
  # use Perou's intrinsic subtype colors instead
  colors = c(NA, "#636363", "red", "pink", "#ffeda0", "#deebf7", "#FFFF33", "#A65628", "#F781BF", "#999999") #positive is dark grey       
  color.names = c("negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB")
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="status", values=colors)
  
  return(clinical.color.scale)
}

# plot_clin: plot clinical data based on order from the corresponding heatmap2 clustering plot
plot_clin = function(matrix, clus_order, clin, figure){
  order = as.data.frame(colnames(matrix[,clus_order]))
  colnames(order) = colnames(clin)[1]
  # add missing samples to clin through merge
  m = merge(order,clin, by = colnames(clin)[1], all.x = TRUE)
  row.names(m)=m[,colnames(clin)[1]]
  # reorder and melt for hmting
  m = m[colnames(matrix[,clus_order]),]
  m.m = melt(m, id=colnames(clin)[1])
  m.m$ExternalIdentifierName<-with(m.m,factor(ExternalIdentifierName,levels = colnames(matrix[,clus_order])))
  color.scale = get.clinical.scale()
  
  colourCount=length(unique(m.m$value))
  p = ggplot(data=m.m, aes(x=ExternalIdentifierName, y=variable)) + geom_tile(aes(fill = value)) + 
    scale_fill_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(colourCount)) +
    xlab("Sample") + ylab("")  + theme_bw() + guides(fill=guide_legend(title="status",nrow=2)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=14), axis.text.y = element_text(colour="black", size=14), legend.position="bottom")
  p = p + color.scale
  p
  ggsave(file=figure, useDingbats=FALSE)
} 

##### pan3can #####
pan3can_pro = read.table(sep = '\t', header=T, row.names=1,file=paste(baseD,"pan3can/pan3can_PRO_formatted_normalized.txt",sep=""))
pan3can_clin = read.table(sep = '\t', header=T, row.names=1,file=paste(baseD,"pan3can/pan3can_ctype.txt",sep=""))
brca_clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file=paste(baseD,"BRCA/BRCA_clinical_summary.txt",sep=""))

clin = t(pan3can_clin)
BRCA=as.vector(row.names(clin[clin=="BRCA",,drop=F]))
OV=as.vector(row.names(clin[clin=="OV",,drop=F]))
CRC=as.vector(row.names(clin[clin=="CRC",,drop=F]))

basal = as.vector(colnames(brca_clin)[brca_clin["pam50",]=="Basal"])

samples=colnames(pan3can_pro)
samples[samples %in% basal]="green"
samples[samples %in% BRCA]="forestgreen"
samples[samples %in% OV]="orange"
samples[samples %in% CRC]="purple"
ITRAQ_m = as.matrix(pan3can_pro)

q99=quantile(ITRAQ_m, probs=0.99, na.rm=T)
q1=quantile(ITRAQ_m, probs=0.01, na.rm=T)
ITRAQ_m2 = matrix(,nrow=dim(ITRAQ_m)[1],ncol=dim(ITRAQ_m)[2])
colnames(ITRAQ_m2)=colnames(ITRAQ_m)
row.names(ITRAQ_m2)=row.names(ITRAQ_m)
for (i in 1:nrow(ITRAQ_m)){
  if ( (sum(ITRAQ_m[i,][!is.na(ITRAQ_m[i,])] > q99) + sum(ITRAQ_m[i,][!is.na(ITRAQ_m[i,])] < q1)) < 1){
    ITRAQ_m2[i,]=ITRAQ_m[i,]
  }
}
ITRAQ_m3 = ITRAQ_m2[rowSums(!is.na(ITRAQ_m2)) >= 320,] # 1305 observations
SD=rowSds(ITRAQ_m3, na.rm=TRUE)
ITRAQ_m4 = ITRAQ_m3[SD>0.7,] # 458 observations

pdf(paste(pd,'pan3can_PRO_naMax10_SD0.7.pdf', sep="_"),width=30)
par(oma=c(3,5,3,5))
ITRAQ_m4_hm = heatmap.2(ITRAQ_m4, trace="none", na.color="white", notecol="black",
                        cexRow=0.8,cexCol=0.8, scale="none",dendrogram='column', ColSideColors = samples,
                        labRow=NA,labCol=NA,col=getPalette, margins=c(5,5))

par(lend = 1)  
legend("bottomleft",    # location of the legend on the heatmap plot
       legend = c("Basal BRCA", "BRCA", "OV", "CRC"), # category labels
       col = c("green","forestgreen", "orange","purple"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)

dev.off()

## plot corresponding clinical panel 
# figure = paste(pd,'ITRAQ_clin_proteome-ratio-norm_naMax10_SD2.15.pdf', sep="_")
# plot_clin(ITRAQ_m2, clus_order1, clin, figure)


##### check OV #####
PRO = read.table(sep = '\t', header=T, row.names=1,file=paste(baseD,"OV/OV_merged_PRO_formatted_normalized.txt",sep=""))
#clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')
# #row.names(clin)=clin$ExternalIdentifierName
# clin=clin[,-9]
# #clin = t(clin[,1:5])
# #write.table(clin, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified_cleaned.txt', quote=F, sep = '\t')
# 
# Basal=as.vector(clin[clin$Intrinsic.subtype=="Basal",1])
# LumB=as.vector(clin[clin$Intrinsic.subtype=="LumB",1])
# Her2=as.vector(clin[clin$Intrinsic.subtype=="HER2-E",1])
# CLDN_low=as.vector(clin[clin$Intrinsic.subtype=="CLDN low",1])
# 
# ### process and cluster ITRAQ proteome
# row.names(ITRAQ) = ITRAQ$Description
# colnames(ITRAQ) = sub("\\..*", "", colnames(ITRAQ))
# ITRAQ_m = as.matrix(ITRAQ[,-c(1,2,3)])
# ITRAQ_m = ITRAQ_m[,-c(17,18,20)]
ITRAQ_m = as.matrix(PRO)

q99=quantile(ITRAQ_m, probs=0.99, na.rm=T)
q1=quantile(ITRAQ_m, probs=0.01, na.rm=T)
ITRAQ_m2 = matrix(,nrow=dim(ITRAQ_m)[1],ncol=dim(ITRAQ_m)[2])
colnames(ITRAQ_m2)=colnames(ITRAQ_m)
row.names(ITRAQ_m2)=row.names(ITRAQ_m)
for (i in 1:nrow(ITRAQ_m)){
   if ( (sum(ITRAQ_m[i,][!is.na(ITRAQ_m[i,])] > q99) + sum(ITRAQ_m[i,][!is.na(ITRAQ_m[i,])] < q1)) < 1){
     ITRAQ_m2[i,]=ITRAQ_m[i,]
   }
}
ITRAQ_m2 = ITRAQ_m2[rowSums(is.na(ITRAQ_m2)) <= 10,] # 2364 observations
SD=rowSds(ITRAQ_m2, na.rm=TRUE)
ITRAQ_m3 = ITRAQ_m2[SD>0.7,] # 706 observations

pdf(paste(pd,'ov_merged_naMax10_SD0.7.pdf', sep="_"),width=30)
par(oma=c(3,5,3,5))
ITRAQ_m3_hm = heatmap.2(ITRAQ_m3, trace="none", na.color="white", notecol="black",
                      cexRow=0.8,cexCol=0.8, scale="none",dendrogram='column', 
                      labRow=NA,col=getPalette, margins=c(5,5))
#clus_order1 = ITRAQ_m2_hm$colInd
dev.off()

## plot corresponding clinical panel 
# figure = paste(pd,'ITRAQ_clin_proteome-ratio-norm_naMax10_SD2.15.pdf', sep="_")
# plot_clin(ITRAQ_m2, clus_order1, clin, figure)
