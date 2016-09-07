#This script takes in a text file which contains gene-correlation coefficient info as table
#Plots the distribution of correlation coefficient and whether the distribution is statistically significantly skewed in
#certain KEGG pathways. 

#USAGE: Rscrip correlation_kegg.R [input text file] [file_description]
# Rscript correlation_kegg.R /Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/cross_correlation/figures/2015-11-11/2015-11-11_KH_BRCA\ RNA\ vs.\ PRO_p_correlation.txt "BRCA_RNA_PRO"
#Input text file: tab delimited where 1st column=gene, 2nd column=correlation coefficient]
#file description: To be added to file name and figure titles. if missing this argument just uses the filePath as its name
library(ggplot2)
library(grid)

args=commandArgs(TRUE)
print(args[1])
filePath=args[1]
if(length(args)>=2){
  figTitle=args[2]
  print(figTitle)
}else{
  figTitle=filePath
  print(figTitle)
}

#Set date and figure location
date=Sys.time()
date = sub(" .*","",date)
figurePath=paste("Figures/", date, "/",sep="")
system(paste("mkdir ", figurePath, sep=""))
date = paste(date, "CJY", sep="_")
#filePath="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/cross_correlation/figures/2015-11-11/2015-11-11_KH_BRCA\ RNA\ vs.\ PRO_p_correlation.txt"
corr=read.table(filePath, sep="\t", header=F, as.is = T); colnames(corr)=c("gene","corr_coeff")
corr$corr_coeff = as.numeric(corr$corr_coeff)

correlation_mean=mean(corr$corr_coeff, na.rm=T)
correlation_median=median(corr$corr_coeff, na.rm=T)
#Draw the histogram of resulting Spearman correlation for all the genes that have observed data for both RSEM RNA and iTRAQ global protein
ggplot(corr, aes(x=corr_coeff))+
  geom_histogram(binwidth = 0.1)+
  geom_density()+
  labs(title = "Spearman Correlation RSEM vs iTRAQ")+
  xlab("Spearman Coefficient")+
  ylab("# of genes")+
  theme_bw()+
  coord_cartesian(xlim=c(-1,1))+
  geom_vline(xintercept = 0, colour="red")+
  geom_vline(xintercept =  correlation_mean)+
  annotate("text", x= correlation_mean, y=200, label=paste("Mean = ", round(correlation_mean,3), sep=""), colour="grey")+
  geom_vline(xintercept =  correlation_median, colour="blue")+
  annotate("text", x= correlation_median+0.2, y=500, label=paste("Median = ", round(correlation_median,3), sep=""), colour="blue")+
  annotate("text", x=0, y=1500, label=paste(round(length(which(corr$corr_coeff>0))/dim(corr)[1]*100, 2), "% positive correlation"), colour="red")

corr$gene<-with(corr,factor(gene,levels = corr[order(corr$corr_coeff, decreasing = T),"gene"] ))

p1<-ggplot(corr, aes(gene,  corr_coeff))+
  geom_bar(aes(gene,  corr_coeff),fill="cyan", stat="identity")+
  theme(axis.text.x=element_blank(), axis.ticks=element_blank())+
  ylab("Spearman Coefficient")+
  coord_cartesian(ylim=c(-1,1))
ggsave(paste(figurePath, figTitle, "_distribution.pdf",sep=""))

rank_diff = vector("list")
plist = vector("list")

pathway_correlation_bar<-function(hsaID, pathway_name){
  detected_genes=intersect(corr$gene, getKeggPathwayGenes(hsaID))
  if(length(detected_genes)>0){
    detection_dataframe=data.frame(detected=rep(0,dim(corr)[1]), gene=corr[order(corr$corr_coeff, decreasing = T),"gene"])
    detection_dataframe[which(detection_dataframe[,"gene"] %in% detected_genes),"detected"]=1
    detection_dataframe$gene=with(detection_dataframe,factor(gene,levels = corr[order(corr$corr_coeff, decreasing = T),"gene"] ))
    
    if (FALSE){
      #     draw a grid for plot layout
      pdf(paste(figurePath, pathway_name, "_kegg_spearman.pdf", sep=""))
      vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(5, 4))) # 5 rows, 4 columns
      
      
      
      p2<-ggplot(data=detection_dataframe, aes(x=gene, y=detected))+
        geom_bar(stat="identity")+  
        theme_bw()+
        theme(axis.text.y=element_text(colour="white"), axis.text.x=element_blank(), axis.ticks=element_blank())+
        xlab(pathway_name)+
        ylab(" ")
      
      print(p1, vp = vplayout(1:4, 1:4))
      print(p2, vp = vplayout(5,1:4))
      dev.off()
    }
    # set this up to also report directionality
    mean_rank_diff = mean(which(detection_dataframe$detected==1), na.rm=T) - mean(1:dim(detection_dataframe)[1])
    p = ks.test(x=which(detection_dataframe$detected==1), y=1:dim(detection_dataframe)[1])$p.value
    result_str = cat(mean_rank_diff, p)
    
    rank_diff[[hsaID]] = mean_rank_diff
    plist[[hsaID]] = p
    return(p)
  }
  else{
    return(NA)
  }
}


getKeggPathwayGenes<-function(hsaID){
  kegg_pathway=readLines(paste("http://rest.kegg.jp/get/",hsaID,sep=""))
  kegg_pathway_genes=kegg_pathway[grep(";",kegg_pathway)]
  kegg_pathway_genes=sub(";.*", "",kegg_pathway_genes )
  kegg_pathway_genes=kegg_pathway_genes[-grep("CLASS",kegg_pathway_genes)]
  kegg_pathway_genes[grep("GENE", kegg_pathway_genes)]=gsub("GENE","",kegg_pathway_genes[grep("GENE", kegg_pathway_genes)])
  
  kegg_pathway_genes=gsub("\\s","",gsub("\\s[0-9]+", "", perl = T, kegg_pathway_genes))
  return(kegg_pathway_genes)
}

#list of all KEGG human pathways
all_kegg_hsa<-readLines("http://rest.kegg.jp/list/pathway/hsa")
all_kegg_hsa<-gsub("path:","",all_kegg_hsa)
all_kegg_hsa<-gsub(" - Homo sapiens \\(human\\)", "", all_kegg_hsa)
all_kegg_hsa<-gsub("\\/", "-", all_kegg_hsa)

ks_test_keggPathway<-sapply(all_kegg_hsa, function(x){
  hsaID=strsplit(x, split = "\t")[[1]][1]
  pathway_name=strsplit(x, split = "\t")[[1]][2]
  print(hsaID)
  print(pathway_name)
  pathway_correlation_bar(hsaID=hsaID, pathway_name=pathway_name)
})

# rank_diff_m = do.call(rbind,rank_diff)
# plist_m = do.call(rbind,plist)

#show the result of KS test in order of increasing P-value. 
ks_test_keggPathway_ordered=ks_test_keggPathway[order(ks_test_keggPathway)]
#adjust the p-values for multiple hypothesis correction. BH: Benjamini Hochberg 
p.adjust(ks_test_keggPathway_ordered, method = "BH")
write.table(ks_test_keggPathway_ordered, file=paste(figTitle, "_ks_test.txt",sep=""), col.names = F, quote=F, append = F, sep = "\t")

