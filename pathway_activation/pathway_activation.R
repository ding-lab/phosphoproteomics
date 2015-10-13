##### pathway_activation.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# outlier analysis pipeline
# called by other scripts
# pathway_activation.R


setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/pathway_activation")
source("/Users/khuang/bin/LIB_exp.R")
#library("reactome.db")
library(KEGGprofile)
library(biomaRt)

# Get list KEGG and REACT for analysis
load("/Users/khuang/bin/2015-08-01_Gene_Set.RData")

KEGG_signaling = KEGG[grep("signaling", names(KEGG))]

##### Can also get other pathways from Broad's GSEA website #####
# test for gene set activation using z-score
KEGGpathway_activation = function (m){
  m.n = deparse(substitute(m))
  # print header
  cat("##### KEGG pathway differential regulation analysis #####\n")
  cat("Detecting differentially expressed KEGG pathways in", m.n, "\n")
  a = factor()
  for (pathway in names(KEGG_signaling)){
    a = c(a,pathway)
  }
  all_fdr = matrix(,nrow=length(a),ncol=0)
  row.names(all_fdr) = a
  all_setM = matrix(,nrow=length(a),ncol=0)
  row.names(all_setM) = a
  
  for (s in colnames(m)){
    sample = m[,s,drop=F]
    sample = sample[!is.na(sample),,drop=F]
    stats = matrix(,nrow=0,ncol=10)
    colnames(stats) = c("hsa","Pathway_name","Num_genes","Num_pathway_genes","All_gene_mean", "Pathway_mean", "T_P", "K-S_P", "Wilcox_P","Up")
    for (pathway in names(KEGG_signaling)){
      hsa = strsplit(pathway, split = "\t")[[1]][1]
      pathway_name=strsplit(pathway, split = "\t")[[1]][2]
      pathway_genes = KEGG[[pathway]]
      inSet = rownames(sample) %in% pathway_genes
      numAllGene = length(inSet)
      numPathwayGene = table(inSet)[2]
      if (is.na(numPathwayGene) || numPathwayGene < 3){next}
      #geneSetP = geneSetTest(inSet, sample)[1]
      setExp = sample[rownames(sample) %in% pathway_genes,]
      ksP = ks.test(x=setExp, y=sample)$p
      WilcoxP = wilcox.test(setExp,sample)$p.value
      T_p = t.test(x=setExp, y=sample)$p.value
      allM = mean(sample)
      setM = mean(setExp)
      #foldChange = log2(setM/allM)
      Up = setM > allM
      a = c(hsa, pathway_name, numAllGene,numPathwayGene,allM, setM, T_p, ksP, WilcoxP, Up)
      stats=rbind(stats,a)
      row.names(stats)=NULL
    }
    
    Wilcox_fdr=p.adjust(stats[,"Wilcox_P"], method="BH")
    stats=as.data.frame(cbind(stats, Wilcox_fdr))
    
    # all samples
    sample_fdr = stats[,"Wilcox_fdr",drop=F] 
    colnames(sample_fdr) = s
    all_fdr = cbind(all_fdr, sample_fdr)
    sample_set = stats[,"Pathway_mean",drop=F]
    colnames(sample_set) = s
    all_setM = cbind(all_setM, sample_set) 
    
    # sample level results    
    stats=stats[order(as.numeric(stats$Wilcox_fdr), stats$Wilcox_P, decreasing=FALSE),]
    tn = paste(pd,m.n,s,"KEGGSig_pathway_activation.txt", sep="_")
    write.table(stats, file=tn, quote=F, sep = '\t', row.names=F)
    cat("Results for", s, "printed to", tn, "\n")
    
  }
  # result for all samples
  # print one table for up and one for down
  tn = paste(pd,m.n,"all_KEGGSig_pathway_activation_FDR.txt", sep="_")
  write.table(all_fdr, file=tn, quote=F, sep = '\t', row.names=T)
  cat("FDR results for", m.n, "printed to", tn, "\n")
  tn = paste(pd,m.n,"all_KEGGSig_pathway_activation_setMean.txt", sep="_")
  write.table(all_setM, file=tn, quote=F, sep = '\t', row.names=T)
  cat("Fold change results for", m.n, "printed to", tn, "\n")
  
  return(list("all_fdr"=all_fdr, "all_setM"=all_setM))
}


Reactome_pathway_activation = function (m){
  m.n = deparse(substitute(m))
  # print header
  cat("##### REACTOME pathway differential regulation analysis #####\n")
  cat("Detecting differentially expressed Reactome pathways in", m.n, "\n")
  
  for (s in colnames(m)){
    sample = m[,s,drop=F]
    sample = sample[!is.na(sample),,drop=F]
    stats = matrix(,,ncol=10)
    colnames(stats) = c("reactID","Pathway name","Num_genes","Num_pathway_genes","All_gene_mean", "Pathway_mean", "T_P", "GeneSet_P", "K-S_P", "Wilcox_P")
    for (pathway in names(REACT)){
      reactID = strsplit(pathway, split = "\t")[[1]][1]
      pathway_name=strsplit(pathway, split = "\t")[[1]][2]
      pathway_genes = REACT[[pathway]]
      inSet = rownames(sample) %in% pathway_genes
      numAllGene = length(inSet)
      numPathwayGene = table(inSet)[2]
      if (is.na(numPathwayGene) || numPathwayGene < 3){next}
      geneSetP = geneSetTest(inSet, sample)[1]
      setExp = sample[rownames(sample) %in% pathway_genes,]
      ksP = ks.test(x=setExp, y=sample)$p
      WilcoxP = wilcox.test(setExp,sample)$p.value
      T_p = t.test(x=setExp, y=sample)$p.value
      allM = mean(sample)
      setM = mean(setExp)
      a = c(reactID, pathway_name, numAllGene,numPathwayGene,allM, setM, T_p, geneSetP, ksP, WilcoxP)
      stats=rbind(stats,a)
      row.names(stats)=NULL
    }
    Wilcox_fdr=p.adjust(stats[,"Wilcox_P"], method="BH")
    stats=as.data.frame(cbind(stats, Wilcox_fdr))
    stats=stats[order(as.numeric(stats$Wilcox_fdr), stats$Wilcox_P, decreasing=FALSE),]
    tn = paste(pd,m.n,s,"REACTpathway_activation.txt", sep="_")
    write.table(stats, file=tn, quote=F, sep = '\t', row.names=F)
    cat("Results for", s, "printed to", tn, "\n")
  }
}

getKeggPathwayGenes<-function(hsaID){
  kegg_pathway=readLines(paste("http://rest.kegg.jp/get/",hsaID,sep=""))
  # KGML: showing the relations: http://rest.kegg.jp/get/hsa04012/kgml 
  kegg_pathway_genes=kegg_pathway[grep(";",kegg_pathway)]
  kegg_pathway_genes=sub(";.*", "",kegg_pathway_genes )
  kegg_pathway_genes=kegg_pathway_genes[-grep("CLASS",kegg_pathway_genes)]
  kegg_pathway_genes[grep("GENE", kegg_pathway_genes)]=gsub("GENE","",kegg_pathway_genes[grep("GENE", kegg_pathway_genes)])
  
  kegg_pathway_genes=gsub("\\s","",gsub("\\s[0-9]+", "", perl = T, kegg_pathway_genes))
  return(kegg_pathway_genes)
  
}

plotKEGG_pathway = function(exp,pathwayNum) {
  exp=exp[rowSums(is.na(exp))==0,,drop=FALSE]
  exp2=rep(0,nrow(exp))
  exp = cbind(exp,exp2)
  # convert hugo gene name rownames to hsaID; bug: this step requires two columns...
  exp = convertId(exp,filters="hgnc_symbol")  
  exp=exp[,1,drop=F]
  limit = max(-min(exp),max(exp))
  col = col_by_value(exp, col = RdBu1024, breaks=seq(-limit,limit,length.out=1025),showColorBar = T)
  
  temp = plot_pathway(exp, type = "bg", bg_col = col, text_col = "black",
                     magnify = 1.2, species = "hsa", database_dir = system.file("extdata", package = "KEGGprofile"),
                     pathway_id = pathwayNum)
}

plotKEGG_pathwayDev = function(exp,pathwayNum) {
  exp=exp[rowSums(is.na(exp))==0,,drop=FALSE]
  exp2=rep(0,nrow(exp))
  exp = cbind(exp,exp2)
  # convert hugo gene name rownames to hsaID; bug: this step requires two columns. That's why the exp2 column is added
  exp = convertId(exp,filters="hgnc_symbol")  
  exp=exp[,-ncol(exp),drop=F]
  limit = max(-min(exp),max(exp))
  col = col_by_value(exp, col = RdBu1024, breaks=seq(-limit,limit,length.out=1025),showColorBar = T)
  
  temp = plot_pathway(exp, type = "bg", bg_col = col, text_col = "black",
                      magnify = 1.2, species = "hsa", database_dir = system.file("extdata", package = "KEGGprofile"),
                      pathway_id = pathwayNum)
}

##### PLOTTING FUNCTION #####
### TODO: add proteomic (scale color?) 
# plot with color scale on all genes
plotWrap = function (sample, pathway){
  whim = ITRAQpho_outlier_zscore.c2[,sample,drop=F]
  plotKEGG_pathway(whim, pathway)
  command= paste("mv hsa",pathway,"_profile_bg.png ",pd,"_",sample,"_",pathway,"_bg.png", sep="")
  system(command)
}
# plot human data with color scale on genes only in the specific pathway
plotWrapH = function (sample, pathway){
  hsaID = paste("hsa",pathway,sep="")
  genes = getKeggPathwayGenes(hsaID)
  whim = BRCA77pho_outlier_zscore.c[,sample,drop=F]
  whim = whim[row.names(whim) %in% genes,,drop=F]
  plotKEGG_pathway(whim, pathway)
  command= paste("mv hsa",pathway,"_profile_bg.png ",pd,"_",sample,"_",pathway,"_bg.png", sep="")
  system(command)
}
# plot WHIM data with color scale on genes only in the specific pathway
plotWrapW = function (sample, pathway){
  hsaID = paste("hsa",pathway,sep="")
  genes = getKeggPathwayGenes(hsaID)
  whim = ITRAQpho_outlier_zscore.c[,sample,drop=F]
  whim = whim[row.names(whim) %in% genes,,drop=F]
  plotKEGG_pathway(whim, pathway)
  command= paste("mv hsa",pathway,"_profile_bg.png ",pd,"_",sample,"_",pathway,"_bg.png", sep="")
  system(command)
}
# plot WHIM data with color scale on genes only in the specific pathway; both proteome and phospho
# the resulting plot looks wrong, check or not use!
plotWrapWDev = function (sample, pathway){
  hsaID = paste("hsa",pathway,sep="")
  genes = getKeggPathwayGenes(hsaID)
  whim_pro = ITRAQ_druggable_pro_zscore[,sample,drop=F]
  whim_pho = ITRAQpho_outlier_zscore.c[,sample,drop=F]
  whim = merge(whim_pro, whim_pho, by = "row.names", all=T)
  row.names(whim) = whim[,1]
  whim = whim[,-1]
  whim = whim[row.names(whim) %in% genes,,drop=F]
  whim[is.na(whim)] = 0 # limitation: set NA to 0 otherwise they don't show up
  plotKEGG_pathwayDev(whim, pathway)
  command= paste("mv hsa",pathway,"_profile_bg.png ",pd,"_",sample,"_",pathway,"_bg.png", sep="")
  system(command)
}





