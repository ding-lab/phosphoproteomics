##### find_pathway.R #####
# Kuan-lin Huang @ WashU 2015 Oct - 2017 Feb
# find activated KEGG pathway through z-score of proteome and phosphoproteomes
# plot the protein and phosphoprotein data to the KEGG pathway

base = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/"
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"
setwd(paste(base,"pathway_activation",sep=""))
source("/Users/khuang/bin/LIB_exp.R")
source(paste(base,"pathway_activation/pathway_activation.R",sep=""))
library(KEGGprofile)

# function to format CDAP proteome data processed by Kuan
format_pro = function(Pro.m){
  colnames(Pro.m) = sub(".Unshared.Log.Ratio","",colnames(Pro.m)) 
  #Pro.m.d = Pro.m[Pro.m$Gene %in% druggable,]
  row.names(Pro.m) = make.names(Pro.m$Gene, unique =T)
  Pro.m.dn = Pro.m[,-1]
  Pro.m.dn = smad(as.matrix(Pro.m.dn))
  return(Pro.m.dn)
}

format_crc = function(Pro.m){
  row.names(Pro.m) = make.names(Pro.m$Gene, uniq=T)
  Gene = Pro.m$Gene
  Pro.m = as.matrix(Pro.m[,-1])
  colnames(Pro.m) = sub(".Unshared.Spectral.Counts","",colnames(Pro.m)) 
  
  # quantile normalization using function from limma and log2 transformation: Nature CRC proteogenomics 2014
  Pro.mn = normalizeQuantiles(Pro.m,ties=T)
  Pro.mnl = log2(Pro.mn)
  Pro.m.n = smad(as.matrix(Pro.mnl))
  
  return(Pro.m.n)
}

score_calc_by_gene = function(m){
  m = as.matrix(m)
  m[!is.finite(m)] = NA
  for (i in 1:nrow(m)){
    m[i,] = (m[i,] - median(m[i,], na.rm=T))/iqr(m[i,], na.rm=T)
  }
  return(m)
}

##### Phospho data #####
### BRCA ###
BRCA_Pho = read.table(row.names = 1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt")
#BRCA_Pho.n = format_pro(BRCA_Pho)
BRCA_Pho = as.matrix(BRCA_Pho)
BRCA_Pho_score = score_calc_by_gene(BRCA_Pho)

# calculate mean phosphorylation level in subtypes ------------------------
clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))

BRCA_Pho_score_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(BRCA_Pho_score)),4), ncol=4, byrow=T))
colnames(BRCA_Pho_score_subtype_mean) <- c("Her2","LumA","LumB","Basal")
for (cohort in c("Her2","LumA","LumB","Basal")) {
  subtype_sample <- na.omit(colnames(clinical)[clinical[1,]==cohort])
  BRCA_Pho_score_subtype_mean[,cohort] <- rowMeans(BRCA_Pho_score[,c(subtype_sample)], na.rm = TRUE)
}
row.names(BRCA_Pho_score_subtype_mean) = row.names(BRCA_Pho_score)

# BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted.txt",sep="")) # use not normalized, doesn't matter for Spearman
# BRCA_expressed_genes = row.names(BRCA_RNA)[rowSums(BRCA_RNA,na.rm=T) > 74*500]

trans_result_f = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/manuscript/supplementary_tables/kinase_substrate_regression_trans_edited.txt"
trans_result = read.table( header=TRUE, sep="\t", file = trans_result_f)
trans_result_sig = trans_result[(trans_result$self & trans_result$FDR_pro_kin < 0.05) | (!trans_result$self & trans_result$FDR_pho_kin < 0.05),]
sig_genes = unique(as.character(unique(trans_result_sig$KINASE,trans_result_sig$SUBSTRATE)))

for (pathway in names(KEGG_signaling)){
  hsa = strsplit(pathway, split = "\t")[[1]][1]
  pathway_name=strsplit(pathway, split = "\t")[[1]][2]
  KEGG[[pathway]] = KEGG[[pathway]][KEGG[[pathway]] %in% sig_genes]
}

BRCA_Pho_pathway = KEGGpathway_activation(BRCA_Pho_score_subtype_mean)

# a brief visualization; need to use ggplot2 and facet in the future
BRCA_Pho_pathway_all_setM = matrix(as.numeric(as.character(unlist(BRCA_Pho_pathway$all_setM))), nrow=45, byrow=T)
BRCA_Pho_pathway_all_fdrM = matrix(as.numeric(as.character(unlist(BRCA_Pho_pathway$all_fdr))), nrow=45, byrow=T)
row.names(BRCA_Pho_pathway_all_setM) = row.names(BRCA_Pho_pathway$all_setM)
colnames(BRCA_Pho_pathway_all_setM) = colnames(BRCA_Pho_pathway$all_setM)
row.names(BRCA_Pho_pathway_all_fdrM) = row.names(BRCA_Pho_pathway$all_setM)
colnames(BRCA_Pho_pathway_all_fdrM) = colnames(BRCA_Pho_pathway$all_setM)

pdf(paste(pd,'BRCA_subtypes_Pho_KEGG_signaling_pathway.pdf', sep="_"), useDingbats=FALSE)
par(oma=c(1,3,1,20))

BRCA_Pho_pathway = heatmap.2(BRCA_Pho_pathway_all_setM, trace="none",na.color="white", notecol="black",Rowv=F,
                             cexRow=0.7,cexCol=1, scale="none",#dendrogram='column',
                             col=getPalette, margins=c(5,5)) #
dev.off()

BRCA_Pho_pathway_all_setM_m = melt(BRCA_Pho_pathway_all_setM)
BRCA_Pho_pathway_all_fdrM_m = melt(BRCA_Pho_pathway_all_fdrM)
colnames(BRCA_Pho_pathway_all_setM_m) = c("Pathway","Sample","Global_phosphorylation")
BRCA_Pho_pathway_all_setM_m$Pathway = as.character(BRCA_Pho_pathway_all_setM_m$Pathway)
BRCA_Pho_pathway_all_setM_m$Sample = as.character(BRCA_Pho_pathway_all_setM_m$Sample)
colnames(BRCA_Pho_pathway_all_fdrM_m) = c("Pathway","Sample","FDR")
BRCA_Pho_pathway_all_fdrM_m$Pathway = as.character(BRCA_Pho_pathway_all_fdrM_m$Pathway)
BRCA_Pho_pathway_all_fdrM_m$Sample = as.character(BRCA_Pho_pathway_all_fdrM_m$Sample)
BRCA_Pho_pathway_all_merge = merge(BRCA_Pho_pathway_all_setM_m, BRCA_Pho_pathway_all_fdrM_m, by=c("Pathway","Sample"))
BRCA_Pho_pathway_all_merge$Sig = FALSE
BRCA_Pho_pathway_all_merge$Sig[BRCA_Pho_pathway_all_merge$FDR < 0.05] = TRUE
BRCA_Pho_pathway_all_merge$hsaID = gsub("\t.*","",BRCA_Pho_pathway_all_merge$Pathway)
BRCA_Pho_pathway_all_merge$Pathway = gsub(".*\t","",BRCA_Pho_pathway_all_merge$Pathway)

BRCA_fn = paste(pd,"BRCA_subtype_pathway_activation.tsv",sep="_")
write.table(BRCA_Pho_pathway_all_merge, row.names = F, quote=F, sep = '\t', file=BRCA_fn)

##### plotting #####
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

source("/Users/khuang/bin/LIB_exp.R") # so that KEGG is sourced again

plotWrapH = function (sample, pathway){
  genes = KEGG[[pathway]]
  pathway_gene = BRCA_Pho_score_subtype_mean[row.names(BRCA_Pho_score_subtype_mean) %in% genes,sample,drop=F]
  pathwayID = gsub("\t.*","",pathway)
  pathwayIDnum = gsub("hsa","",pathwayID)
  plotKEGG_pathway(pathway_gene, pathwayIDnum)
  command= paste("mv ",pathwayID,"_profile_bg.png ",pd,"_",sample,"_",pathwayID,"_bg.png", sep="")
  system(command)
}

paths = c("hsa04110\tCell cycle","hsa04630\tJak-STAT signaling pathway","hsa04310\tWnt signaling pathway",
          "hsa04115\tp53 signaling pathway","hsa04150\tmTOR signaling pathway","hsa04151\tPI3K-Akt signaling pathway",
          "hsa04010\tMAPK signaling pathway","hsa04012\tErbB signaling pathway","hsa04014\tRas signaling pathway")
for (subtype in colnames(BRCA_Pho_score_subtype_mean)){
  for (pathway in paths){
    plotWrapH( subtype, pathway)
  # plotWrapH(subtype, "04151") # PI3K-AKT signaling
  # plotWrapH(subtype, "04012") # ERBB2 signaling
  # plotWrapH(subtype, "04064") # NFkB signaling
  # plotWrapH(subtype, "04010") # MAPK signaling
  # plotWrapH(subtype, "hsa04110\tCell cycle") # cell cycle
  # plotWrapH(subtype, "04115") # p53 signaling
  # plotWrapH(subtype, "05200") # pathways in cancer
  }
}