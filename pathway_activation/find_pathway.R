##### plot_outlier_summary.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run pathway analysis for 3 cancer types and plot the result
# find activated KEGG pathway through z-score of proteome and phosphoproteomes
# plot the protein and phosphoprotein data to the KEGG pathway

base = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/"
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"
setwd(paste(base,"pathway_activation",sep=""))
source("/Users/khuang/bin/LIB_exp.R")
source(paste(base,"pathway_activation/pathway_activation.R",sep=""))
source_date = "2015-10-09/2015-10-09"

rbind_gtable <- function(...){
  
  gtl = lapply(list(...), ggplotGrob)
  
  bind2 <- function (x, y) 
  {
    stopifnot(ncol(x) == ncol(y))
    if (nrow(x) == 0) 
      return(y)
    if (nrow(y) == 0) 
      return(x)
    y$layout$t <- y$layout$t + nrow(x)
    y$layout$b <- y$layout$b + nrow(x)
    x$layout <- rbind(x$layout, y$layout)
    x$heights <- gtable:::insert.unit(x$heights, y$heights)
    x$rownames <- c(x$rownames, y$rownames)
    x$widths <- grid::unit.pmax(x$widths, y$widths)
    x$grobs <- append(x$grobs, y$grobs)
    x
  }
  
  Reduce(bind2, gtl)
}
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

##### Phospho data #####
### BRCA ###
BRCA_Pho = read.table(row.names = 1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt")
#BRCA_Pho.n = format_pro(BRCA_Pho)
BRCA_Pho = as.matrix(BRCA_Pho)
BRCA_Pho_outliers=find_outlier(BRCA_Pho, "BRCA_PHO")#, filter = F, plot=F)
BRCA_Pho_outliers_zscore=BRCA_Pho_outliers$outlier_zscore
BRCA_Pho_pathway = KEGGpathway_activation(BRCA_Pho_outliers_zscore)

# a brief visualization; need to use ggplot2 and facet in the future
BRCA_Pho_pathway_all_setM = matrix(as.numeric(as.character(unlist(BRCA_Pho_pathway$all_setM))), nrow=44, byrow=T)
BRCA_Pho_pathway_all_fdrM = matrix(as.numeric(as.character(unlist(BRCA_Pho_pathway$all_fdr))), nrow=44, byrow=T)
row.names(BRCA_Pho_pathway_all_setM) = row.names(BRCA_Pho_pathway$all_setM)
colnames(BRCA_Pho_pathway_all_setM) = colnames(BRCA_Pho_pathway$all_setM)
row.names(BRCA_Pho_pathway_all_fdrM) = row.names(BRCA_Pho_pathway$all_setM)
colnames(BRCA_Pho_pathway_all_fdrM) = colnames(BRCA_Pho_pathway$all_setM)

# pdf(paste(pd,'BRCA_Pho_KEGG_signaling_pathway.pdf', sep="_"), width=15,height=10, useDingbats=FALSE)
# par(oma=c(1,3,1,20))
# 
# BRCA_Pho_pathway = heatmap.2(BRCA_Pho_pathway_all_setM, trace="none",na.color="white", notecol="black",Rowv=F,
#                              cexRow=1.2,cexCol=0.6, scale="none",#dendrogram='column',
#                              col=getPalette, margins=c(5,5)) #
# dev.off()

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

plots = list()

BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep=""))
brcaGenes = c("TP53", "ERBB2", "PIK3CA", "CDH1", "GATA3", "MAP3K1")
BRCA_mut_g = as.matrix(BRCA_mut[row.names(BRCA_mut) %in% brcaGenes,])
BRCA_mut_g_m = melt(BRCA_mut_g)

colnames(BRCA_mut_g_m) = c("Gene","Sample","Status")
BRCA_mut_g_m$Status = as.character(BRCA_mut_g_m$Status)
BRCA_mut_g_m$Status[BRCA_mut_g_m$Status %in% c("wt","silent","intronic")] = "wt"
BRCA_mut_g_m$Status[!(BRCA_mut_g_m$Status %in% c("wt","missense"))] = "truncation"


p = ggplot(data=BRCA_mut_g_m)
p = p + geom_tile(aes(x=Sample, y=Gene, fill=Status))
p = p + scale_fill_manual(values=c("blue","green",NA))
p = p + theme_bw() + theme_nogrid()
p = p + theme(text = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              strip.text = element_text(size = 8),
              panel.background = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position="right")
plots[[1]] = p

p = ggplot(data=BRCA_Pho_pathway_all_merge)
p = p + geom_point(aes(x=Sample, y=Pathway, fill=as.numeric(Global_phosphorylation), size=-log10(FDR), color=ifelse(Sig, "black",NA)),pch=21) 
p = p + scale_fill_gradientn(name= "Global Phosphorylation", na.value=NA, colours=RdBu1024, limit=c(-1.5,1.5))
p = p + scale_colour_manual(values=c("black",NA))
p = p + theme_bw() + theme_nogrid()
p = p + theme(text = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              strip.text = element_text(size = 8),
              panel.background = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position="right")
plots[[2]] = p



gp = do.call(rbind_gtable, plots)
# print the integrated plot
grid.newpage()
cal_width = 20
fn = paste(pd, 'merged_mut_pathway.pdf',sep ="_")
pdf(fn, height=15, width=20,useDingbats = F)
grid.draw(gp)
dev.off()

### OV PNNL ###
OV_PNNL_Pho = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt")
#OV_PNNL_Pho.n = format_pro(OV_PNNL_Pho)
OV_PNNL_Pho_pathway = KEGGpathway_activation(OV_PNNL_Pho)

# a brief visualization; need to use ggplot2 and facet in the future
OV_PNNL_Pho_pathway_all_setM = matrix(as.numeric(as.character(unlist(OV_PNNL_Pho_pathway$all_setM))), nrow=44, byrow=T)
row.names(OV_PNNL_Pho_pathway_all_setM) = row.names(OV_PNNL_Pho_pathway$all_setM)
colnames(OV_PNNL_Pho_pathway_all_setM) = colnames(OV_PNNL_Pho_pathway$all_setM)

pdf(paste(pd,'OV_PNNL_Pho_KEGG_signaling_pathway.pdf', sep="_"), width=15,height=10, useDingbats=FALSE)
par(oma=c(1,3,1,20))

OV_PNNL_Pho_pathway = heatmap.2(OV_PNNL_Pho_pathway_all_setM, trace="none",na.color="white", notecol="black",Rowv=F,
                             cexRow=1.2,cexCol=0.6, scale="none",#dendrogram='column',
                             col=getPalette, margins=c(5,5)) #
dev.off()

# ##### PROTEOME DATA #####
# # read z-score table directly
# #fn = paste(base, "druggable_outlier/figures/",source_date,"_KH_BRCA\ druggable\ proteome_outlier_score_table.txt", sep="")
# #PRO = read.table(row.names=1,header=TRUE, sep="\t",file=fn)
# BRCA_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/BRCA/TCGA_Breast_BI_Proteome_CDAP.r2/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
# BRCA_Pro.n = format_pro(BRCA_Pro)
# BRCA_Pro_pathway = KEGGpathway_activation(BRCA_Pro.n)
# # how many are significant; need to coordinate with up or down...
# # BRCA_Pro_pathway_all_fdr = data.frame(matrix(unlist(BRCA_Pro_pathway$all_fdr), nrow=44, byrow=T),stringsAsFactors=FALSE)
# # row.names(BRCA_Pro_pathway_all_fdr) = row.names(BRCA_Pro_pathway$all_fdr)
# # colnames(BRCA_Pro_pathway_all_fdr) = colnames(BRCA_Pro_pathway$all_fdr)
# # for (i in 1:nrow(BRCA_Pro_pathway_all_fdr)){
# #   cat(row.names(BRCA_Pro_pathway_all_fdr)[i], ": ", sum(BRCA_Pro_pathway_all_fdr[i,] <=0.05, na.rm=T), "\n")
# # }
# 
# # a brief visualization; need to use ggplot2 and facet in the future
# BRCA_Pro_pathway_all_setM = matrix(as.numeric(as.character(unlist(BRCA_Pro_pathway$all_setM))), nrow=44, byrow=T)
# row.names(BRCA_Pro_pathway_all_setM) = row.names(BRCA_Pro_pathway$all_setM)
# colnames(BRCA_Pro_pathway_all_setM) = colnames(BRCA_Pro_pathway$all_setM)
# 
# pdf(paste(pd,'BRCA_Pro_KEGG_signaling_pathway.pdf', sep="_"), width=15,height=10, useDingbats=FALSE)
# par(oma=c(1,3,1,20))
# 
# BRCA_Pro_pathway = heatmap.2(BRCA_Pro_pathway_all_setM, trace="none",na.color="white", notecol="black",Rowv=F,
#                           cexRow=1.2,cexCol=0.6, scale="none",#dendrogram='column',
#                           col=getPalette, margins=c(5,5)) #
# dev.off()
# 
# ### OV_PNNL ###
# OV_PNNL_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/PNNL_Proteome_CDAP.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
# OV_PNNL_Pro.n = format_pro(OV_PNNL_Pro)
# OV_PNNL_Pro_pathway = KEGGpathway_activation(OV_PNNL_Pro.n)
# # how many are significant; need to coordinate with up or down...
# OV_PNNL_Pro_pathway_all_fdr = data.frame(matrix(unlist(OV_PNNL_Pro_pathway$all_fdr), nrow=44, byrow=T),stringsAsFactors=FALSE)
# row.names(OV_PNNL_Pro_pathway_all_fdr) = row.names(OV_PNNL_Pro_pathway$all_fdr)
# colnames(OV_PNNL_Pro_pathway_all_fdr) = colnames(OV_PNNL_Pro_pathway$all_fdr)
# for (i in 1:nrow(OV_PNNL_Pro_pathway_all_fdr)){
#   cat(row.names(OV_PNNL_Pro_pathway_all_fdr)[i], ": ", sum(OV_PNNL_Pro_pathway_all_fdr[i,] <=0.05, na.rm=T), "\n")
# }
# 
# # a brief visualization; need to use ggplot2 and facet in the future
# OV_PNNL_Pro_pathway_all_setM = matrix(as.numeric(as.character(unlist(OV_PNNL_Pro_pathway$all_setM))), nrow=44, byrow=T)
# row.names(OV_PNNL_Pro_pathway_all_setM) = row.names(OV_PNNL_Pro_pathway$all_setM)
# colnames(OV_PNNL_Pro_pathway_all_setM) = colnames(OV_PNNL_Pro_pathway$all_setM)
# 
# pdf(paste(pd,'OV_PNNL_Pro_KEGG_signaling_pathway.pdf', sep="_"), width=15,height=10, useDingbats=FALSE)
# par(oma=c(1,3,1,20))
# 
# OV_PNNL_Pro_pathway = heatmap.2(OV_PNNL_Pro_pathway_all_setM, trace="none",na.color="white", notecol="black",Rowv=F,
#                              cexRow=1.2,cexCol=0.6, scale="none",#dendrogram='column',
#                              col=getPalette, margins=c(5,5)) #
# dev.off()
# 
# # ### OV_JHU ###
# # OV_JHU_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/OV/JHU_Proteome_CDAP.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv_hugoified.unshared_log_ratio.txt")
# # OV_JHU_Pro.n = format_pro(OV_JHU_Pro)
# # OV_JHU_Pro_pathway = KEGGpathway_activation(OV_JHU_Pro.n)
# # # how many are significant; need to coordinate with up or down...
# # OV_JHU_Pro_pathway_all_fdr = data.frame(matrix(unlist(OV_JHU_Pro_pathway$all_fdr), nrow=44, byrow=T),stringsAsFactors=FALSE)
# # row.names(OV_JHU_Pro_pathway_all_fdr) = row.names(OV_JHU_Pro_pathway$all_fdr)
# # colnames(OV_JHU_Pro_pathway_all_fdr) = colnames(OV_JHU_Pro_pathway$all_fdr)
# # for (i in 1:nrow(OV_JHU_Pro_pathway_all_fdr)){
# #   cat(row.names(OV_JHU_Pro_pathway_all_fdr)[i], ": ", sum(OV_JHU_Pro_pathway_all_fdr[i,] <=0.05, na.rm=T), "\n")
# # }
# # 
# # # a brief visualization; need to use ggplot2 and facet in the future
# # OV_JHU_Pro_pathway_all_setM = matrix(as.numeric(as.character(unlist(OV_JHU_Pro_pathway$all_setM))), nrow=44, byrow=T)
# # row.names(OV_JHU_Pro_pathway_all_setM) = row.names(OV_JHU_Pro_pathway$all_setM)
# # colnames(OV_JHU_Pro_pathway_all_setM) = colnames(OV_JHU_Pro_pathway$all_setM)
# # 
# # pdf(paste(pd,'OV_JHU_Pro_KEGG_signaling_pathway.pdf', sep="_"), width=15,height=10, useDingbats=FALSE)
# # par(oma=c(1,3,1,20))
# # 
# # OV_JHU_Pro_pathway = heatmap.2(OV_JHU_Pro_pathway_all_setM, trace="none",na.color="white", notecol="black",Rowv=F,
# #                                 cexRow=1.2,cexCol=0.6, scale="none",#dendrogram='column',
# #                                 col=getPalette, margins=c(5,5)) #
# # dev.off()
# # 
# # ##### doesn't work yet, rough #####
# # ### CRC ###
# # CRC_Pro = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/201507_pancan_proteome_CDAP_r2/CRC/VU_Proteome_CDAP.r2/TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv_hugoified.unshared_spectal_counts.txt")
# # CRC_Pro.n = format_crc(CRC_Pro)
# # CRC_Pro_pathway = KEGGpathway_activation(CRC_Pro.n)
# # # how many are significant; need to coordinate with up or down...
# # CRC_Pro_pathway_all_fdr = data.frame(matrix(unlist(CRC_Pro_pathway$all_fdr), nrow=44, byrow=T),stringsAsFactors=FALSE)
# # row.names(CRC_Pro_pathway_all_fdr) = row.names(CRC_Pro_pathway$all_fdr)
# # colnames(CRC_Pro_pathway_all_fdr) = colnames(CRC_Pro_pathway$all_fdr)
# # for (i in 1:nrow(CRC_Pro_pathway_all_fdr)){
# #   cat(row.names(CRC_Pro_pathway_all_fdr)[i], ": ", sum(CRC_Pro_pathway_all_fdr[i,] <=0.05, na.rm=T), "\n")
# # }
# # 
# # # a brief visualization; need to use ggplot2 and facet in the future
# # CRC_Pro_pathway_all_setM = matrix(as.numeric(as.character(unlist(CRC_Pro_pathway$all_setM))), nrow=44, byrow=T)
# # row.names(CRC_Pro_pathway_all_setM) = row.names(CRC_Pro_pathway$all_setM)
# # colnames(CRC_Pro_pathway_all_setM) = colnames(CRC_Pro_pathway$all_setM)
# # 
# # pdf(paste(pd,'CRC_Pro_KEGG_signaling_pathway.pdf', sep="_"), width=15,height=10, useDingbats=FALSE)
# # par(oma=c(1,3,1,20))
# # 
# # CRC_Pro_pathway = heatmap.2(CRC_Pro_pathway_all_setM, trace="none",na.color="white", notecol="black",Rowv=F,
# #                                cexRow=1.2,cexCol=0.6, scale="none",#dendrogram='column',
# #                                col=getPalette, margins=c(5,5)) #
# # dev.off()
# # 
# # 
# # #Reactome_pathway_activation(ITRAQ_outlier_zscore)
# # 
# # # ##### PHOHSPHO DATA ##### 
# # # ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
# # # row.names(ITRAQpho) = ITRAQpho$gene.site
# # # colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
# # # # 56651 phosphosites
# # # ITRAQpho=ITRAQpho[,-c(1,2)]
# # # # get rid of TaxIR, HumIR, WHIM13.1
# # # ITRAQpho = ITRAQpho[,-c(17,18,20)]
# # # ITRAQpho.na = ITRAQpho[rowSums(is.na(ITRAQpho)) <= 14,] #35838 phosphosites
# # # rm(ITRAQpho)
# # # row.names(ITRAQpho.na) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho.na)), unique=T)
# # # row.names(ITRAQpho.na) = make.names(sub(" _.*","",row.names(ITRAQpho.na)), unique=T)
# # # row.names(ITRAQpho.na) = make.names(sub("_.*","",row.names(ITRAQpho.na)), unique=T)
# # # ITRAQpho.na = as.matrix(ITRAQpho.na)
# # # # use the adjusted z-score to find outliers
# # # ITRAQpho_outliers=find_outlier(ITRAQpho.na, "ITRAQ phosphoproteome")#, filter = F, plot=F)
# # # ITRAQpho_outlier_zscore=ITRAQpho_outliers$outlier_zscore
# # # genes = sub("\\..*", "", row.names(ITRAQpho_outlier_zscore))
# # # ITRAQpho_outlier_zscore.c = collapseRows(ITRAQpho_outlier_zscore, rowGroup=genes, rowID=row.names(ITRAQpho_outlier_zscore))$datETcollapsed
# # # KEGGpathway_activation(ITRAQpho_outlier_zscore.c)
# # # #Reactome_pathway_activation(ITRAQpho_outlier_zscore.c)
# # 
# # ##### plotting for WHIM #####
# # # for (whim in colnames(ITRAQpho_outlier_zscore.c)){
# # # #   plotWrapW(whim, "03460") # Fanconi anemia
# # # #   plotWrapW(whim, "04014") # Ras signaling
# # #   plotWrapW(whim, "04915") # ER signaling
# # # #   plotWrapW(whim, "04151") # PI3K-AKT signaling
# # # #   plotWrapW(whim, "04012") # ERBB2 signaling
# # # #   plotWrapW(whim, "04064") # NFkB signaling
# # # #   plotWrapW(whim, "04010") # MAPK signaling
# # # #   plotWrapW(whim, "04110") # cell cycle
# # # #   plotWrapW(whim, "04115") # p53 signaling
# # # #   plotWrapW(whim, "05200") # pathways in cancer
# # # }
# # 
# # # ##### BRCA77 phosphoproteome #####
# # # BRCA77pho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
# # # # collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# # # # 33239 phosphosites
# # # row.names(BRCA77pho) = BRCA77pho$Gene.site
# # # BRCA77pho = BRCA77pho[,-c(1,2)]
# # # BRCA77pho.na = BRCA77pho[rowSums(is.na(BRCA77pho)) <= 67,] #33239 phosphosites
# # # rm(BRCA77pho)
# # # row.names(BRCA77pho.na) = make.names(sub("-NP_\\d+_"," ",row.names(BRCA77pho.na)), unique=T)
# # # row.names(BRCA77pho.na) = make.names(sub(" _.*","",row.names(BRCA77pho.na)), unique=T)
# # # row.names(BRCA77pho.na) = make.names(sub("_.*","",row.names(BRCA77pho.na)), unique=T)
# # # BRCA77pho.m = as.matrix(BRCA77pho.na)
# # # BRCA77pho_outliers=find_outlier(BRCA77pho.m, "BRCA ITRAQ phosphoproteomes")
# # # BRCA77pho_outlier_zscore=BRCA77pho_outliers$outlier_zscore
# # # BRCA77pho_outlier=BRCA77pho_outliers$outlier
# # # genes = sub("\\..*", "", row.names(BRCA77pho_outlier_zscore))
# # # BRCA77pho_outlier_zscore.c = collapseRows(BRCA77pho_outlier_zscore, rowGroup=genes, rowID=row.names(BRCA77pho_outlier_zscore))$datETcollapsed
# # # KEGGpathway_activation(BRCA77pho_outlier_zscore.c)
# # # #Reactome_pathway_activation(BRCA77pho_outlier_zscore.c)
# # 
# # # for (h in colnames(BRCA77pho_outlier_zscore.c)){
# # #   plotWrapH(h, "04151") # PI3K-AKT signaling
# # #   plotWrapH(h, "04012") # ERBB2 signaling
# # #   plotWrapH(h, "04064") # NFkB signaling
# # #   plotWrapH(h, "04010") # MAPK signaling
# # #   plotWrapH(h, "04110") # cell cycle
# # #   plotWrapH(h, "04115") # p53 signaling
# # #   plotWrapH(h, "05200") # pathways in cancer
# # # }
# # 
# # # plot violin for top pathways?
# # #       ### plot dens
# # #       p = ggplot(data=sample2, aes(x=factor(inSet), y =S))
# # #       p = p + geom_violin(aes(fill=factor(inSet)), alpha=0.5) + geom_jitter(height=0, alpha=0.2)
# # #       p = p + theme_bw() +ylab("normalized protein expression ratio")
# # #       p
