##### plot_pathway.R #####
# Kuan-lin Huang @ WashU 2017 Feb
# find activated KEGG pathway through z-score of proteome and phosphoproteomes
# plot the protein and phosphoprotein data to the KEGG pathway

base = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/"
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"
setwd(paste(base,"pathway_activation",sep=""))
source("/Users/khuang/bin/LIB_exp.R")
source(paste(base,"pathway_activation/pathway_activation.R",sep=""))
source_date = "figures/2017-02-15/2017-02-15_KH"

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

score_calc_by_gene = function(m){
  m = as.matrix(m)
  m[!is.finite(m)] = NA
  for (i in 1:nrow(m)){
    m[i,] = (m[i,] - median(m[i,], na.rm=T))/iqr(m[i,], na.rm=T)
  }
  return(m)
}

### BRCA ###
BRCA_Pho = read.table(row.names = 1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt")
BRCA_Pho = as.matrix(BRCA_Pho)
BRCA_Pho_score = score_calc_by_gene(BRCA_Pho)

BRCA_fn = paste(source_date,"BRCA_cross_pathway_activation.tsv",sep="_")
BRCA_Pho_pathway_all_merge = read.table(header = T, quote="", sep = '\t', file=BRCA_fn)

##### plotting #####

plots = list()
brcaGenes = c("TP53", "ERBB2", "PIK3CA","PIK3CB", "MTOR","CDK1", "CDK2", "MAP3K1","MAPK3","MAPK9","MAPK1","MAP3K3")
BRCA_Pho_score_g = BRCA_Pho_score[row.names(BRCA_Pho_score) %in% brcaGenes,]
BRCA_Pho_score_g_m = melt(BRCA_Pho_score_g)
colnames(BRCA_Pho_score_g_m) = c("Gene","Sample","Phosphorylation_score")

BRCA_Pho_score_g_m$Phosphorylation_score_plot = BRCA_Pho_score_g_m$Phosphorylation_score
BRCA_Pho_score_g_m$Phosphorylation_score_plot[BRCA_Pho_score_g_m$Phosphorylation_score_plot > 2]=2
BRCA_Pho_score_g_m$Phosphorylation_score_plot[BRCA_Pho_score_g_m$Phosphorylation_score_plot < -2]=-2

p = ggplot(data=BRCA_Pho_score_g_m)
p = p + geom_tile(aes(x=Sample, y=Gene, fill=as.numeric(Phosphorylation_score_plot), color = ifelse(as.numeric(Phosphorylation_score) > 1, "black",NA))) 
p = p + scale_fill_gradientn(name= "Phosphorylation score", na.value=NA, colours=RdBu1024, limit=c(-2,2))
p = p + scale_colour_manual(values=c("black",NA))
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(text = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              strip.text = element_text(size = 8),
              # panel.background = element_blank(),
              # panel.grid.major.x = element_blank(),
              # panel.grid.minor = element_blank(),
              legend.position="right")
plots[[1]] = p

# select specific pathways for plotting
sele_paths = c("MAPK signaling pathway","ErbB signaling pathway","Ras signaling pathway","p53 signaling pathway","PI3K-Akt signaling pathway"
               ,"mTOR signaling pathway","Cell cycle")
BRCA_Pho_pathway_all_merge_p = BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway %in% sele_paths,]

p = ggplot(data=BRCA_Pho_pathway_all_merge_p)
p = p + geom_point(aes(x=Sample, y=Pathway, fill=as.numeric(Global_phosphorylation), size=-log10(FDR), color=ifelse(Sig, "black",NA)),pch=21) 
p = p + scale_fill_gradientn(name= "Global Phosphorylation", na.value=NA, colours=RdBu1024, limit=c(-1.5,1.5))
p = p + scale_colour_manual(values=c("black",NA))
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(text = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              strip.text = element_text(size = 8),
              # panel.background = element_blank(),
              # #panel.grid.major.x = element_blank(),
              # panel.grid.minor = element_blank(),
              legend.position="right")
plots[[2]] = p

gp = do.call(rbind_gtable, plots)
# print the integrated plot
grid.newpage()
fn = paste(pd, 'BRCA_merged_pho_pathway.pdf',sep ="_")
pdf(fn, height=3, width=15,useDingbats = F)
grid.draw(gp)
dev.off()

### OV PNNL ###
OV_PNNL_Pho = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt")
OV_Pho = as.matrix(OV_PNNL_Pho)
OV_Pho_score = score_calc_by_gene(OV_Pho)

OV_fn = paste(source_date,"OV_cross_pathway_activation.tsv",sep="_")
OV_Pho_pathway_all_merge = read.table(header = T, quote="", sep = '\t', file=OV_fn)

##### plotting #####

plots = list()
OVGenes = c("TP53", "ERBB2", "MTOR","CDK1", "CDK2", "MAP3K1","RB1","PTEN")
OV_Pho_score_g = OV_Pho_score[row.names(OV_Pho_score) %in% OVGenes,]
OV_Pho_score_g_m = melt(OV_Pho_score_g)
colnames(OV_Pho_score_g_m) = c("Gene","Sample","Phosphorylation_score")

OV_Pho_score_g_m$Phosphorylation_score_plot = OV_Pho_score_g_m$Phosphorylation_score
OV_Pho_score_g_m$Phosphorylation_score_plot[OV_Pho_score_g_m$Phosphorylation_score_plot > 2]=2
OV_Pho_score_g_m$Phosphorylation_score_plot[OV_Pho_score_g_m$Phosphorylation_score_plot < -2]=-2

p = ggplot(data=OV_Pho_score_g_m)
p = p + geom_tile(aes(x=Sample, y=Gene, fill=as.numeric(Phosphorylation_score_plot), color = ifelse(as.numeric(Phosphorylation_score) > 1, "black",NA))) 
p = p + scale_fill_gradientn(name= "Phosphorylation score", na.value=NA, colours=RdBu1024, limit=c(-2,2))
p = p + scale_colour_manual(values=c("black",NA))
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(text = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              strip.text = element_text(size = 8),
              # panel.background = element_blank(),
              # panel.grid.major.x = element_blank(),
              # panel.grid.minor = element_blank(),
              legend.position="right")
plots[[1]] = p

# select specific pathways for plotting
sele_paths = c("MAPK signaling pathway","ErbB signaling pathway","Ras signaling pathway","p53 signaling pathway","PI3K-Akt signaling pathway"
               ,"mTOR signaling pathway","Cell cycle")
OV_Pho_pathway_all_merge_p = OV_Pho_pathway_all_merge[OV_Pho_pathway_all_merge$Pathway %in% sele_paths,]

p = ggplot(data=OV_Pho_pathway_all_merge_p)
p = p + geom_point(aes(x=Sample, y=Pathway, fill=as.numeric(Global_phosphorylation), size=-log10(FDR), color=ifelse(Sig, "black",NA)),pch=21) 
p = p + scale_fill_gradientn(name= "Global Phosphorylation", na.value=NA, colours=RdBu1024, limit=c(-1.5,1.5))
p = p + scale_colour_manual(values=c("black",NA))
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(text = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              strip.text = element_text(size = 8),
              #panel.background = element_blank(),
              #panel.grid.major.x = element_blank(),
              #panel.grid.minor = element_blank(),
              legend.position="right")
plots[[2]] = p


gp = do.call(rbind_gtable, plots)
# print the integrated plot
grid.newpage()
fn = paste(pd, 'OV_merged_pho_pathway.pdf',sep ="_")
pdf(fn, height=3, width=15,useDingbats = F)
grid.draw(gp)
dev.off()


##### merged kind of plot####
BRCA_merge = merge(BRCA_Pho_score_g_m,BRCA_Pho_pathway_all_merge_p,by="Sample")
OV_merge = merge(OV_Pho_score_g_m,OV_Pho_pathway_all_merge_p,by="Sample")

plot_gene_vs_pathway = function(merge_data, gene, pathway, name){
  #merge_data = BRCA_merge; gene = "CDK1"; pathway="Cell cycle"; name="BRCA"
  merge_data_sele = BRCA_merge[BRCA_merge$Gene==gene & BRCA_merge$Pathway==pathway,]
  merge_data_sele$outlier = merge_data_sele$Phosphorylation_score > 1
  
  f.test = fisher.test(table(merge_data_sele$outlier, merge_data_sele$Sig))
  
  getPalette = colorRampPalette(c("#fed976","#800026"))
  
  p = ggplot(data=merge_data_sele)
  p = p + geom_point(aes(x=Phosphorylation_score, y=Global_phosphorylation, color=ifelse(Sig, "black",NA), fill=-log10(FDR)), shape=21) 
  p = p + scale_fill_gradientn(name= "-log(FDR)", na.value=NA, colours=getPalette(100))
  p = p + scale_colour_manual(values=c("grey",NA))
  p = p + theme_bw()
  p = p + geom_vline(xintercept = 1, alpha=0.3)
  p = p + labs(x=paste(gene,"phosphorylation score"), y = paste(pathway, "global phosphorylation"))
  p = p + theme(legend.position="bottom")
  p
  fn = paste(pd, name,gene,pathway,'pho_pathway.pdf',sep ="_")
  ggsave(fn, w=5,h=6,useDingbats = F)
  return(f.test)
}

plot_gene_vs_pathway(BRCA_merge,"ERBB2","ErbB signaling pathway","BRCA")
plot_gene_vs_pathway(BRCA_merge,"PIK3CA","PI3K-Akt signaling pathway","BRCA")
plot_gene_vs_pathway(BRCA_merge,"MAPK3","MAPK signaling pathway","BRCA")
plot_gene_vs_pathway(BRCA_merge,"MAPK9","MAPK signaling pathway","BRCA")
plot_gene_vs_pathway(BRCA_merge,"TP53","p53 signaling pathway","BRCA")
plot_gene_vs_pathway(BRCA_merge,"CDK1","Cell cycle","BRCA")
plot_gene_vs_pathway(BRCA_merge,"MAP3K3","MAPK signaling pathway","BRCA")
