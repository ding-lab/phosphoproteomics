# Yige Wu @ WashU 2017 Feb
# look at correlation between kinase/substrate-based score and protein phos level in cancer genes


# choose kinase or phosphotase, outlier threshold and least sample number-------------------------
protein <- "kinase"
cancer <- "BRCA"
sig <- 0.1

# library -----------------------------------------------------------------
library(reshape)
library(stringr)
library(readr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))
source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# input correlation table and gene list -------------------------------------------------
# overlap <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/kinase_BRCA_reverse_marking_median_scored.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
overlap <- read_delim("~/Box Sync/pan3can_shared_data/analysis_results/tables/BRCA_ks_score_validated_ks_only_median_scored_kinase_table.txt", 
                      "\t", escape_double = FALSE, col_types = cols(n_sub_pho = col_double(),n_up_pho = col_double(), sub_phos_ctrans = col_double(), sub_phos_ctransscore = col_double(), up_phos_c = col_double(), up_phos_cscore = col_double()), trim_ws = TRUE)
ks_phos_corr <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/",protein,"_",cancer,"_kin:sub-based_phos_score_correlate_protein_phos_score.txt", sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
ks_phos_corr_sig <- ks_phos_corr[(!is.na(ks_phos_corr$s_fdr) & ks_phos_corr$s_fdr<=sig) | (!is.na(ks_phos_corr$k_fdr) & ks_phos_corr$k_fdr<=sig),]
load(paste(baseD,"pan3can_shared_data/analysis_results/2015-08-01_Gene_Set.RData",sep = ""))
kegg <- names(KEGG)

# BRCA_Pho_pathway_all_merge <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/pathway_activation/2017-02-15_KH_",cancer,"_cross_pathway_activation.tsv",sep = "") , 
#                             "\t", escape_double = FALSE, trim_ws = TRUE)
# BRCA_Pho_pathway_all_merge <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/pathway_activation/2017-02-18_KH_",cancer,"_cross_pathway_activation.tsv",sep = "") , 
#                                          "\t", escape_double = FALSE, trim_ws = TRUE)
BRCA_Pho_pathway_all_merge <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/pathway_activation/2017-02-23_KH_",cancer,"_cross_pathway_activation.tsv",sep = "") , 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)
cancer_genes <- read_delim(paste(baseD,"pan3can_shared_data/reference_files/Volgestin2013Science_125genes.txt", sep=""), 
                           "\t", escape_double = FALSE, col_types = cols(`Ocogene score*` = col_character()), 
                           trim_ws = TRUE)

BRCA_pho_f = paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv",sep="")
pho_data = read.delim(BRCA_pho_f)
pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))
pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")


out_thres <- 0.5
overlap$is.phos.outlier <- (overlap$kin_phos_cscore>=out_thres) & (!is.na(overlap$kin_phos_cscore))
overlap$is.up.pho.outlier <- (overlap$up_phos_cscore>=out_thres) & (!is.na(overlap$up_phos_cscore))
overlap$is.sub.pho.outlier <- (overlap$sub_phos_ctransscore>=out_thres) & (!is.na(overlap$sub_phos_ctransscore))
overlap$kinpho_subpho <- paste(overlap$is.up.pho.outlier, overlap$is.sub.pho.outlier, sep = "_")


for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  #for (path in "Cell cycle") {
  
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  pho_path$sample <- reorder(pho_path$Sample, pho_path$FDR)
  
  genes <- KEGG[[which(grepl(path,kegg))]]
  # genes_overlap <- intersect(genes, ks_phos_corr_sig$gene)
  genes_overlap <- genes
  rows <- c()
  for (gene in genes_overlap) {
    temp <- which(overlap$kinase == gene)
    rows <- c(rows,temp)
  }
  pho_table <- overlap[rows,c("kinase","sample","kin_phos_cscore")]
  pho_table$sig <- pho_table$kin_phos_cscore>1
  colnames(pho_table)[2] <- "Sample"
  pho_table <- merge(pho_table,pho_path[,c("Sample","FDR")],all.x = T)
  pho_table$sample <- reorder(pho_table$Sample,pho_table$FDR)
  
  p = ggplot(data=pho_path)
  p = p + geom_point(aes(x=sample, y=Pathway, fill=as.numeric(Global_phosphorylation), size=-log10(FDR), color=ifelse(Sig, "black",NA)),pch=21) 
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
  plots[[1]] = p
  
  lim = max(abs(max(pho_table$kin_phos_cscore, na.rm = T)),abs(min(pho_table$kin_phos_cscore, na.rm = T)))
  p = ggplot(data=pho_table)
  p = p + geom_tile(aes(x=sample, y=kinase, fill=kin_phos_cscore, colour=sig))
  p = p + scale_color_manual(values=c("TRUE" = "black","FALSE" = NA))
  p = p + scale_fill_gradientn(name= "phosphorylation score", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
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
  #cal_width = 20
  path_name_print <- gsub("\t", "_", path, fixed = TRUE)
  
  fn = paste(baseD,'pan3can_shared_data/analysis_results/pathway/',path_name_print,"_",cancer,"_",protein,'.pdf',sep ="")
  
  #fn = paste(pd, 'merged_mut_pathway.pdf',sep ="_")
  pdf(fn, height=6, width=15,useDingbats = F)
  grid.draw(gp)
  dev.off()
}

# violin ------------------------------------------------------------------
for (gene in cancer_genes$`Gene Symbol`) {
  # pho_split <- overlap[overlap$kinase==gene & overlap$is.phos.outlier,]
  pho_split <- overlap[overlap$kinase==gene,]
  p <- ggplot(pho_split, aes(kinpho_subpho, kin_phos_cscore))
  p = p + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
  p = p + geom_jitter(height = 0)
  p = p + theme_bw()
  p = p + theme_nogrid()
  p
  fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/phos_outlier_violin/',gene,'_',out_thres,'IQR.pdf',sep ="")
  ggsave(file=fn, height=3, width=4.5)
}

getPalette = colorRampPalette(c("#ffffcc","#fd8d3c","#800026"))
outlier.colors=c("NA", "#000000")

genes_ks_outlier <- unique(overlap$kinase[overlap$is.up.pho.outlier | overlap$is.sub.pho.outlier])

for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  #for (path in "Cell cycle") {
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  genes <- KEGG[[which(grepl(path,kegg))]]
  
  # genes_overlap <- intersect(genes, cancer_genes$`Gene Symbol`)
  genes_overlap <- intersect(genes, genes_ks_outlier)
  
  rows <- c()
  for (gene in genes_overlap) {
    temp <- which(overlap$kinase == gene)
    rows <- c(rows,temp)
  }
  pho_split <- overlap[rows,]
  colnames(pho_path)[2] <- "sample"
  pho_split <- merge(pho_split, pho_path[,c("sample","Global_phosphorylation", "Sig")], all.x = T)
  
  if (length(which(!is.na(pho_split$kin_phos_cscore))) > 0) {
    lim = max(abs(max(pho_split$kin_phos_cscore, na.rm = T)),abs(min(pho_split$kin_phos_cscore, na.rm = T)))
    
    p <- ggplot(pho_split, aes(x=kinpho_subpho,y= Global_phosphorylation, colour = kin_phos_cscore))
    p = p + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
    # p = p + geom_point(size = 2, alpha=0.8, stroke = 0) 
    p = p + geom_jitter(height = 0)
    p = p + scale_colour_gradientn(name= "protein_phos_score", na.value="grey40", colours=getPalette(100)) # , limit=c(-lim,lim)
    # p = p + geom_text(data=subset(pho_split, aa != "wt" ),
    #                   aes(x= kinpho_subpho, y = kin_phos_cscore ,label=aa), hjust = 0, nudge_x = -1.2 ,nudge_y = 0.2, size = 2)
    p = p + facet_wrap(~ kinase, nrow = 8 )
    p = p + theme_bw()
    #p = p + coord_flip()
    p
    fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/phos_outlier_violin/',path,'_',out_thres,'IQR.pdf',sep ="")
    ggsave(file=fn, height=20, width=20)
  }
}


# correlation genes for each pathway --------------------------------------

pathway <- c()
gene <- c()
pvalue <- c()
cor <- c()
count <- 0
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  #for (path in "Cell cycle") {
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  genes <- KEGG[[which(grepl(path,kegg))]]
  
  rows <- c()
  for (g in genes) {
    temp <- which(overlap$kinase == g)
    rows <- c(rows,temp)
  }
  
  pho_path_genes <- overlap[rows,]
  colnames(pho_path)[2] <- "sample"
  pho_path_genes <- merge(pho_path_genes, pho_path[,c("sample","Global_phosphorylation","FDR","Sig")], all.x = T)
  
  for (g in unique(pho_path_genes$kinase)) {
    corr_stat = try(cor.test(pho_path_genes$up_phos_cscore[pho_path_genes$kinase==g],pho_path_genes$Global_phosphorylation[pho_path_genes$kinase==g] , method = "pearson"), silent=T)
    if (!is(corr_stat,"try-error")){ # if the correlation test is carried out successfully
      count <- count + 1
      pathway[count] <- path
      gene[count] <- g
      cor[count] <- corr_stat$estimate
      pvalue[count] <- corr_stat$p.value
    }
  }
}
fdr <- p.adjust(pvalue, method = "fdr")
pathway_gene_pho_cor <- data.frame(pathway,gene,fdr,cor,pvalue)

# Global_phosphorylation ~ is.phos.outlier + is.up.pho.outlier + is.sub.pho.outlier,pvalues adjusted across pathway --------
pathway <- c()
protein <- c()
P_phos_outlier <- c(); P_kin_outlier <- c(); P_sub_outlier <- c()
coef_phos_outlier <- c(); coef_kin_outlier <- c(); coef_sub_outlier <- c()
count <- 0
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  #for (path in "MAPK signaling pathway") {
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  genes <- KEGG[[which(grepl(path,kegg))]]
  
  rows <- c()
  for (gene in genes) {
    temp <- which(overlap$kinase == gene)
    rows <- c(rows,temp)
  }
  pho_split <- overlap[rows,]
  colnames(pho_path)[2] <- "sample"
  pho_split <- merge(pho_split, pho_path[,c("sample","Global_phosphorylation", "Sig")], all.x = T)
  
  for (gene in unique(pho_split$kinase)) {
    pho_gene <- pho_split[pho_split$kinase==gene,]
    fit <- glm(Global_phosphorylation ~ is.phos.outlier + is.up.pho.outlier + is.sub.pho.outlier, data = pho_gene)
    
    count <- count + 1
    pathway[count] <- path
    protein[count] <- gene
    
    pvalues <- coef(summary(fit))
    
    if ( length(which(rownames(pvalues)=="is.phos.outlierTRUE")) > 0 ) {
      P_phos_outlier[count] <- pvalues["is.phos.outlierTRUE",4]
    } else {
      P_phos_outlier[count] <- NA
    }
    if ( length(which(rownames(pvalues)=="is.up.pho.outlierTRUE")) > 0 ) {
      P_kin_outlier[count] <- pvalues["is.up.pho.outlierTRUE",4]
    } else {
      P_kin_outlier[count] <- NA
    }
    if ( length(which(rownames(pvalues)=="is.sub.pho.outlierTRUE")) > 0 ) {
      P_sub_outlier[count] <- pvalues["is.sub.pho.outlierTRUE",4]
    } else {
      P_sub_outlier[count] <- NA
    }
    
    coef_phos_outlier[count] <- fit$coefficients[2];
    coef_kin_outlier[count] <- fit$coefficients[3]; 
    coef_sub_outlier[count] <- fit$coefficients[4]
  }
  # fdr_temp <- p.adjust(pvalue[pathway==path], method = "fdr")
  # fdr <- c(fdr,fdr_temp)
}
fdr_phos_outlier <- p.adjust(P_phos_outlier, method = "fdr")
fdr_kin_outlier <- p.adjust(P_kin_outlier, method = "fdr")
fdr_sub_outlier <- p.adjust(P_sub_outlier, method = "fdr")

pathway_pho_reg <- data.frame(pathway,protein,
                              fdr_phos_outlier,fdr_kin_outlier,fdr_sub_outlier,
                              coef_phos_outlier,coef_kin_outlier,coef_sub_outlier,
                              P_phos_outlier,P_kin_outlier,P_sub_outlier)

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_phos_regression_adjust_pvalue_across_pathways_BRCA.txt", sep="")
write.table(pathway_pho_reg, file=tn, quote=F, sep = '\t', row.names = FALSE)

# Global_phosphorylation ~ is.phos.outlier + is.up.pho.outlier + is.sub.pho.outlier,pvalues adjusted within pathways --------
pathway <- c()
protein <- c()
P_phos_outlier <- c(); P_kin_outlier <- c(); P_sub_outlier <- c()
coef_phos_outlier <- c(); coef_kin_outlier <- c(); coef_sub_outlier <- c()
fdr_phos_outlier <- c(); fdr_kin_outlier <- c(); fdr_sub_outlier <- c()
count <- 0
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  #for (path in "MAPK signaling pathway") {
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  genes <- KEGG[[which(grepl(path,kegg))]]
  
  rows <- c()
  for (gene in genes) {
    temp <- which(overlap$kinase == gene)
    rows <- c(rows,temp)
  }
  pho_split <- overlap[rows,]
  colnames(pho_path)[2] <- "sample"
  pho_split <- merge(pho_split, pho_path[,c("sample","Global_phosphorylation", "Sig")], all.x = T)
  
  for (gene in unique(pho_split$kinase)) {
    pho_gene <- pho_split[pho_split$kinase==gene,]
    fit <- glm(Global_phosphorylation ~ is.phos.outlier + is.up.pho.outlier + is.sub.pho.outlier, data = pho_gene)
    
    count <- count + 1
    pathway[count] <- path
    protein[count] <- gene
    
    pvalues <- coef(summary(fit))
    
    if ( length(which(rownames(pvalues)=="is.phos.outlierTRUE")) > 0 ) {
      P_phos_outlier[count] <- pvalues["is.phos.outlierTRUE",4]
    } else {
      P_phos_outlier[count] <- NA
    }
    if ( length(which(rownames(pvalues)=="is.up.pho.outlierTRUE")) > 0 ) {
      P_kin_outlier[count] <- pvalues["is.up.pho.outlierTRUE",4]
    } else {
      P_kin_outlier[count] <- NA
    }
    if ( length(which(rownames(pvalues)=="is.sub.pho.outlierTRUE")) > 0 ) {
      P_sub_outlier[count] <- pvalues["is.sub.pho.outlierTRUE",4]
    } else {
      P_sub_outlier[count] <- NA
    }
    
    coef_phos_outlier[count] <- fit$coefficients[2];
    coef_kin_outlier[count] <- fit$coefficients[3]; 
    coef_sub_outlier[count] <- fit$coefficients[4]
  }
  fdr_phos_outlier <- c(fdr_phos_outlier,p.adjust(P_phos_outlier[pathway==path], method = "fdr"))
  fdr_kin_outlier <- c(fdr_kin_outlier,p.adjust(P_kin_outlier[pathway==path], method = "fdr"))
  fdr_sub_outlier <- c(fdr_sub_outlier,p.adjust(P_sub_outlier[pathway==path], method = "fdr"))
}

pathway_pho_reg <- data.frame(pathway,protein,
                              fdr_phos_outlier,fdr_kin_outlier,fdr_sub_outlier,
                              coef_phos_outlier,coef_kin_outlier,coef_sub_outlier,
                              P_phos_outlier,P_kin_outlier,P_sub_outlier)

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_phos_regression_adjust_pvalue_within_pathways_BRCA.txt", sep="")
write.table(pathway_pho_reg, file=tn, quote=F, sep = '\t', row.names = FALSE)

rows <- c()
for (gene in cancer_genes$`Gene Symbol`) {
  temp <- which(pathway_pho_reg$protein==gene)
  rows <- c(rows,temp)
}
pathway_pho_reg_cancer <- pathway_pho_reg[rows,]
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_phos_regression_adjust_pvalue_cancer_genes_within_pathways_BRCA.txt", sep="")
write.table(pathway_pho_reg_cancer, file=tn, quote=F, sep = '\t', row.names = FALSE)



# Global_phosphorylation ~ phos_level + kin_pho_level + sub_pho_level,pvalues adjusted within pathways --------
kegg <- names(KEGG)

pathway <- c()
protein <- c()
P_phos_level <- c(); P_kin_level <- c(); P_sub_level <- c()
coef_phos_level <- c(); coef_kin_level <- c(); coef_sub_level <- c()
fdr_phos_level <- c(); fdr_kin_level <- c(); fdr_sub_level <- c()
count <- 0
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
#for (path in "ErbB signaling pathway") {
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  genes <- KEGG[[which(grepl(path,kegg))]]
  
  rows <- c()
  for (gene in genes) {
    temp <- which(overlap$kinase == gene)
    rows <- c(rows,temp)
  }
  pho_split <- overlap[rows,]
  colnames(pho_path)[2] <- "sample"
  pho_split <- merge(pho_split, pho_path[,c("sample","Global_phosphorylation", "Sig")], all.x = T)
  
  for (gene in unique(pho_split$kinase)) {
    pho_gene <- data.frame(pho_split[pho_split$kinase==gene,])
    if (!all(is.na(pho_gene$up_phos_c)) & !all(is.na(pho_gene$sub_phos_ctrans))) {
      fit <- glm(Global_phosphorylation ~ kin_phos_c + up_phos_c + sub_phos_ctrans, data = pho_gene)
    } else if (!all(is.na(pho_gene$sub_phos_ctrans))) {
      fit <- glm(Global_phosphorylation ~ kin_phos_c + sub_phos_ctrans, data = pho_gene)
    } else if (!all(is.na(pho_gene$up_phos_c))){
      fit <- glm(Global_phosphorylation ~ kin_phos_c + up_phos_c, data = pho_gene)
    } else if (!all(is.na(pho_gene$kin_phos_c))) {
      fit <- glm(Global_phosphorylation ~ kin_phos_c, data = pho_gene)
    }
    
    count <- count + 1
    pathway[count] <- path
    protein[count] <- gene
    
    pvalues <- coef(summary(fit))
    
    if ( length(which(rownames(pvalues)=="kin_phos_c")) > 0 ) {
      P_phos_level[count] <- pvalues["kin_phos_c",4]
    } else {
      P_phos_level[count] <- NA
    }
    if ( length(which(rownames(pvalues)=="up_phos_c")) > 0 ) {
      P_kin_level[count] <- pvalues["up_phos_c",4]
    } else {
      P_kin_level[count] <- NA
    }
    if ( length(which(rownames(pvalues)=="sub_phos_ctrans")) > 0 ) {
      P_sub_level[count] <- pvalues["sub_phos_ctrans",4]
    } else {
      P_sub_level[count] <- NA
    }
    coef_phos_level[count] <- fit$coefficients["kin_phos_c"];
    coef_kin_level[count] <- fit$coefficients["up_phos_c"]; 
    coef_sub_level[count] <- fit$coefficients["sub_phos_ctrans"]
  }
  fdr_phos_level <- c(fdr_phos_level,p.adjust(P_phos_level[pathway==path], method = "fdr"))
  fdr_kin_level <- c(fdr_kin_level,p.adjust(P_kin_level[pathway==path], method = "fdr"))
  fdr_sub_level <- c(fdr_sub_level,p.adjust(P_sub_level[pathway==path], method = "fdr"))
}

pathway_pho_level_reg <- data.frame(pathway,protein,
                                    fdr_phos_level,fdr_kin_level,fdr_sub_level,
                                    coef_phos_level,coef_kin_level,coef_sub_level,
                                    P_phos_level,P_kin_level,P_sub_level)

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_phos_level_regression_adjust_pvalue_within_pathways_BRCA.txt", sep="")
write.table(pathway_pho_level_reg, file=tn, quote=F, sep = '\t', row.names = FALSE)

rows <- c()
for (gene in cancer_genes$`Gene Symbol`) {
  temp <- which(pathway_pho_reg$protein==gene)
  rows <- c(rows,temp)
}
pathway_pho_reg_cancer <- pathway_pho_reg[rows,]
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_phos_regression_adjust_pvalue_cancer_genes_within_pathways_BRCA.txt", sep="")
write.table(pathway_pho_reg_cancer, file=tn, quote=F, sep = '\t', row.names = FALSE)

# pathway_pho&phophosite_pho_corr -----------------------------------------

pathway <- c()
protein <- c()
rsd <- c()
pvalue <- c()
cor <- c()
fdr <- c()
count <- 0
#for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
for (path in "ErbB signaling pathway") {
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  genes <- KEGG[[which(grepl(path,kegg))]]
  genes_overlap <- intersect(genes,unique(pho_rsd_split$SUBSTRATE))
  for (g in genes_overlap) {
    rsd_rows <- which(pho_rsd_split$SUBSTRATE==g)
    for (i in rsd_rows) {
      corr_stat = try(cor.test(pho_path$Global_phosphorylation,
                               as.numeric(pho_data[i,-1]),
                               method = "pearson"), silent=T)
      if (!is(corr_stat,"try-error")){ # if the correlation test is carried out successfully
        count <- count + 1
        pathway[count] <- path
        protein[count] <- g
        rsd[count] <- pho_rsd_split$SUB_MOD_RSD[i]
        cor[count] <- corr_stat$estimate
        pvalue[count] <- corr_stat$p.value
      }
    }
  }
  fdr <- c(fdr,p.adjust(pvalue[pathway==path], method = "fdr"))
}
#fdr <- p.adjust(pvalue, method = "fdr")
pathway_site_pho_cor <- data.frame(pathway,protein,rsd,fdr,cor,pvalue)

pathway_site_pho_cor_ordr <- c()
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  temp <- pathway_site_pho_cor[pathway_site_pho_cor$pathway==path,]
  temp <- temp[order(temp$pvalue),]
  pathway_site_pho_cor_ordr <- rbind(pathway_site_pho_cor_ordr,temp)
}
pathway_site_pho_cor_ordr$cor_type[pathway_site_pho_cor_ordr$cor>0] <- "positively-correlated"
pathway_site_pho_cor_ordr$cor_type[pathway_site_pho_cor_ordr$cor<0] <- "negatively-correlated"

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_phos_&_phosphosite_corr_adjust_pvalue_within_pathways_BRCA.txt", sep="")
write.table(pathway_site_pho_cor_ordr, file=tn, quote=F, sep = '\t', row.names = FALSE)

# tailored_outlier_overlap ------------------------------------------------
outlier_overlap <- overlap[,c("kinase","sample","up_phos_cscore","sub_phos_ctransscore")]
colnames(outlier_overlap)[1] <- "protein"
center <- 0.5

for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  #for (path in "Cell cycle") {
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  pho_path$sample <- reorder(pho_path$Sample, pho_path$Global_phosphorylation)
  
  p = ggplot(data=pho_path)
  p = p + geom_point(aes(x=sample, y=Global_phosphorylation, fill=Sig),pch=21)
  p = p + facet_grid(Pathway~.,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
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
  
  cor_table <- pathway_site_pho_cor_ordr[pathway_site_pho_cor_ordr$pathway==path,]
  proteins <- as.vector(cor_table$protein)
  rsds <- as.vector(cor_table$rsd)
  
  rows <- c()
  len <- min(nrow(cor_table),top)
  for (i in 1:len) {
    temp <- which(pho_rsd_split$SUBSTRATE == proteins[i] & pho_rsd_split$SUB_MOD_RSD == rsds[i])
    rows <- c(rows,temp)
  }
  
  phosite_score <- pho_data[rows,-1]
  for (i in 1:nrow(phosite_score)) {
    temp <- phosite_score[i,]
    IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
    phosite_score[i,] = ( temp - quantile(temp, probs=center, na.rm=T))/IQR
  }
  
  pho_table <- cbind(pho_rsd_split[rows,c("SUBSTRATE","SUB_MOD_RSD")],phosite_score)
  
  pho_table_m <- melt(pho_table,id = c("SUBSTRATE","SUB_MOD_RSD"))
  colnames(pho_table_m)[1:3] <- c("protein","rsd","Sample")
  pho_table_m <- merge(pho_table_m,pho_path[,c("Global_phosphorylation","Sample")],all.x = T)
  pho_table_m <- merge(pho_table_m,cor_table[,c("protein","rsd","cor_type","pvalue")], all.x = T)
  pho_table_m$sample <- reorder(pho_table_m$Sample,pho_table_m$Global_phosphorylation)
  pho_table_m$phosphosite <- paste(pho_table_m$protein,pho_table_m$rsd,sep = ":")
  pho_table_m$Phosphosite <- reorder(pho_table_m$phosphosite,-pho_table_m$pvalue)
  pho_table_m <- merge(pho_table_m,outlier_overlap, all.x = T)
  pho_table_m$is.kin.outlier <- (!is.na(pho_table_m$up_phos_cscore)) & ((pho_table_m$cor_type=="positively-correlated" & pho_table_m$up_phos_cscore > out_thres) | (pho_table_m$cor_type=="negatively-correlated" & pho_table_m$up_phos_cscore < - out_thres))
  pho_table_m$is.sub.outlier <- (!is.na(pho_table_m$sub_phos_ctransscore)) & ((pho_table_m$cor_type=="positively-correlated" & pho_table_m$sub_phos_ctransscore > out_thres) | (pho_table_m$cor_type=="negatively-correlated" & pho_table_m$sub_phos_ctransscore < - out_thres))
  pho_table_m$kinorsub <- pho_table_m$is.kin.outlier | pho_table_m$is.sub.outlier
  pho_table_m$is.self.outlier <- (!is.na(pho_table_m$value)) & ((pho_table_m$cor_type=="positively-correlated" & pho_table_m$value > 1) | (pho_table_m$cor_type=="negatively-correlated" & pho_table_m$value < - 1))
  pho_table_m$self_kinorsub <- paste(pho_table_m$is.self.outlier, pho_table_m$kinorsub,sep = "_")
  
  lim = max(abs(max(pho_table_m$value, na.rm = T)),abs(min(pho_table_m$value, na.rm = T)))
  p = ggplot(data=pho_table_m)
  p = p + geom_tile(aes(x=sample, y=Phosphosite, fill=value, color=self_kinorsub, width=0.7, height=0.7), size = 0.5)
  p = p + scale_color_manual(values=c("FALSE_FALSE" = NA,"FALSE_TRUE" = NA ,"TRUE_FALSE" = "yellow", "TRUE_TRUE" = "red"))
  p = p + scale_fill_gradientn(name= "phosphorylation level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
  p = p + facet_grid(protein~.,scales = "free_y", space = "free_y", drop=T)#, , space = "free_y",scales = "free_y")#, space = "free", scales = "free")
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
  p
  plots[[2]] = p
  
  gp = do.call(rbind_gtable, plots)
  # print the integrated plot
  grid.newpage()
  #cal_width = 20
  
  fn = paste(baseD,'pan3can_shared_data/analysis_results/pathway/tailored_outlier/',path,"_correlated_phosphosites_outlier_marked_top_",len,"_",cancer,"_kinase_table.pdf",sep ="")
  pdf(fn, height=6, width=15,useDingbats = F)
  grid.draw(gp)
  dev.off()
}


# positive_outlier_overlap ------------------------------------------------
outlier_overlap <- overlap[,c("kinase","sample","up_phos_cscore","sub_phos_ctransscore")]
colnames(outlier_overlap)[1] <- "protein"
center <- 0.5

for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  #for (path in "Cell cycle") {
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  pho_path$sample <- reorder(pho_path$Sample, pho_path$Global_phosphorylation)
  
  p = ggplot(data=pho_path)
  p = p + geom_point(aes(x=sample, y=Global_phosphorylation, fill=Sig),pch=21)
  p = p + facet_grid(Pathway~.,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
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
  
  cor_table <- pathway_site_pho_cor_ordr[pathway_site_pho_cor_ordr$pathway==path,]
  proteins <- as.vector(cor_table$protein)
  rsds <- as.vector(cor_table$rsd)
  
  rows <- c()
  len <- min(nrow(cor_table),top)
  for (i in 1:len) {
    temp <- which(pho_rsd_split$SUBSTRATE == proteins[i] & pho_rsd_split$SUB_MOD_RSD == rsds[i])
    rows <- c(rows,temp)
  }
  
  phosite_score <- pho_data[rows,-1]
  for (i in 1:nrow(phosite_score)) {
    temp <- phosite_score[i,]
    IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
    phosite_score[i,] = ( temp - quantile(temp, probs=center, na.rm=T))/IQR
  }
  
  pho_table <- cbind(pho_rsd_split[rows,c("SUBSTRATE","SUB_MOD_RSD")],phosite_score)
  
  pho_table_m <- melt(pho_table,id = c("SUBSTRATE","SUB_MOD_RSD"))
  colnames(pho_table_m)[1:3] <- c("protein","rsd","Sample")
  pho_table_m <- merge(pho_table_m,pho_path[,c("Global_phosphorylation","Sample")],all.x = T)
  pho_table_m <- merge(pho_table_m,cor_table[,c("protein","rsd","cor_type","pvalue")], all.x = T)
  pho_table_m$sample <- reorder(pho_table_m$Sample,pho_table_m$Global_phosphorylation)
  pho_table_m$phosphosite <- paste(pho_table_m$protein,pho_table_m$rsd,sep = ":")
  pho_table_m$Phosphosite <- reorder(pho_table_m$phosphosite,-pho_table_m$pvalue)
  pho_table_m <- merge(pho_table_m,outlier_overlap, all.x = T)
  pho_table_m$is.kin.outlier <- (!is.na(pho_table_m$up_phos_cscore)) & (pho_table_m$up_phos_cscore > out_thres)
  pho_table_m$is.sub.outlier <- (!is.na(pho_table_m$sub_phos_ctransscore)) & (pho_table_m$sub_phos_ctransscore > out_thres)
  pho_table_m$kinorsub <- pho_table_m$is.kin.outlier | pho_table_m$is.sub.outlier
  pho_table_m$is.self.outlier <- (!is.na(pho_table_m$value)) & (pho_table_m$value > 1)
  pho_table_m$self_kinorsub <- paste(pho_table_m$is.self.outlier, pho_table_m$kinorsub,sep = "_")
  
  lim = max(abs(max(pho_table_m$value, na.rm = T)),abs(min(pho_table_m$value, na.rm = T)))
  p = ggplot(data=pho_table_m)
  p = p + geom_tile(aes(x=sample, y=Phosphosite, fill=value, color=self_kinorsub, width=0.7, height=0.7), size = 0.5)
  p = p + scale_color_manual(values=c("FALSE_FALSE" = NA,"FALSE_TRUE" = NA ,"TRUE_FALSE" = "yellow", "TRUE_TRUE" = "red"))
  p = p + scale_fill_gradientn(name= "phosphorylation level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
  p = p + facet_grid(protein~.,scales = "free_y", space = "free_y", drop=T)#, , space = "free_y",scales = "free_y")#, space = "free", scales = "free")
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
  p
  plots[[2]] = p
  
  gp = do.call(rbind_gtable, plots)
  # print the integrated plot
  grid.newpage()
  #cal_width = 20
  
  fn = paste(baseD,'pan3can_shared_data/analysis_results/pathway/positive_outlier/',path,"_correlated_phosphosites_outlier_marked_top_",len,"_",cancer,"_kinase_table.pdf",sep ="")
  pdf(fn, height=6, width=15,useDingbats = F)
  grid.draw(gp)
  dev.off()
}

# tailored_outlier_overlap,mark_kinorsub_only_also ------------------------------------------------
outlier_overlap <- overlap[,c("kinase","sample","up_phos_cscore","sub_phos_ctransscore")]
colnames(outlier_overlap)[1] <- "protein"
center <- 0.5

for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  #for (path in "Cell cycle") {
  pho_path <- BRCA_Pho_pathway_all_merge[BRCA_Pho_pathway_all_merge$Pathway==path,]
  pho_path$sample <- reorder(pho_path$Sample, pho_path$Global_phosphorylation)
  
  p = ggplot(data=pho_path)
  p = p + geom_point(aes(x=sample, y=Global_phosphorylation, fill=Sig),pch=21)
  p = p + facet_grid(Pathway~.,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
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
  
  cor_table <- pathway_site_pho_cor_ordr[pathway_site_pho_cor_ordr$pathway==path,]
  proteins <- as.vector(cor_table$protein)
  rsds <- as.vector(cor_table$rsd)
  
  rows <- c()
  len <- min(nrow(cor_table),top)
  for (i in 1:len) {
    temp <- which(pho_rsd_split$SUBSTRATE == proteins[i] & pho_rsd_split$SUB_MOD_RSD == rsds[i])
    rows <- c(rows,temp)
  }
  
  phosite_score <- pho_data[rows,-1]
  for (i in 1:nrow(phosite_score)) {
    temp <- phosite_score[i,]
    IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
    phosite_score[i,] = ( temp - quantile(temp, probs=center, na.rm=T))/IQR
  }
  
  pho_table <- cbind(pho_rsd_split[rows,c("SUBSTRATE","SUB_MOD_RSD")],phosite_score)
  
  pho_table_m <- melt(pho_table,id = c("SUBSTRATE","SUB_MOD_RSD"))
  colnames(pho_table_m)[1:3] <- c("protein","rsd","Sample")
  pho_table_m <- merge(pho_table_m,pho_path[,c("Global_phosphorylation","Sample")],all.x = T)
  pho_table_m <- merge(pho_table_m,cor_table[,c("protein","rsd","cor_type","pvalue")], all.x = T)
  pho_table_m$sample <- reorder(pho_table_m$Sample,pho_table_m$Global_phosphorylation)
  pho_table_m$phosphosite <- paste(pho_table_m$protein,pho_table_m$rsd,sep = ":")
  pho_table_m$Phosphosite <- reorder(pho_table_m$phosphosite,-pho_table_m$pvalue)
  pho_table_m <- merge(pho_table_m,outlier_overlap, all.x = T)
  pho_table_m$is.kin.outlier <- (!is.na(pho_table_m$up_phos_cscore)) & ((pho_table_m$cor_type=="positively-correlated" & pho_table_m$up_phos_cscore > 1) | (pho_table_m$cor_type=="negatively-correlated" & pho_table_m$up_phos_cscore < -1))
  pho_table_m$is.sub.outlier <- (!is.na(pho_table_m$sub_phos_ctransscore)) & ((pho_table_m$cor_type=="positively-correlated" & pho_table_m$sub_phos_ctransscore > 1) | (pho_table_m$cor_type=="negatively-correlated" & pho_table_m$sub_phos_ctransscore < -1))
  pho_table_m$kinorsub <- pho_table_m$is.kin.outlier | pho_table_m$is.sub.outlier
  pho_table_m$is.self.outlier <- (!is.na(pho_table_m$value)) & ((pho_table_m$cor_type=="positively-correlated" & pho_table_m$value > 1) | (pho_table_m$cor_type=="negatively-correlated" & pho_table_m$value < - 1))
  pho_table_m$self_kinorsub <- paste(pho_table_m$is.self.outlier, pho_table_m$kinorsub,sep = "_")
  
  lim = max(abs(max(pho_table_m$value, na.rm = T)),abs(min(pho_table_m$value, na.rm = T)))
  p = ggplot(data=pho_table_m)
  p = p + geom_tile(aes(x=sample, y=Phosphosite, fill=value, color=self_kinorsub, width=0.7, height=0.7), size = 0.5)
  p = p + scale_color_manual(values=c("FALSE_FALSE" = NA,"FALSE_TRUE" = "yellow" ,"TRUE_FALSE" = "orange", "TRUE_TRUE" = "red"))
  p = p + scale_fill_gradientn(name= "phosphorylation level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
  p = p + facet_grid(protein~.,scales = "free_y", space = "free_y", drop=T)#, , space = "free_y",scales = "free_y")#, space = "free", scales = "free")
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
  p
  plots[[2]] = p
  
  gp = do.call(rbind_gtable, plots)
  # print the integrated plot
  grid.newpage()
  #cal_width = 20
  
  fn = paste(baseD,'pan3can_shared_data/analysis_results/pathway/tailored_outlier/mark_kinorsub_only_also/',path,"_correlated_phosphosites_outlier_marked_top_",len,"_",cancer,"_kinase_table.pdf",sep ="")
  pdf(fn, height=6, width=15,useDingbats = F)
  grid.draw(gp)
  dev.off()
}

# mark pathway phosphoryaltion with subtypes ------------------------------
for (subtype in c("Her2","Basal","LumA","LumB")) {
  sam_subtype <- colnames(clinical)[clinical[1,]==subtype & !is.na(clinical[1,])]
  for (sam in sam_subtype) {
    pho_pathway_table$subtype[pho_pathway_table$Sample==sam] <- subtype
  }
}

# Global_phosphorylation ~ phos_level + kin_pho_level + sub_pho_level,divided by subtypes,pvalues adjusted within pathways --------
kegg <- names(KEGG)
least_samples <- 5
pathway_pho_level_reg_subtypes <- c()
for (subtype in c("Her2","Basal","LumA","LumB")) {
  sam_subtype <- colnames(clinical)[clinical[1,]==subtype & !is.na(clinical[1,])]
  
  pathway <- c()
  protein <- c()
  P_phos_level <- c(); P_kin_level <- c(); P_sub_level <- c()
  coef_phos_level <- c(); coef_kin_level <- c(); coef_sub_level <- c()
  fdr_phos_level <- c(); fdr_kin_level <- c(); fdr_sub_level <- c()
  count <- 0
  for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
    #for (path in "ErbB signaling pathway") {
    pho_path_all <- pho_pathway_table[pho_pathway_table$Pathway==path,]
    rownames(pho_path_all) <- pho_path_all$Sample
    pho_path <- pho_path_all[sam_subtype,]
    colnames(pho_path)[2] <- "sample"
    
    genes <- KEGG[[which(grepl(path,kegg))]]
    rows <- c()
    for (gene in genes) {
      temp <- which(overlap$kinase == gene)
      rows <- c(rows,temp)
    }
    pho_split_all <- overlap[rows,]
    
    rows <- c()
    for (sam in sam_subtype) {
      temp <- which(pho_split_all$sample==sam)
      rows <- c(rows,temp)
    }
    pho_split <- pho_split_all[rows,]
    
    pho_split <- merge(pho_split, pho_path[,c("sample","Global_phosphorylation", "Sig")], all.x = T)
    
    for (gene in unique(pho_split$kinase)) {
      pho_gene <- data.frame(pho_split[pho_split$kinase==gene,])
      if (length(which(complete.cases(pho_gene[,c("kin_phos_c","up_phos_c","sub_phos_ctrans","Global_phosphorylation")]))) >= least_samples) {
        fit <- glm(Global_phosphorylation ~ kin_phos_c + up_phos_c + sub_phos_ctrans, data = pho_gene); success <- 1
      } else if (length(which(complete.cases(pho_gene[,c("kin_phos_c","sub_phos_ctrans","Global_phosphorylation")]))) >= least_samples) {
        fit <- glm(Global_phosphorylation ~ kin_phos_c + sub_phos_ctrans, data = pho_gene); success <- 1
      } else if (length(which(complete.cases(pho_gene[,c("kin_phos_c","up_phos_c","Global_phosphorylation")]))) >= least_samples) {
        fit <- glm(Global_phosphorylation ~ kin_phos_c + up_phos_c, data = pho_gene); success <- 1
      } else if (length(which(complete.cases(pho_gene[,c("kin_phos_c","Global_phosphorylation")]))) >= least_samples) {
        fit <- glm(Global_phosphorylation ~ kin_phos_c, data = pho_gene); success <- 1
      } else {
        success <- 0
      }
      if (success == 1) {
        count <- count + 1
        pathway[count] <- path
        protein[count] <- gene
        
        pvalues <- coef(summary(fit))
        
        if ( length(which(rownames(pvalues)=="kin_phos_c")) > 0 ) {
          P_phos_level[count] <- pvalues["kin_phos_c",4]
        } else {
          P_phos_level[count] <- NA
        }
        if ( length(which(rownames(pvalues)=="up_phos_c")) > 0 ) {
          P_kin_level[count] <- pvalues["up_phos_c",4]
        } else {
          P_kin_level[count] <- NA
        }
        if ( length(which(rownames(pvalues)=="sub_phos_ctrans")) > 0 ) {
          P_sub_level[count] <- pvalues["sub_phos_ctrans",4]
        } else {
          P_sub_level[count] <- NA
        }
        coef_phos_level[count] <- fit$coefficients["kin_phos_c"];
        coef_kin_level[count] <- fit$coefficients["up_phos_c"]; 
        coef_sub_level[count] <- fit$coefficients["sub_phos_ctrans"]
      }
    }
    fdr_phos_level <- c(fdr_phos_level,p.adjust(P_phos_level[pathway==path], method = "fdr"))
    fdr_kin_level <- c(fdr_kin_level,p.adjust(P_kin_level[pathway==path], method = "fdr"))
    fdr_sub_level <- c(fdr_sub_level,p.adjust(P_sub_level[pathway==path], method = "fdr"))
  }
  pathway_pho_level_reg <- data.frame(pathway,protein,
                                      fdr_phos_level,fdr_kin_level,fdr_sub_level,
                                      coef_phos_level,coef_kin_level,coef_sub_level,
                                      P_phos_level,P_kin_level,P_sub_level)
  
  pathway_pho_level_reg$subtype <- subtype
  pathway_pho_level_reg_subtypes <- rbind(pathway_pho_level_reg_subtypes, pathway_pho_level_reg)
}

pathway_pho_level_reg_sub_ordr <- c()
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  temp <- pathway_pho_level_reg_subtypes[pathway_pho_level_reg_subtypes$pathway==path,]
  temp <- temp[order(temp$P_phos_level),]
  pathway_pho_level_reg_sub_ordr <- rbind(pathway_pho_level_reg_sub_ordr,temp)
}

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_phos_collapsed_subtype_regression_adjust_pvalue_within_pathways_BRCA.txt", sep="")
write.table(pathway_pho_level_reg_sub_ordr, file=tn, quote=F, sep = '\t', row.names = FALSE)





rows <- c()
for (gene in cancer_genes$`Gene Symbol`) {
  temp <- which(pathway_pho_reg$protein==gene)
  rows <- c(rows,temp)
}
pathway_pho_reg_cancer <- pathway_pho_reg[rows,]
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_phos_regression_adjust_pvalue_cancer_genes_within_pathways_BRCA.txt", sep="")
write.table(pathway_pho_reg_cancer, file=tn, quote=F, sep = '\t', row.names = FALSE)


# pathway_pho&phophosite_pho_corr, divided by subtypes -----------------------------------------
clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))
pho_pathway_table <- data.frame(BRCA_Pho_pathway_all_merge)

pathway_site_pho_cor_subtypes <- c()
for (subtype in c("Her2","Basal","LumA","LumB")) {
  sam_subtype <- colnames(clinical)[clinical[1,]==subtype & !is.na(clinical[1,])]
  
  pathway <- c()
  protein <- c()
  rsd <- c()
  pvalue <- c()
  cor <- c()
  fdr <- c()
  count <- 0
  for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
    #for (path in "ErbB signaling pathway") {
    pho_path_all <- pho_pathway_table[pho_pathway_table$Pathway==path,]
    rownames(pho_path_all) <- pho_path_all$Sample
    pho_path <- pho_path_all[sam_subtype,]
    genes <- KEGG[[which(grepl(path,kegg))]]
    genes_overlap <- intersect(genes,unique(pho_rsd_split$SUBSTRATE))
    for (g in genes_overlap) {
      rsd_rows <- which(pho_rsd_split$SUBSTRATE==g)
      for (i in rsd_rows) {
        corr_stat = try(cor.test(pho_path$Global_phosphorylation,
                                 as.numeric(pho_data[i,sam_subtype]),
                                 method = "pearson"), silent=T)
        if (!is(corr_stat,"try-error")){ # if the correlation test is carried out successfully
          count <- count + 1
          pathway[count] <- path
          protein[count] <- g
          rsd[count] <- pho_rsd_split$SUB_MOD_RSD[i]
          cor[count] <- corr_stat$estimate
          pvalue[count] <- corr_stat$p.value
        }
      }
    }
    fdr <- c(fdr,p.adjust(pvalue[pathway==path], method = "fdr"))
  }
  pathway_site_pho_cor <- data.frame(pathway,protein,rsd,fdr,cor,pvalue)
  pathway_site_pho_cor$subtype <- subtype
  pathway_site_pho_cor_subtypes <- rbind(pathway_site_pho_cor_subtypes, pathway_site_pho_cor)
}

pathway_site_pho_cor_sub_ordr <- c()
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  temp <- pathway_site_pho_cor_subtypes[pathway_site_pho_cor_subtypes$pathway==path,]
  temp <- temp[order(temp$pvalue),]
  pathway_site_pho_cor_sub_ordr <- rbind(pathway_site_pho_cor_sub_ordr,temp)
}
pathway_site_pho_cor_sub_ordr$cor_type[pathway_site_pho_cor_sub_ordr$cor>0] <- "positively-correlated"
pathway_site_pho_cor_sub_ordr$cor_type[pathway_site_pho_cor_sub_ordr$cor<0] <- "negatively-correlated"

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_phos_&_phosphosite_corr_subtypes_adjust_pvalue_within_pathways_BRCA.txt", sep="")
write.table(pathway_site_pho_cor_sub_ordr, file=tn, quote=F, sep = '\t', row.names = FALSE)

# Global_phosphorylation ~ phos_level(site) + kin_pho_level + sub_pho_level,divided by subtypes,pvalues adjusted within pathways --------
pathway_site_kinsub_cor_subtypes <- c()
for (subtype in c("Her2","Basal","LumA","LumB")) {
  sam_subtype <- colnames(clinical)[clinical[1,]==subtype & !is.na(clinical[1,])]
  
  pathway <- c()
  protein <- c()
  rsd <- c()
  P_phos_level <- c(); P_kin_level <- c(); P_sub_level <- c()
  coef_phos_level <- c(); coef_kin_level <- c(); coef_sub_level <- c()
  fdr_phos_level <- c(); fdr_kin_level <- c(); fdr_sub_level <- c()
  
  count <- 0
  for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
    #for (path in "ErbB signaling pathway") {
    pho_path_all <- pho_pathway_table[pho_pathway_table$Pathway==path,]
    rownames(pho_path_all) <- pho_path_all$Sample
    pho_path <- pho_path_all[sam_subtype,]
    Global_phosphorylation <- pho_path$Global_phosphorylation
    genes <- KEGG[[which(grepl(path,kegg))]]
    genes_overlap <- intersect(genes,unique(pho_rsd_split$SUBSTRATE))
    
    for (g in genes_overlap) {
      rsd_rows <- which(pho_rsd_split$SUBSTRATE==g)
      scores_all <- data.frame(overlap[overlap$kinase==g,c("sample","up_phos_c","sub_phos_ctrans")])
      rownames(scores_all) <- scores_all$sample
      scores <- scores_all[sam_subtype,]
      
      for (i in rsd_rows) {
        pho_site <- as.numeric(pho_data[i,sam_subtype])
        pho_gene <- data.frame(scores, pho_site, Global_phosphorylation)
        
        if (length(which(complete.cases(pho_gene[,c("pho_site","up_phos_c","sub_phos_ctrans","Global_phosphorylation")]))) >= least_samples) {
          fit <- glm(Global_phosphorylation ~ pho_site + up_phos_c + sub_phos_ctrans, data = pho_gene); success <- 1
        } else if (length(which(complete.cases(pho_gene[,c("pho_site","sub_phos_ctrans","Global_phosphorylation")]))) >= least_samples) {
          fit <- glm(Global_phosphorylation ~ pho_site + sub_phos_ctrans, data = pho_gene); success <- 1
        } else if (length(which(complete.cases(pho_gene[,c("pho_site","up_phos_c","Global_phosphorylation")]))) >= least_samples) {
          fit <- glm(Global_phosphorylation ~ pho_site + up_phos_c, data = pho_gene); success <- 1
        } else if (length(which(complete.cases(pho_gene[,c("pho_site","Global_phosphorylation")]))) >= least_samples) {
          fit <- glm(Global_phosphorylation ~ pho_site, data = pho_gene); success <- 1
        } else {
          success <- 0
        }
        if (success == 1) {
          count <- count + 1
          pathway[count] <- path
          protein[count] <- g
          rsd[count] <- pho_rsd_split$SUB_MOD_RSD[i]
          
          pvalues <- coef(summary(fit))
          
          if ( length(which(rownames(pvalues)=="pho_site")) > 0 ) {
            P_phos_level[count] <- pvalues["pho_site",4]
          } else {
            P_phos_level[count] <- NA
          }
          if ( length(which(rownames(pvalues)=="up_phos_c")) > 0 ) {
            P_kin_level[count] <- pvalues["up_phos_c",4]
          } else {
            P_kin_level[count] <- NA
          }
          if ( length(which(rownames(pvalues)=="sub_phos_ctrans")) > 0 ) {
            P_sub_level[count] <- pvalues["sub_phos_ctrans",4]
          } else {
            P_sub_level[count] <- NA
          }
          coef_phos_level[count] <- fit$coefficients["pho_site"];
          coef_kin_level[count] <- fit$coefficients["up_phos_c"]; 
          coef_sub_level[count] <- fit$coefficients["sub_phos_ctrans"]
        }
      }
    }
    fdr_phos_level <- c(fdr_phos_level,p.adjust(P_phos_level[pathway==path], method = "fdr"))
    fdr_kin_level <- c(fdr_kin_level,p.adjust(P_kin_level[pathway==path], method = "fdr"))
    fdr_sub_level <- c(fdr_sub_level,p.adjust(P_sub_level[pathway==path], method = "fdr"))
  }
  pathway_site_kinsub_cor <- data.frame(pathway,protein,rsd,
                                        fdr_phos_level,fdr_kin_level,fdr_sub_level,
                                        coef_phos_level,coef_kin_level,coef_sub_level,
                                        P_phos_level,P_kin_level,P_sub_level)
  pathway_site_kinsub_cor$subtype <- subtype
  pathway_site_kinsub_cor_subtypes <- rbind(pathway_site_kinsub_cor_subtypes, pathway_site_kinsub_cor)
}

pathway_site_kinsub_cor_sub_ordr <- c()
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  temp <- pathway_site_kinsub_cor_subtypes[pathway_site_kinsub_cor_subtypes$pathway==path,]
  temp <- temp[order(temp$fdr_phos_level),]
  pathway_site_kinsub_cor_sub_ordr <- rbind(pathway_site_kinsub_cor_sub_ordr,temp)
}

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/pathway_site_kinsub_phos_reg_by_subtypes_adjust_pvalues_within_pathways_BRCA.txt", sep="")
write.table(pathway_site_kinsub_cor_sub_ordr, file=tn, quote=F, sep = '\t', row.names = FALSE)

# pathway_pho&phophosite_pho_corr, divided by subtypes -----------------------------------------
clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))
pho_pathway_table <- data.frame(BRCA_Pho_pathway_all_merge)

pathway_site_pho_cor_subtypes <- c()
for (subtype in c("Her2")) {
  sam_subtype <- colnames(clinical)[clinical[1,]==subtype & !is.na(clinical[1,])]
  
  pathway <- c()
  protein <- c()
  rsd <- c()
  pvalue <- c()
  cor <- c()
  fdr <- c()
  count <- 0
  #for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  for (path in "ErbB signaling pathway") {
    pho_path_all <- pho_pathway_table[pho_pathway_table$Pathway==path,]
    rownames(pho_path_all) <- pho_path_all$Sample
    pho_path <- pho_path_all[sam_subtype,]
    genes <- KEGG[[which(grepl(path,kegg))]]
    genes_overlap <- intersect(genes,unique(pho_rsd_split$SUBSTRATE))
    for (g in genes_overlap) {
      rsd_rows <- which(pho_rsd_split$SUBSTRATE==g)
      for (i in rsd_rows) {
        corr_stat = try(cor.test(pho_path$Global_phosphorylation,
                                 as.numeric(pho_data[i,sam_subtype]),
                                 method = "pearson"), silent=T)
        if (!is(corr_stat,"try-error")){ # if the correlation test is carried out successfully
          count <- count + 1
          pathway[count] <- path
          protein[count] <- g
          rsd[count] <- pho_rsd_split$SUB_MOD_RSD[i]
          cor[count] <- corr_stat$estimate
          pvalue[count] <- corr_stat$p.value
        }
      }
    }
    fdr <- c(fdr,p.adjust(pvalue[pathway==path], method = "fdr"))
  }
  pathway_site_pho_cor <- data.frame(pathway,protein,rsd,fdr,cor,pvalue)
  pathway_site_pho_cor$subtype <- subtype
  pathway_site_pho_cor_subtypes <- rbind(pathway_site_pho_cor_subtypes, pathway_site_pho_cor)
}

pathway_site_pho_cor_sub_ordr <- c()
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway)) {
  temp <- pathway_site_pho_cor_subtypes[pathway_site_pho_cor_subtypes$pathway==path,]
  temp <- temp[order(temp$pvalue),]
  pathway_site_pho_cor_sub_ordr <- rbind(pathway_site_pho_cor_sub_ordr,temp)
}
pathway_site_pho_cor_sub_ordr$cor_type[pathway_site_pho_cor_sub_ordr$cor>0] <- "positively-correlated"
pathway_site_pho_cor_sub_ordr$cor_type[pathway_site_pho_cor_sub_ordr$cor<0] <- "negatively-correlated"

