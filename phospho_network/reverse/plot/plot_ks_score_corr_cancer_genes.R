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
overlap <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/kinase_BRCA_reverse_marking_median_scored.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
ks_phos_corr <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/",protein,"_",cancer,"_kin:sub-based_phos_score_correlate_protein_phos_score.txt", sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
cancer_genes <- read_delim(paste(baseD,"pan3can_shared_data/reference_files/Volgestin2013Science_125genes.txt", sep=""), 
                           "\t", escape_double = FALSE, col_types = cols(`Ocogene score*` = col_character()), 
                           trim_ws = TRUE)
gene_list <- cancer_genes
rows <- c()
for (gene in cancer_genes$`Gene Symbol`) {
  temp <- which(ks_phos_corr$gene==gene)
  rows <- c(rows,temp)
}
ks_phos_corr4list <- ks_phos_corr[rows,]
ks_phos_corr4list_sig <- ks_phos_corr4list[(!is.na(ks_phos_corr4list$s_fdr) & ks_phos_corr4list$s_fdr<=sig) | (!is.na(ks_phos_corr4list$k_fdr) & ks_phos_corr4list$k_fdr<=sig),]

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/cancer_genes_fdr_",sig,"_corr_ks-based_score4",cancer,"_from_",protein,"table.txt", sep="")
write.table(ks_phos_corr4list_sig, file=tn, quote=F, sep = '\t', row.names = FALSE)


# scatterplot the significant ones for k_score----------------------------------------
gene4k <- as.vector(ks_phos_corr4list$gene[!is.na(ks_phos_corr4list$k_fdr) & ks_phos_corr4list$k_fdr<=sig])
rows <- c()
for (gene in gene4k) {
  temp <- which(overlap$kinase==gene)
  rows <- c(rows,temp)
}
tablek <- overlap[rows,]
lim = max(abs(max(tablek$sub_phos_ctransscore, na.rm = T)),abs(min(tablek$sub_phos_ctransscore, na.rm = T)))
p = ggplot(data = tablek,aes(x= up_phos_cscore, y = kin_phos_cscore, shape = mut, colour = sub_phos_ctransscore , alpha = 1))
p = p + scale_colour_gradientn(name= "sub_phos_ctransscore", na.value="grey40", colours=getPalette(100)) # , limit=c(-lim,lim)
p = p + geom_point(size = 2, alpha=0.8, stroke = 0) 
p = p + geom_text(data=subset(tablek, aa != "wt" ),
                  aes(x= up_phos_cscore, y = kin_phos_cscore,label=aa), hjust = 0, nudge_x = -0.5 ,nudge_y = 0.2, size = 2) +
  facet_wrap(~ kinase, nrow = 4)
p <- p + xlab("kinase-based phosphorylation score")+ylab("protein phosphorylation score")
p = p + theme_bw()
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/',cancer,'_cancer_genes_kin-based_score~protein_phos_score_fdr_',sig,'.pdf',sep ="")
ggsave(file=fn, height=8, width=8)

# scatterplot the significant ones for s_score----------------------------------------
gene4s <- as.vector(ks_phos_corr4list$gene[!is.na(ks_phos_corr4list$s_fdr) & ks_phos_corr4list$s_fdr<=sig])
rows <- c()
for (gene in gene4s) {
  temp <- which(overlap$kinase==gene)
  rows <- c(rows,temp)
}
tables <- overlap[rows,]


lim = max(abs(max(tables$up_phos_cscore, na.rm = T)),abs(min(tables$up_phos_cscore, na.rm = T)))
p = ggplot(data = tables,aes(x= sub_phos_ctransscore , y = kin_phos_cscore, shape = mut, colour = up_phos_cscore , alpha = 1))
p = p + scale_colour_gradientn(name= "up_phos_cscore", na.value="grey40", colours=getPalette(100)) # , limit=c(-lim,lim)
p = p + geom_point(size = 2, alpha=0.8, stroke = 0) 
# p = p + geom_text(data=subset(tables, aa != "wt" & !is.na(aa) ),
#                   aes(x= sub_phos_ctransscore, y = kin_phos_cscore,label=aa), hjust = 0, nudge_x = -0.5 ,nudge_y = 0.2, size = 2)
p = p + facet_wrap(~ kinase, nrow = 3)
p <- p + xlab("substrate-based phosphorylation score")+ylab("protein phosphorylation score")
p = p + theme_bw()
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/',cancer,'_cancer_genes_sub-based_score~protein_phos_score_fdr_',sig,'.pdf',sep ="")
ggsave(file=fn, height=4, width=4)

