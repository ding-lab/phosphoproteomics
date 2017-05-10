# adopted from Yige Wu @ WashU 2016 Dec
# Kuan @wustl 2017
# barplot the validation statistics for regression results

# choose kinase/phosphotase, cancer , significance level -----------------------------------------------
protein <- "kinase"
sig <- 0.05
# cancer <- "BRCA"
cancer <- "OV"


# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)
library(readr)

# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

valid_count_thres=0

# input valid_trans -------------------------------------------------------
valid_trans_brca <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/BRCA_kinase_trans_regression_validation_statistics.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
valid_trans_brca$cancer <- "BRCA"
valid_trans_ov <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/OV_kinase_trans_regression_validation_statistics.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
valid_trans_ov$cancer <- "OV"
valid_trans <- data.frame(rbind(valid_trans_brca,valid_trans_ov))
valid_trans$none <- valid_trans$all_count-valid_trans$pos_sig

# barplot ranking validation count for trans----------------------------------------
top <- 10
valid_trans0 <- valid_trans[!is.na(valid_trans$valid_ratio) & valid_trans$valid_count > 1 ,c("kinase","valid_count","cancer","pos_sig","none")]
valid_trans0 <- valid_trans0[order(valid_trans0$valid_count, decreasing = T),]
valid_trans0brca <- valid_trans0[valid_trans0$cancer=="BRCA",]; valid_trans0ov <- valid_trans0[valid_trans0$cancer=="OV",]
valid_trans_top <- rbind(valid_trans0brca[1:(min(top,nrow(valid_trans0brca))),],valid_trans0ov[1:(min(top,nrow(valid_trans0ov))),])
table_trans <- melt(valid_trans_top,id=c("kinase","valid_count","cancer"))
colnames(table_trans) <- c("kinase","valid_count","cancer","coef_FDR","count")
table_trans$KINASE <- reorder(table_trans$kinase, table_trans$valid_count)
table_trans$coef_FDR = as.character(table_trans$coef_FDR)
table_trans$SELF <- "trans"

p <- ggplot()
p <- p + geom_bar(data=table_trans, aes(y = count, x = KINASE, fill = coef_FDR ), stat="identity",
                  position='stack')
p <- p + facet_grid(SELF~cancer, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p <- p + scale_fill_manual(values = c("pos_sig" = "red","none" = "grey"))
p <- p + theme_bw() + scale_y_log10()
p <- p + xlab(protein)+ylab("number of substrate phosphosites")
p <- p + theme(axis.title=element_text(size=10))
p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/2can_',protein,'_trans_validation_count_ranking_per_gene_top_least_',valid_count_thres,'valid.pdf',sep ="")
ggsave(file=fn, height=3, width=4)

# input valid_cis ---------------------------------------------------------
valid_cis_brca <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/BRCA_kinase_cis_regression_validation_statistics.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
valid_cis_brca$cancer <- "BRCA"
valid_cis_ov <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/OV_kinase_cis_regression_validation_statistics.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
valid_cis_ov$cancer <- "OV"
valid_cis <- data.frame(rbind(valid_cis_brca,valid_cis_ov))
valid_cis$none <- valid_cis$all_count-valid_cis$pos_sig

# barplot ranking validation count for cis----------------------------------------
temp <- valid_cis[!is.na(valid_cis$valid_ratio),c("kinase","valid_count","cancer","pos_sig","none")]
table_cis <- melt(temp,id=c("kinase","valid_count","cancer"))
colnames(table_cis) <- c("kinase","valid_count","cancer","coef_FDR","count")
table_cis$KINASE <- reorder(table_cis$kinase, table_cis$valid_count)
table_cis$coef_FDR = as.character(table_cis$coef_FDR)
table_cis_top = table_cis[(table_cis$cancer=="BRCA" & table_cis$valid_count >= 8) | (table_cis$cancer=="OV" & table_cis$valid_count >= 2),]
table_cis_top$log_count = log10(table_cis_top$count)
table_cis_top$log_count[is.infinite(table_cis_top$log_count)] = 0
table_cis_top$SELF <- "cis"

p <- ggplot()
p <- p + geom_bar(data=table_cis_top, aes(y = count, x = KINASE, fill = coef_FDR ), stat="identity",
                  position='stack')
p <- p + facet_grid(SELF~cancer, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p <- p + scale_fill_manual(values = c("pos_sig" = "red","none" = "grey"))
p <- p + theme_bw()
p <- p + xlab(protein)+ylab("number of substrate phosphosites")
p <- p + theme(axis.title=element_text(size=10))
#p <- p + coord_flip() 
p <- p + theme(axis.text.x = element_text(colour="black", size=8,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/2can_',protein,'_cis_validation_count_ranking_per_gene_top_least_',valid_count_thres,'valid.pdf',sep ="")
ggsave(file=fn, height=3, width=4)

