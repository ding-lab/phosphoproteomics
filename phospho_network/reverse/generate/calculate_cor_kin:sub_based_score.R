# Yige Wu @ WashU 2017 Feb
# do correlation between kinase-based phosphorylation score/substrate-based phosphorylation score and protein phosphorylation score
# I think we should follow the model 

# choose kinase or phosphotase, outlier threshold and least sample number-------------------------
protein <- "kinase"
cancer <- "BRCA"

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

# input overlap file ------------------------------------------------------
#overlap <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/kinase_BRCA_reverse_marking_out_thers_1.5.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
# overlap <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/kinase_BRCA_reverse_marking_median_out_thers_1.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
overlap <- read_delim("~/Box Sync/pan3can_shared_data/analysis_results/tables/BRCA_ks_score_validated_ks_only_median_scored_kinase_table.txt", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)

# initiate data.frame put in the cor.stat ---------------------------------
gene_list <- unique(overlap$kinase)
k_estimate <- vector(mode = "numeric", length = length(gene_list)) + NA
k_pvalue <- k_estimate;s_estimate <- k_estimate; s_pvalue <- k_estimate

# loop over genes in overlap file -----------------------------------------
for (i in 1:length(gene_list)) {
  gene <- gene_list[i]
  g_table <- overlap[overlap$kinase==gene,]
  # test upstream kinase phos correlation
  corr_statk = try(cor.test(g_table$up_phos_c,g_table$kin_phos_c, method = "pearson"), silent=T)
  if (!is(corr_statk,"try-error")){ # if the correlation test is carried out successfully
    k_estimate[i] <- corr_statk$estimate
    k_pvalue[i] <- corr_statk$p.value
  }
  # test downstream substrate phos correlation
  corr_stats = try(cor.test(g_table$sub_phos_ctrans,g_table$kin_phos_c, method = "pearson"), silent=T)
  if (!is(corr_stats,"try-error")){ # if the correlation test is carried out successfully
    s_estimate[i] <- corr_stats$estimate
    s_pvalue[i] <- corr_stats$p.value
  }
}

# integrate & adjust pvalue-----------------------------------------------
ks_phos_corr <- data.frame(gene_list,k_estimate,k_pvalue,s_estimate,s_pvalue)
colnames(ks_phos_corr)[1] <- "gene"

ks_phos_corr$k_fdr <- p.adjust(ks_phos_corr$k_pvalue,method = "fdr")
ks_phos_corr$s_fdr <- p.adjust(ks_phos_corr$s_pvalue,method = "fdr")
num_ks <- data.frame(unique(overlap[,c("kinase","n_up_pho","n_sub_pho")]))
colnames(num_ks) <- c("gene","num_k","num_s")
ks_phos_corr <- merge(ks_phos_corr,num_ks, all.x = T)

# write table -------------------------------------------------------------
# tn = paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/",protein,"_",cancer,"_kin:sub-based_phos_score_correlate_protein_phos_score.txt", sep="")
# tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_",cancer,"_kin:sub-based_phos_score_correlate_protein_phos_score.txt", sep="")
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_",cancer,"_validated_kin:sub-based_phos_score_correlate_protein_phos_score.txt", sep="")
write.table(ks_phos_corr, file=tn, quote=F, sep = '\t', row.names = FALSE)

ks_phos_corr_sig <- ks_phos_corr[(!is.na(ks_phos_corr$s_fdr) & ks_phos_corr$s_fdr<=sig) | (!is.na(ks_phos_corr$k_fdr) & ks_phos_corr$k_fdr<=sig),]


# examine intersted proteins --------------------------------
k_s_table = read.delim(paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
gene_list <- c("ABL2","BCR","BRAF","BRCA2","CDK1","CDK13","CDK7","CHEK2","CHEK1","EGFR","MAPK1","MAPK3","NUP153","PAK2","PLK1","PTK2","RBL1","STAT3")
time <- Sys.time()
sink(file = paste(baseD,"pan3can_shared_data/analysis_results/tables/ks_composition4ks_score_sig_correlated_cancer_genes",time,".txt",sep = ""))
for(g in gene_list) {
  temp <- ks_phos_corr[ks_phos_corr$gene==g,]
  cat(paste(g, ":\n", 
            "     kinase-based score: coef = ",temp$k_estimate,", fdr = ",temp$k_fdr,", ",temp$num_k," kinase(s) among :\n",sep = ""))
  cat("                                                                                                 ")
  cat(paste(as.vector(unique(k_s_table$GENE[k_s_table$SUB_GENE==g & k_s_table$GENE!=g])), collapse = " "))
  cat(paste("\n\n","     substrate-based score: coef = ",temp$s_estimate,", fdr = ",temp$s_fdr,", ",temp$num_s," substrates(s) among :\n",sep = ""))
  cat("                                                                                                 ")
  cat(paste(as.vector(unique(k_s_table$SUB_GENE[k_s_table$GENE==g &  k_s_table$SUB_GENE!=g])), collapse = " "))
  cat("\n\n")
}
sink()


