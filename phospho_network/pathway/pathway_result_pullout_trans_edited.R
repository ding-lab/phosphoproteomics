## pull out all the relevant regression result of a certain pathway

# load data for kinase and phosphotase ------------------------------------
baseD <-"/Users/yigewu/Box Sync/"
sig <- 0.05

table_phosphotase <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/phosphotase_substrate_regression_trans_edited.txt",sep = ""))
table_phosphotase$protein <- "phosphotase"

table_kinase <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/kinase_substrate_regression_trans_edited.txt",sep = ""))
table_kinase$protein <- "kinase"

load("~/Box Sync/pan3can_shared_data/analysis_results/2015-08-01_Gene_Set.RData")
BRCA_Pho_pathway_all_merge <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/pathway_activation/2017-02-15_KH_BRCA_cross_pathway_activation.tsv",sep = "") , 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)

# integrate kinase & phosphotase data and shrink the gene list ------------
table <- rbind(table_kinase,table_phosphotase)
scale_FDR <- vector(mode = "numeric", length = nrow(table))
scale_FDR[table$self] <- -log10(table$FDR_pro_kin[table$self])
scale_FDR[!table$self] <- -log10(table$FDR_pho_kin[!table$self])
table$scale_FDR <- scale_FDR

table$sig <- (table$self & table$FDR_pro_kin <= sig ) | (!table$self & table$FDR_pho_kin <= sig)
table$coef_pos <- (table$self & table$coef_pro_kin > 0 ) | (!table$self & table$coef_pho_kin > 0 )
kinase_list <- as.vector(unique(table$KINASE))
substrate_list <- as.vector(unique(table$SUBSTRATE))

# loop around KEGG pathways,output all the trans ----------------------------------------------------
kegg <- names(KEGG)
#for (i in 1:length(kegg) ) {
for (path in unique(BRCA_Pho_pathway_all_merge$Pathway) ) {
  genes <- KEGG[[which(grepl(path,kegg))]]
 
  gene_kin <- intersect(kinase_list,genes)
  gene_sub <- intersect(substrate_list,genes)
  
  # loop around the genes in the pathway ------------------------------------
  table_path <- c()
  for (k in gene_kin) {
    for (s in gene_sub) {
      temp <- table[table$KINASE==k & table$SUBSTRATE==s & !table$self,c("KINASE","SUBSTRATE","SUB_MOD_RSD"
                                                                         ,"FDR_pro_sub","FDR_pho_kin","coef_pro_sub","coef_pho_kin",
                                                                         "Cancer","protein","scale_FDR","sig","coef_pos")]
      # for each kinase:substrate pair, keep all significant records but only keep one record with least FDR of insignificant records
      if (nrow(temp) > 0) {
        temp_sig <- temp[temp$sig,]
        temp_insig <- temp[!temp$sig,]; temp_insig <- temp_insig[order(temp_insig$FDR_pho_kin),]
        table_path <- rbind(table_path,temp_sig,temp_insig[1,])
      }
    }
  }
  table_path <- table_path[!is.na(table_path$KINASE),]
  if ( !is.null(table_path) ) {
    tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/tables/KEGG/",path,"_trans_sig_",sig,".txt", sep="")
    write.table(table_path, file=tn, quote=F, sep = '\t', row.names = FALSE)
  }
}

# loop around REACT pathways ----------------------------------------------------
react <- names(REACT)
react_name <- gsub("\t", "_", react, fixed = TRUE)
for (i in 1:length(react) ) {
  path <- react[i]
  genes <- REACT[[path]]
  
  gene_kin <- intersect(kinase_list,genes)
  gene_sub <- intersect(substrate_list,genes)
  
  # loop around the genes in the pathway ------------------------------------
  table_path <- c()
  for (k in gene_kin) {
    for (s in gene_sub) {
      temp <- table[table$KINASE==k & table$SUBSTRATE==s,]
      table_path <- rbind(table_path,temp)
    }
  }
  for (self in c(TRUE, FALSE)) {
    table_temp <- table_path[table_path$self==self,]
    # write results -----------------------------------------------------------
    tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/tables/REACT/",react_name[i],"_self_",self,"_fdr_",sig,".txt", sep="")
    write.table(table_temp, file=tn, quote=F, sep = '\t', row.names = FALSE)
  }
}

# output all trans --------------------------------------------------------
ks_pairs <- unique(table[,c("KINASE","SUBSTRATE")])
gene_kin <- as.vector(ks_pairs$KINASE)
gene_sub <- as.vector(ks_pairs$SUBSTRATE)
table_path <- c()
for (i in 1:nrow(ks_pairs)) {
  k <- gene_kin[i]
  s <- gene_sub[i]
  temp <- table[table$KINASE==k & table$SUBSTRATE==s & !table$self,c("KINASE","SUBSTRATE","SUB_MOD_RSD"
                                                                     ,"FDR_pro_sub","FDR_pho_kin","coef_pro_sub","coef_pho_kin",
                                                                     "Cancer","protein","scale_FDR","sig","coef_pos")]
  # for each kinase:substrate pair, keep all significant records but only keep one record with least FDR of insignificant records
  if (nrow(temp) > 0) {
    temp_sig <- temp[temp$sig,]
    temp_insig <- temp[!temp$sig,]; temp_insig <- temp_insig[order(temp_insig$FDR_pho_kin),]
    table_path <- rbind(table_path,temp_sig,temp_insig[1,])
  }
}
table_path <- table_path[!is.na(table_path$KINASE),]
if ( !is.null(table_path) ) {
  tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/tables/KEGG/all_trans_sig_",sig,".txt", sep="")
  write.table(table_path, file=tn, quote=F, sep = '\t', row.names = FALSE)
}


