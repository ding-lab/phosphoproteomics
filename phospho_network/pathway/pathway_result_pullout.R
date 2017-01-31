## pull out all the relevant regression result of a certain pathway

# load data for kinase and phosphotase ------------------------------------
baseD <-"/Users/yigewu/Box Sync/"
load(paste(baseD,"pan3can_shared_data/analysis_results/regression/Rdata/phosphotase.RData",sep = ""))
table_photase <- table_2can
table_photase$protein <- "phosphotase"

load(paste(baseD,"pan3can_shared_data/analysis_results/regression/Rdata/kinase.RData",sep = ""))
table_kinase <- table_2can
table_kinase$protein <- "kinase"
table_kinase <- table_kinase[,colnames(table_photase)]


# integrate kinase & phosphotase data and shrink the gene list ------------
table <- rbind(table_kinase,table_photase)
mod <- "pho_sub~pro_kin"
sig <- 0.1
table <- table[table$model==mod & table$FDR_pro_kin <= sig,]
table$scale_FDR <- -log10(table$FDR_pro_kin)
kinase_list <- as.vector(unique(table$KINASE))
substrate_list <- as.vector(unique(table$SUBSTRATE))

# loop around pathways ----------------------------------------------------
kegg <- names(KEGG)
kegg_name <- gsub("\t", "_", kegg, fixed = TRUE)
# initiate the gene list for the pathway ------------------------------------------------------
for (i in 1:length(kegg) ) {
  path <- kegg[i]
  genes <- KEGG[[path]]
  
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
    tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/tables/",kegg_name[i],"_self_",self,"_fdr_",sig,".txt", sep="")
    write.table(table_temp, file=tn, quote=F, sep = '\t', row.names = FALSE)
  }
}

