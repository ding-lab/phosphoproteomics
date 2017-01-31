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
mod <- "pho_sub~pro_kin"
table <- rbind(table_kinase,table_photase)

# add rtk variable
rtk <- c("ERBB3","ARAF","AURKA","EGFR","FGF4","FGFR1","FGFR2","FGFR3","FGFR4","IGF1R","KIT","NTRK1","PDGFRA","PDGFRB","PRKCE","RET","SRC","SYK","EPHA3","ERBB4")
table$rtk <- vector(mode = "logical", length = nrow(table))

for( kinase in rtk) {
  temp <- which(table$KINASE==kinase)
  if ( length(temp) > 0 ) {
    table$rtk[temp] <- TRUE
  }
}

# write intergrated table
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/kinase_phosphotase_regression.txt", sep="")
write.table(table, file=tn, quote=F, sep = '\t', row.names = FALSE)

# make the table you want to write ----------------------------------------
self <- TRUE
sig <- 0.1
table <- table[table$model == mod & table$FDR_pro_kin <= sig & table$self==self,]
table$scale_FDR <- -log10(table$FDR_pro_kin)

# write results -----------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/self_",self,"_fdr_",sig,".txt", sep="")
write.table(table, file=tn, quote=F, sep = '\t', row.names = FALSE)




