# Yige Wu @ WashU 2017 Jan
# overlap regression result with PDB active sites

# directory and library ---------------------------------------------------
# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"

library(stringr)
library(ggplot2)
library(readr)

setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# choose kinase or phosphotase, significance level, outlier threshold and least sample number-------------------------
sig <- 0.05 # significance level
out_thres <- 1.5 #threshold for outlier
least_samples <- 5# least number of samples with complete data for each model
protein <- "kinase"


# input PDB file and processed regression result ---------------------------------------
PDB_site <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/PDB_site.maf",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/kinase_substrate_regression.txt",sep = ""))
table1 <- table_2can[table_2can$model=="pho_sub~pro_kin",]
rsd_list <- str_split_fixed(PDB_site$amino_acid_change,"[p.XZ]",5)[,4]
sub_list <- as.vector(PDB_site$Hugo_Symbol)
# bubble chart, appoint phosphosites ---------------------
for (cis in c("cis","trans")) { # loop around cis and trans
  t0 <- table1[table1$SELF==cis & table1$FDR_pro_kin <= sig,]
  t0_rsd <- str_split_fixed(t0$SUB_MOD_RSD,"[STY]",3)[,2]
  rows <- c()
  for(i in 1:nrow(t0)){
    r <- unlist(which(t0$SUBSTRATE[i]==sub_list & t0_rsd[i] ==rsd_list))
    if (length(r) > 0) {
      rows <- c(rows,i)
    }
  }
  table <- t0[rows,]
  if (nrow(table) > 0){
    print((length(rows))/nrow(t0))
    table$sig <- (table$FDR_pro_kin <= sig)
    ## actual plotting
    lim = max(max(table$coef_pro_kin),min(table$coef_pro_kin))
    p = ggplot(table,aes(x=SELF, y=pair))# make this the original ethni
    #  p = p + facet_grid(KINASE~., drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
    p = p + facet_grid(KINASE~Cancer, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
    p = p + geom_point(aes(fill=coef_pro_kin, size =-log10(FDR_pro_kin), color=ifelse(sig, "black",NA)),pch=21) 
    p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
    p = p + scale_colour_manual(values=c("black",NA))
    p = p + theme_bw() #+ theme_nogrid()
    p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
    p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
    p
    fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/regression_PDB_overlap_',cis,'_sig_',sig,'.pdf',sep ="")
    ggsave(file=fn, height=22, width=5.5, useDingbats=FALSE)
  }
}


# prepare the pairwise file and PDB file ----------------------------------
pairwise <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/PDB_cptac.brca.ov_sites.pairwise.singleprotein.collapsed",sep = ""), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
pairwise <- data.frame(pairwise[,c(1,2,6,7,12,13,14,15)])
colnames(pairwise) <- c("Gene1","RSD1","Gene2","RSD2","dis_linear","dis_3d","PDB","pvalue_3d")
rsd1_pair <-  str_split_fixed(pairwise$RSD1,"[p.XSTYCZ]",5)[,4]
rsd2_pair <-  str_split_fixed(pairwise$RSD2,"[p.XSTYCZ]",5)[,4]
pairwise$dis_lin <- abs(as.numeric(rsd2_pair)-as.numeric(rsd1_pair))
gene1_pair <- as.vector(pairwise$Gene1)

rsd_pdb <- str_split_fixed(PDB_site$amino_acid_change,"[p.XSTYCZ]",5)[,4]
cat_pdb <- str_split_fixed(PDB_site$Tumor_Sample_Barcode,":",2)[,1]

# mark the active sites for PDB_site and mark pairwise data ---------------
rsd1_cat <- vector(mode = "character", length = nrow(pairwise)); rsd2_cat <- rsd1_cat
for( i in 1:nrow(pairwise) ) {
  temp <- cat_pdb[PDB_site$Hugo_Symbol==gene1_pair[i] & rsd_pdb==rsd1_pair[i]]
  if (length(temp) > 0) {
    rsd1_cat[i] <- temp
  }
  temp <- cat_pdb[PDB_site$Hugo_Symbol==gene1_pair[i] & rsd_pdb==rsd2_pair[i]]
  if (length(temp) > 0) {
    rsd2_cat[i] <- temp
  }
}
pairwise$rsd1_cat <- rsd1_cat; pairwise$rsd2_cat <- rsd2_cat

# mark the pariwise data with regression result: significant/not significant/unknown ---------
# significant means having at least one significant regression result
table <- table_2can[table_2can$model=="pho_sub~pro_kin" & table_2can$Cancer=="BRCA",]
rsd_reg <- str_split_fixed(table$SUB_MOD_RSD,"[STY]",3)[,2]
rsd1_sig <- vector(mode = "numeric", length = nrow(pairwise)); rsd2_sig <- rsd1_sig
for( i in 1:nrow(pairwise) ) {
  temp1 <- which((table$SUBSTRATE==gene1_pair[i]) & (rsd_reg==rsd1_pair[i]) & (table$FDR_pro_kin <= sig) )
  rsd1_sig[i] <- length(temp1)
  
  temp2 <- which((table$SUBSTRATE==gene1_pair[i]) & (rsd_reg==rsd2_pair[i]) & (table$FDR_pro_kin <= sig) )
  rsd2_sig[i] <- length(temp2)
}
pairwise$rsd1_sig <- rsd1_sig; pairwise$rsd2_sig <- rsd2_sig

length(unique(pairwise$Gene1[rsd1_cat=="ACT_SITE" | rsd2_cat=="ACT_SITE"]))
length(unique(pairwise$Gene1[rsd1_sig > 0 | rsd2_sig >0]))

# look at the distance between significant phosphosites and active sites --------
pairwise$dis_sig <- (pairwise$rsd1_cat=="ACT_SITE" & pairwise$rsd2_sig >0 ) | (pairwise$rsd2_cat=="ACT_SITE" & pairwise$rsd1_sig >0 )
pro_sig <- unique(pairwise$Gene1[pairwise$dis_sig])
pairwise_act <- c()
for( gene in pro_sig ) {
  pairwise_act <- rbind(pairwise_act, pairwise[pairwise$Gene1==gene,])
}
tn = paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/BRCA_pairwise_distance_with_active_site_and_sig_",sig,"_regression.txt", sep="")
write.table(pairwise_act, file=tn, quote=F, sep = '\t', row.names = FALSE)
