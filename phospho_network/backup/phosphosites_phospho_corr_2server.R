##### cross_lvl_corr.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# cross level correlation between CNV, RNA, and Proteome data

# directory and library ---------------------------------------------------
install.packages("stringr")
library(stringr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer/"

# input regardless to cancer type, choose between kinase or phosphotase-------------------------------------------------------------------
### read in the kinase/substrate table/ phosphorylation data ###
K_S_f = paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep="")
k_s_table = read.delim(K_S_f)

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

unique_kinase <- unique(k_s_table$GENE)
unique_substrate <- unique(k_s_table$SUB_GENE)
least_samples <- 5

# input -------------------------------------------------------------------
cancer <- "BRCA"
BRCA_pho_f = paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv",sep="")
pho_data = read.delim(BRCA_pho_f)

# ordering the columns by sample name
pho_data <- pho_data[,order(names(pho_data))]
pho_main <- pho_data[,1:77]

#split the SUBSTRATE and SUB_MOD_RSD in the first column
pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))

#covert the SUB_MOD_RSD from lowercase to uppercase
pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")


# initiate ----------------------------------------------------------------
## substrate, rsd1,rsd2, pvalue, coef, genomic.pos1, genomic.pos2, 
template <- data.frame(t(vector(mode = "numeric", length = 9 )))
colnames(template) <- c("SUBSTRATE","RSD1","RSD2","fdr","pvalue","coef","pos1","pos2","size")
phosite_corr <- c()

# loop around each protein ------------------------------------------------
start.time <- Sys.time()
for (gene in unique_substrate) {
  #for (gene in "EGFR") {
  rows <- which(pho_rsd_split$SUBSTRATE==gene)
  n <- length(rows)
  if (n > 1) {
    template1 <- template
    template1$SUBSTRATE <- gene
    
    for (i in rows[-n]) {
      rsd1 <- pho_rsd_split$SUB_MOD_RSD[i];#pos1 <- pho_data$genomic_pos[i]
      pho_temp1 <- pho_main[i,]
      pho_temp1_norm <- range01(unlist(pho_temp1),na.rm = T)
      
      template2 <- template1
      template2$RSD1 <- rsd1;#template2$pos1 <- pos1
      
      for (j in (i+1):(max(rows)) ) {
        rsd2 <- pho_rsd_split$SUB_MOD_RSD[j];#pos2 <- pho_data$genomic_pos[j]
        
        pho_temp2 <- pho_main[j,]
        pho_temp2_norm <- range01(unlist(pho_temp2),na.rm = T)
        
        #prepare regression data for model1
        data1 <- data.frame(pho_temp1_norm,pho_temp2_norm)
        
        size <- nrow(data1[complete.cases(data1),])
        if( size > least_samples ){#more than 2 complete dataset
          temp <- template2
          temp$RSD2 <- rsd2; #temp$pos2 <- pos2;
          temp$size <- size
          fit1 <- glm(pho_temp1_norm ~ pho_temp2_norm,data = data1, family=gaussian())
          
          temp$pvalue <- c(coef(summary(fit1))[2,4])
          temp$coef <- fit1$coefficients[2]
          phosite_corr <- rbind(phosite_corr,temp)
        }
      }
    }
  }
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# adjust for FDR ----------------------------------------------------------
phosite_corr <- phosite_corr[phosite_corr$pos1 != as.character(phosite_corr$pos2),]
phosite_corr$fdr <-p.adjust(phosite_corr$pvalue,method = "fdr")
phosite_corr$pair <- paste(phosite_corr$SUBSTRATE,phosite_corr$RSD1,phosite_corr$RSD2,sep = ":")
phosite_corr$fdr_log10 <- -log10(phosite_corr$fdr)

## write out tables
tn = "/gscuser/ayigewu/regression/kinase_substrate_phosphosite_correlation.txt"
write.table(phosite_corr, file=tn, quote=F, sep = '\t', row.names = FALSE)



