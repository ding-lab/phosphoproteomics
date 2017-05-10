# Yige Wu @ WashU 2017 Feb
# adopted by Kuan Huang @ Washu 2017 Apr
# look at correlations of kinase and downstream substrates phosphorylation status within BRCA subtypes

# choose kinase/phosphotase, significance level, outlier threshold and least sample number-------------------------
least_samples <- 10# least number of samples with complete data for each model
protein <- "kinase"
# protein <- "phosphotase"

# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)


#for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))
source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# input k_s_table according to kinase or phosphotase-------------------------------------------------------------------
if ( protein == "kinase" ) {
  ### read in the kinase/substrate table/ phosphorylation data ###
  k_s_table_phosphosite = read.delim(paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
  k_s_table_network = read.delim(paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphoNetworks/comKSI.csv",sep=""))
  k_s_table_phosphosite_sum = k_s_table_phosphosite[,c("GENE","SUB_GENE")]
  colnames(k_s_table_phosphosite_sum) = c("Kinase","Substrate")
  k_s_table_phosphosite_sum$Score=0
  k_s_table = rbind(k_s_table_phosphosite_sum,k_s_table_network)
}

if ( protein == "phosphotase" ) {
  ### read in the phosphotase/substrate table/ phosphorylation data ### 
  k_s_table <- read.csv(paste(baseD,"pan3can_shared_data/Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.csv",sep = ""))
  colnames(k_s_table) <- c("Phosphatase_UniProtAC_human","GENE","Substrate_UniProtAC_ref","SUB_GENE","Substrate_Type","DephosphoSite","BioassayType", "PubMed_ID_rev")
}

kinase_trans <- as.vector(unique(k_s_table$Kinase[as.vector(k_s_table$Kinase)!=as.vector(k_s_table$Substrate)]))
kinase_cis <- as.vector(unique(k_s_table$Kinase[as.vector(k_s_table$Kinase)==as.vector(k_s_table$Substrate)]))

HUMAN_pro_f = paste(baseD,"pan3can_shared_data/BRCA/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt",sep="")
pro_data <- read.delim(HUMAN_pro_f)
HUMAN_pho_f = paste(baseD,"pan3can_shared_data/BRCA/BRCA77_unimodal_phosphoproteome-ratio-norm_wGpos_cleaned.txt",sep="")
pho_data = read.delim(HUMAN_pho_f)
## read in grouped phosphorylation data!
HUMAN_pho_g = paste(baseD,"pan3can_shared_data/BRCA/BRCA77_unimodal_phosphoproteome-ratio-norm_collapsed.txt",sep="")
pho_gdata = read.delim(HUMAN_pho_g)
clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))

#split the SUBSTRATE and SUB_MOD_RSD in the first column
pho_rsd_split <- data.frame(str_split_fixed(pho_data$Gene.site, ":", 3))

#covert the SUB_MOD_RSD from lowercase to uppercase
pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")

# remove duplicate and identical phosphoyrlation levels for different transcript
# ps: they're not important genes
dup_pho <- data.frame(table(paste(pho_rsd_split$SUBSTRATE,pho_rsd_split$SUB_MOD_RSD,sep = ":")))
dup_pho <- dup_pho[dup_pho$Freq>1,]
dup_pho <- data.frame(str_split_fixed(dup_pho$Var1, ":", 2))
dup_pro <- as.vector(dup_pho$X1); dup_rsd <- as.vector(dup_pho$X2)

remove_rows <- c()
for (i in 1:nrow(dup_pho)) {
  dup_rows <- which(pho_rsd_split$SUBSTRATE==dup_pro[i] & pho_rsd_split$SUB_MOD_RSD==dup_rsd[i])
  remove_rows <- c(remove_rows,dup_rows[-1])
}
pho_rsd_split <- pho_rsd_split[-remove_rows,]
pho_data <- pho_data[-remove_rows,]
transcripts <- as.vector(pho_rsd_split$transcript)

# initiate ----------------------------------------------------------------
# calculate the length of cis table
ncis <- 0
for (gene in kinase_cis) {
  ncis <- ncis + length(which(pho_rsd_split$SUBSTRATE==gene))
}

# calculate the length of trans table
ntrans <- 0
for (gene in kinase_trans) {
  subs <- k_s_table$Substrate[k_s_table$Kinase==gene & k_s_table$Substrate!=gene]
  for ( sub in unique(subs)) {
    ntrans <- ntrans + length(which(pho_rsd_split$SUBSTRATE==sub))
  }
}

# looping subtype -----------------------------------------------------------
for (cohort in c("Basal","Her2","LumA","LumB")) { 
  subtype_sample <- as.vector(na.omit(colnames(clinical)[clinical[1,]==cohort]))
  pro_main <- pro_data[,gsub("[0-9][0-9]TCGA","01A",colnames(pro_data)) %in% subtype_sample]
  pho_main <- pho_data[,gsub("[0-9][0-9]TCGA","01A",colnames(pho_data)) %in% subtype_sample]
  pho_gmain <- pho_gdata[,gsub("[0-9][0-9]TCGA","01A",colnames(pho_gdata)) %in% subtype_sample]
  
  
  # looping over kinases for cis pairs -----------------------------------------------------------------
  # initiating the table for cis
  vec_char <- vector(mode = "character", length = ncis)
  vec_num <- vector(mode = "numeric", length = ncis) + NA
  KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
  FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
  coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
  Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
  Size <- vec_num;P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
  
  count <- 0
  for (kinase in kinase_cis){
    #for (kinase in "ERBB2"){#test
    # find protein expression level for the kinase
    pro_kin <- pro_main[pro_data$X == kinase,]
    
    substrate <- kinase
    if(nrow(pro_kin) != 0){
      s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
      
      # go through all the phosphosites
      for (i in s_pho_table) {
        
        # find phosphorylation level
        pho_sub <- pho_main[i,]
        
        #prepare regression data for model1
        data1 <- data.frame(t(rbind(pho_sub,pro_kin)))
        colnames(data1) <- c("pho_sub","pro_kin")
        
        size <- nrow(data1[complete.cases(data1),])
        if( size > least_samples ){
          fit1 <- lm(pho_sub ~ pro_kin,data = data1)
          
          count <- count + 1
          KINASE[count] <- kinase
          SUBSTRATE[count] <- substrate
          SUB_MOD_RSD[count] <- pho_rsd_split$SUB_MOD_RSD[i]
          transcript[count] <- transcripts[i]
          P_pro_kin[count] <- c(coef(summary(fit1))[2,4])
          coef_pro_kin[count] <- fit1$coefficients[2]
          Size[count] <- size
        }
      }
    }
  }
  
  table_cis <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,
                          FDR_pro_kin,FDR_pro_sub,FDR_pho_kin,
                          coef_pro_kin,coef_pro_sub,coef_pho_kin,
                          Cancer,transcript,model,Size,
                          P_pro_kin,P_pro_sub,P_pho_kin)
  
  # looping over kinases for trans pairs -----------------------------------------------------------------
  # initiating the table for trans
  vec_char <- vector(mode = "character", length = ntrans)
  vec_num <- vector(mode = "numeric", length = ntrans) + NA
  KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
  FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
  coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
  Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
  Size <- vec_num;P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
  
  count <- 0
  for (kinase in kinase_trans){
    pho_kinase_g <- pho_gmain[pho_gdata$X == kinase,]
    
    if ( nrow(pho_kinase_g) > 0 ){
      k_sub <- unique(k_s_table$Substrate[k_s_table$Kinase == kinase & k_s_table$Substrate != kinase])
      
      for (substrate in k_sub){# for each substrate for one kinase
        # find its phosphosites-row numbers
        s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
        
        # go through all the phosphosites
        for (i in s_pho_table) {
          
          # find phosphorylation level
          pho_sub <- pho_main[i,]
          
          # find substrate expressio level and normaize
          pro_sub <- pro_main[pro_data$X == substrate,]
          
          if(nrow(pro_sub) != 0){
            #prepare regression data for model2
            data2 <- data.frame(t(rbind(pho_sub,pro_sub,pho_kinase_g)))
            colnames(data2) <- c("pho_sub","pro_sub","pho_kinase_g")
            
            size <- nrow(data2[complete.cases(data2),])
            if(size > least_samples ){
              fit2 <- lm(pho_sub ~ pro_sub + pho_kinase_g, data = data2)
              
              count <- count + 1
              KINASE[count] <- kinase
              SUBSTRATE[count] <- substrate
              SUB_MOD_RSD[count] <- pho_rsd_split$SUB_MOD_RSD[i]
              transcript[count] <- transcripts[i]
              P_pro_sub[count] <- c(coef(summary(fit2))[2,4])
              coef_pro_sub[count] <- fit2$coefficients[2]
              P_pho_kin[count] <- c(coef(summary(fit2))[3,4])
              coef_pho_kin[count] <- fit2$coefficients[3]
              Size[count] <- size
            }
          }
        }
      }
    }
  }
  table_trans <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,
                            FDR_pro_kin,FDR_pro_sub,FDR_pho_kin,
                            coef_pro_kin,coef_pro_sub,coef_pho_kin,
                            Cancer,transcript,model,Size,
                            P_pro_kin,P_pro_sub,P_pho_kin)
  
  # integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
  table_cis$model <- "pho_sub~pro_kin"
  tablecis <- table_cis[!is.na(table_cis$P_pro_kin),]
  
  table_trans$model <- "pho_sub~pro_sub+pho_kin"
  tabletrans <- table_trans[!is.na(table_trans$P_pro_sub),]
  
  table <- rbind(tablecis,tabletrans)
  table$cohort <- cohort
  
  if ( cohort == "Basal" ) {
    table_Basal <- table
  } 
  if ( cohort == "Her2" ) {
    table_Her2 <- table
  } 
  if ( cohort == "LumA" ) {
    table_LumA <- table
  }
  if ( cohort == "LumB" ) {
    table_LumB <- table
  }
}



# combine tables -------------------------------------------------------
temp <- rbind(table_Basal,table_Her2,table_LumA,table_LumB)
table_2can <- unique(temp) ## because there're duplicate phophorylation levels for a residual
table_2can$self <- as.character(table_2can$KINASE) == as.character(table_2can$SUBSTRATE)
table_2can$SELF <- "trans"; table_2can$SELF[table_2can$self] <- "cis"
table_2can$Cancer <- "BRCA"
name = c("pro_kin","pro_sub","pho_kin")

## adjust p-values to FDR
for (cohort in c("Basal","Her2","LumA","LumB")) {
  #for (cancer in "BRCA") {
  for(self in c(TRUE,FALSE)) {
    for(coln in name) {#adjust pvalues for each variable
      row <- (table_2can$self==self) & (table_2can$cohort==cohort)
      table_2can[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_2can[row,paste("P_",coln,sep = "")],method = "fdr")
    }
  }
}
table_2can$pair <- paste(table_2can$KINASE,table_2can$SUBSTRATE,table_2can$SUB_MOD_RSD,sep = ":")

# write out tables --------------------------------------------------------

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_BRCA_subtype_trans_edited.txt", sep="")
write.table(table_2can, file=tn, quote=F, sep = '\t', row.names = FALSE)
