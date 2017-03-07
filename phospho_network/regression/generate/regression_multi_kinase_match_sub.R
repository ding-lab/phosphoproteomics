# Yige Wu @ WashU 2017 Feb
# look at correlations of kinase and downstream substrates phosphorylation status
# pho_sub~pro_sub+pho_kin(callapsed,multiple)

# choose kinase/phosphotase, significance level, outlier threshold and least sample number-------------------------
least_samples <- 5# least number of samples with complete data for each model
protein <- "kinase"
# protein <- "phosphotase"

# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)


# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))
source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# input k_s_table according to kinase or phosphotase-------------------------------------------------------------------
if ( protein == "kinase" ) {
  ### read in the kinase/substrate table/ phosphorylation data ###
  k_s_table = read.delim(paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
}

if ( protein == "phosphotase" ) {
  ### read in the phosphotase/substrate table/ phosphorylation data ### 
  k_s_table <- read.csv(paste(baseD,"pan3can_shared_data/Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.csv",sep = ""))
  colnames(k_s_table) <- c("Phosphatase_UniProtAC_human","GENE","Substrate_UniProtAC_ref","SUB_GENE","Substrate_Type","DephosphoSite","BioassayType", "PubMed_ID_rev")
}

substrate_trans <- as.vector(unique(k_s_table$SUB_GENE[as.vector(k_s_table$GENE)!=as.vector(k_s_table$SUB_GENE)]))
kinase_trans <- as.vector(unique(k_s_table$GENE[as.vector(k_s_table$GENE)!=as.vector(k_s_table$SUB_GENE)]))
table_single <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))

# looping cancer -----------------------------------------------------------
for (cancer in c("BRCA","OV")) { 
# input according to cancer type-------------------------------------------------------------------
  if (cancer == "BRCA") {
    # BRCA
    BRCA_pro_f = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized_max10NA.txt",sep="")
    pro_data <- read.delim(BRCA_pro_f)
    BRCA_pho_f = paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv",sep="")
    pho_data = read.delim(BRCA_pho_f)
    ## read in grouped phosphorylation data!
    BRCA_pho_g = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized_max10NA.txt",sep="")
    pho_gdata = read.delim(BRCA_pho_g)
  }
  
  if ( cancer == "OV" ) {
    OV_pho_f = paste(baseD,"pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA.tsv",sep="")
    pho_data = read.delim(OV_pho_f)
    OV_pho_g = paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized_max10NA.txt",sep="")
    pho_gdata = read.delim(OV_pho_g)
    OV_pro_f = paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized_max10NA.txt",sep="")
    pro_data <- read.delim(OV_pro_f)
    colnames(pro_data) <- str_split_fixed(colnames(pro_data),"_",2)[,1]
    pro_data <- pro_data[,colnames(pho_data)]
  }

  # ordering the columns by sample name
  pro_data <- pro_data[,order(names(pro_data))]
  pho_data <- pho_data[,order(names(pho_data))]
  pho_gdata <- pho_gdata[,order(names(pho_gdata))]#order the grouped phospho data
  
  #split the SUBSTRATE and SUB_MOD_RSD in the first column
  pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))
  
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
  colx <- which(colnames(pro_data)=="X")
# initiate ----------------------------------------------------------------
# calculate the length of trans table
  substrate_reg <- intersect(intersect(substrate_trans, pro_data$X),pho_rsd_split$SUBSTRATE) 
  # use pro data and phosphosite data to narrow down the substrates to be examined(2119 to 1364)
  ntrans <- 0
  for (gene in kinase_trans) {
    subs <- k_s_table$SUB_GENE[k_s_table$GENE==gene & k_s_table$SUB_GENE!=gene]
    for ( sub in unique(subs)) {
      ntrans <- ntrans + length(which(pho_rsd_split$SUBSTRATE==sub))
    }
  }
  
  phosphosites <- unique(table_single[table_single$Cancer==cancer & !table_single$self,c("SUBSTRATE","SUB_MOD_RSD","transcript")])
  substrates <- as.vector(phosphosites$SUBSTRATE)
  rsds <- as.vector(phosphosites$SUB_MOD_RSD)
  for (i in 1:nrow(phosphosites)) {
    phosphosites$pho_size[i] <- length(which(!is.na(pho_data[pho_rsd_split$SUBSTRATE==substrates[i] & pho_rsd_split$SUB_MOD_RSD==rsds[i],-colx])))
    phosphosites$num_k[i] <- nrow(table_single[table_single$Cancer==cancer & !table_single$self  & table_single$SUBSTRATE==substrates[i] & table_single$SUB_MOD_RSD==rsds[i],])
  }
  phosphosites_multi <- phosphosites[phosphosites$num_k > 1,]
  
  substrates <- as.vector(phosphosites_multi$SUBSTRATE)
  rsds <- as.vector(phosphosites_multi$SUB_MOD_RSD)
  
# looping over phosphosites for trans pairs -----------------------------------------------------------------
  # initiating the table for trans
  vec_char <- vector(mode = "character", length = ntrans)
  vec_num <- vector(mode = "numeric", length = ntrans) + NA
  KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
  FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
  coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
  Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
  Size <- vec_num;num_k <- vec_num;P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
  
  count <- 0
  examine1 <- vector(mode = "logical", length = nrow(phosphosites_multi))
  for (i in 1:nrow(phosphosites_multi)) {
    substrate <- substrates[i]
    rsd <- rsds[i]
    
    pho_row <- which(pho_rsd_split$SUBSTRATE==substrate & pho_rsd_split$SUB_MOD_RSD==rsd)
    tscp <- transcripts[pho_row]
    pho_sub <- t(pho_data[pho_row,-colx])
    pro_sub <- t(pro_data[pro_data$X == substrate,-colx])
    sub_k <- as.vector(unique(table_single$KINASE[table_single$Cancer==cancer & !table_single$self & table_single$SUBSTRATE==substrate & table_single$SUB_MOD_RSD==rsd]))
    num_k <- length(sub_k)
    
    rows <- c()
    for (k in sub_k) {
      temp <- which(pho_gdata$X==k)
      rows <- c(rows,temp)
    }
    pho_kins <- t(pho_gdata[rows,-colx])
    data1 <- data.frame(pro_sub,pho_kins)
    colnames(data1)[1] <- "pro_sub"
    
    data2 <- data.frame(pho_sub,data1)
    colnames(data2)[1] <- "pho_sub"
    
    size <- nrow(data2[complete.cases(data2),])
    if(size > least_samples & size >= num_k + 3){
      fit2 <- lm(pho_sub ~ ., data = data2)
      pvalues <- coef(summary(fit2))
      coefs <- fit2$coefficients
      
      examine1[i] <- TRUE
      
      for (j in 1:num_k) {
        count <- count + 1
        KINASE[count] <- sub_k[j]
        SUBSTRATE[count] <- substrate
        SUB_MOD_RSD[count] <- rsd
        transcript[count] <- transcripts[i]
        P_pro_sub[count] <- pvalues[2,4]
        coef_pro_sub[count] <- fit2$coefficients[2]
        P_pho_kin[count] <- pvalues[j+2,4]
        coef_pho_kin[count] <- coefs[j+2]
        Size[count] <- size
      }
    }
  }
  table_trans <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,
                            FDR_pro_kin,FDR_pro_sub,FDR_pho_kin,
                            coef_pro_kin,coef_pro_sub,coef_pho_kin,
                            Cancer,transcript,model,Size,
                            P_pro_kin,P_pro_sub,P_pho_kin)
  phosphosites_multi$examine1 <- examine1
  
  
# integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
  table_trans$model <- "pho_sub~pro_sub+pho_kins"
  tabletrans <- table_trans[!is.na(table_trans$P_pro_sub),]
  
  # combine table
  table <- table_trans

  # mark cancer and self-regulation
  table$Cancer <- cancer
  
  #choose one command
  if ( cancer == "BRCA" ) {
    table_BRCA <- table
  } 
  if ( cancer == "OV" ) {
    table_OV <- table
  }
}

# combine table from BRCA and OV and adjust -------------------------------------------------------
## combine table from BRCA and OV
temp <- rbind(table_BRCA,table_OV)
table_2can <- unique(temp) ## because there're duplicate phophorylation levels for a residual
table_2can$self <- as.character(table_2can$KINASE) == as.character(table_2can$SUBSTRATE)
table_2can$SELF <- "trans"; table_2can$SELF[table_2can$self] <- "cis"
name = c("pro_kin","pro_sub","pho_kin")

## adjust p-values to FDR
for (cancer in c("BRCA","OV")) {
  for(self in c(TRUE,FALSE)) {
    for(coln in name) {#adjust pvalues for each variable
      row <- (table_2can$self==self) & (table_2can$Cancer==cancer)
      table_2can[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_2can[row,paste("P_",coln,sep = "")],method = "fdr")
    }
  }
}

table_2can$pair <- paste(table_2can$KINASE,table_2can$SUBSTRATE,table_2can$SUB_MOD_RSD,sep = ":")
# write out tables --------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/multi_",protein,"_substrate_regression.txt", sep="")
write.table(table_2can, file=tn, quote=F, sep = '\t', row.names = FALSE)
