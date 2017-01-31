# Yige Wu @ WashU 2016 Nov
# look at correlations of kinase and downstream substrates phosphorylation status

# choose kinase/phosphotase, significance level, outlier threshold and least sample number-------------------------
sig <- 0.05 # significance level
out_thres <- 1.5 #threshold for outlier
least_samples <- 5# least number of samples with complete data for each model
protein <- "kinase"

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

# input regardless to cancer type, choose between kinase or phosphotase-------------------------------------------------------------------
if ( protein == "kinase" ) {
  ### read in the kinase/substrate table/ phosphorylation data ###
  k_s_table = read.delim(paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
}

if ( protein == "phosphotase" ) {
  ### read in the phosphotase/substrate table/ phosphorylation data ### 
  k_s_table <- read.csv(paste(baseD,"pan3can_shared_data/Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.csv",sep = ""))
  colnames(k_s_table) <- c("Phosphatase_UniProtAC_human","GENE","Substrate_UniProtAC_ref","SUB_GENE","Substrate_Type","DephosphoSite","BioassayType", "PubMed_ID_rev")
}

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

unique_kinase <- unique(k_s_table$GENE)

# looping cancer,run following sections all at once before plot module -----------------------------------------------------------
for (cancer in c("BRCA","OV")) {
# input according to cancer type-------------------------------------------------------------------
  if (cancer == "BRCA") {
    # BRCA
    BRCA_pro_f = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized_max10NA.txt",sep="")
    pro_data <- read.delim(BRCA_pro_f)
    BRCA_pho_f = paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv",sep="")
    pho_data = read.delim(BRCA_pho_f)
    ## read in grouped phosphorylation data!
    BRCA_pho_g = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_by_PRO_formatted_normalized_max10NA.txt",sep="")
    pho_gdata = read.delim(BRCA_pho_g)
    colx <- 78 # the column of protein name
  }
  
  if ( cancer == "OV" ) {
    #OV
    cancer = "OV"
    OV_pho_f = paste(baseD,"pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA.tsv",sep="")
    pho_data = read.delim(OV_pho_f)
    OV_pho_g = paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_by_PRO_formatted_normalized_max10NA.txt",sep="")
    pho_gdata = read.delim(OV_pho_g)
    OV_pro_f = paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized_max10NA.txt",sep="")
    pro_data <- read.delim(OV_pro_f)
    colnames(pro_data) <- str_split_fixed(colnames(pro_data),"_",2)[,1]
    pro_data <- pro_data[,colnames(pho_data)]
    colx <- 1 # the column of protein name
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
  
# initiate ----------------------------------------------------------------
  #initiating the lists to store info for each model
  list1 <- list()
  list2 <- list()
  list3 <- list()
  list4 <- list()
  
# looping for model1,2,3 -----------------------------------------------------------------
  for (kinase in unique_kinase){
    #for (kinase in "ERBB2"){#test
    # find protein expression level for the kinase
    pro_kinase <- pro_data[pro_data$X == kinase,-colx]
    
    if(nrow(pro_kinase) != 0){
      # find grouped phosphorylation level for kinase
      pho_kinase_g <- pho_gdata[pho_gdata$X == kinase,-colx]
      pho_kin_g_norm <- range01(unlist(pho_kinase_g),na.rm = T)
      
      # find its substrate set
      k_k_s_table = k_s_table[k_s_table$GENE == kinase,]
      k_sub <- unique(k_k_s_table$SUB_GENE)
      
      for (substrate in k_sub){# for each substrate for one kinase
        # find its phosphosites-row numbers
        s_pho_table <- which(pho_rsd_split[,1]==substrate)
        
        # go through all the phosphosites
        for (i in s_pho_table) {
          
          # find phosphorylation level
          sub_mod_rsd <- pho_rsd_split[i,3]
          pho_sub <- pho_data[i,-colx]
          
          #normalize phospho level of substrate and protein expression level of the kinase to 0 to 1
          pho_sub_norm <- range01(unlist(pho_sub),na.rm = T)
          pro_kin_norm <- range01(unlist(pro_kinase),na.rm = T)
          
          #prepare regression data for model1
          data1 <- data.frame(pho_sub_norm,pro_kin_norm)
          
          size <- nrow(data1[complete.cases(data1),])
          if( size > least_samples ){#more than 2 complete dataset
            # fit regression model1: pho_substrate ~ a*pro_kinase + k
            fit1 <- glm(pho_sub_norm ~ pro_kin_norm,data = data1, family=gaussian())
            
            # record the kinase name, kinase expression level, substrate name, SUB_MOD_RSD, phophorylation level into the list1
            list1[[length(list1)+1]] <- list(KINASE=kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,Pvalue=c(coef(summary(fit1))[2,4]),Coef_pro_kin=fit1$coefficients[2],Size=size)
            
            # find substrate expressio level and normaize
            pro_sub <- pro_data[pro_data$X == substrate,-colx]
            pro_sub_norm <- range01(unlist(pro_sub),na.rm = T)
            
            if(nrow(pro_sub) != 0){
              #prepare regression data for model2
              data2 <- data.frame(data1,pro_sub_norm)
              
              size <- nrow(data2[complete.cases(data2),])
              if(size > least_samples ){#because we need pvalue to be not NAN
                # fit regressio model2: sub-P ~ a*kinase + b*sub + k
                fit2 <- glm(pho_sub_norm ~ pro_kin_norm + pro_sub_norm, data = data2, family = gaussian())
                pvalue2 <- coef(summary(fit2))[-1,4]
                if(length(pvalue2) == 1){#this means pro_sub and pro_kinase have linear relation, most likely substrate=kinase
                  pvalue2 <- rbind(pvalue2,NA)
                }
                
                # recording like model1 + substrate expression level
                list2[[length(list2)+1]] <- list(KINASE=kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,Pvalue= pvalue2,Coef_pro_kin=fit2$coefficients[2], Coef_pro_sub=fit2$coefficients[3],Size=size)
                
                if(nrow(pho_kinase_g) != 0){
                  #prepaper regression data for model3
                  data3 <- data.frame(data2,pho_kin_g_norm)
                  size <- nrow(data3[complete.cases(data3),])
                  
                  if( size > least_samples){
                    # sub-P ~ a*pro_kinase + b*pro_sub + c*pho_kinase_grouped
                    fit3 <- glm(pho_sub_norm ~ pro_kin_norm + pro_sub_norm + pho_kin_g_norm, data = data3, family = gaussian())
                    
                    pvalue3 <- coef(summary(fit3))[-1,4]
                    if(length(pvalue3)==2){#this means pro_sub and pro_kinase have linear relation, most likely substrate=kinase
                      pvalue3 <- rbind(pvalue3[1],NA,pvalue3[2])
                    }
                    # recording like model1 + substrate expression level
                    list3[[length(list3)+1]] <- list(KINASE=kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,Pvalue= pvalue3,Coef_pro_kin=fit3$coefficients[2], Coef_pro_sub=fit3$coefficients[3],Coef_pho_kin=fit3$coefficients[4],Size=size )
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
# looping for model4: pho_sub~pho_kin_grouped -------------------------------------------------
  for (kinase in unique_kinase){
    #for (kinase in "ERBB2"){#test
    # find protein expression level for the kinase
    
    # find grouped phosphorylation level for kinase
    pho_kinase_g <- pho_gdata[pho_gdata$X == kinase,-colx]
    if (nrow(pho_kinase_g) != 0){
      pho_kin_g_norm <- range01(unlist(pho_kinase_g),na.rm = T)
      # find its substrate set
      k_k_s_table = k_s_table[k_s_table$GENE == kinase,]
      k_sub <- unique(k_k_s_table$SUB_GENE)
      
      for (substrate in k_sub){# for each substrate for one kinase
        # find its phosphosites-row numbers
        s_pho_table <- which(pho_rsd_split[,1]==substrate)
        
        # go through all the phosphosites
        for (i in s_pho_table) {
          
          # find phosphorylation level
          sub_mod_rsd <- pho_rsd_split[i,3]
          pho_sub <- pho_data[i,-colx]
          
          #normalize phospho level of substrate and protein expression level of the kinase to 0 to 1
          pho_sub_norm <- range01(unlist(pho_sub),na.rm = T)
          pro_kin_norm <- range01(unlist(pro_kinase),na.rm = T)
          
          #prepare regression data for model1
          data1 <- data.frame(pho_sub_norm,pho_kin_g_norm)
          size <- length(which(complete.cases(data1)))
          if( size > least_samples ){#more than 2 complete dataset
            # fit regression model1: pho_substrate ~ a*pro_kinase + k
            fit1 <- glm(pho_sub_norm ~ pho_kin_g_norm,data = data1, family=gaussian())
            
            # record the kinase name, kinase expression level, substrate name, SUB_MOD_RSD, phophorylation level into the list1
            list4[[length(list4)+1]] <- list(KINASE=kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,Pvalue=c(coef(summary(fit1))[2,4]),Coef_pho_kin=fit1$coefficients[2],Size=size)
          }
        }
      }
    }
  }
  
# integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
  # construct table for model1
  # columns: KINASE, SUBSTRATE, SUB_MOD_RSD, Pvalue, size, pair, self
  KINASE <- sapply(list1, "[[", "KINASE")
  SUBSTRATE <- sapply(list1, "[[", "SUBSTRATE")
  SUB_MOD_RSD  <- sapply(list1, "[[", "SUB_MOD_RSD")
  size <- sapply(list1, "[[", "Size")
  table1 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,size)
  table1$model <- "pho_sub~pro_kin"
  table1$P_pro_kin <- sapply(list1, "[[", "Pvalue")
  table1$P_pro_sub <- NA
  table1$P_pho_kin <- NA
  table1$coef_pro_kin <- sapply(list1, "[[", "Coef_pro_kin")
  table1$coef_pro_sub <- NA
  table1$coef_pho_kin <- NA
  table1$pair <- paste(table1$KINASE,":",table1$SUBSTRATE,":",table1$SUB_MOD_RSD,sep="")
  
  # construct table for model2
  KINASE <- sapply(list2, "[[", "KINASE")
  SUBSTRATE <- sapply(list2, "[[", "SUBSTRATE")
  SUB_MOD_RSD  <- sapply(list2, "[[", "SUB_MOD_RSD")
  size <- sapply(list2, "[[", "Size")
  table2 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,size)
  table2$model <- "pho_sub~pro_kin+pro_sub"
  Pvalue <- t(sapply(list2, "[[", "Pvalue"))
  table2$P_pro_kin <- Pvalue[,1]
  table2$P_pro_sub <- Pvalue[,2]
  table2$P_pho_kin <- NA
  table2$coef_pro_kin <- sapply(list2, "[[", "Coef_pro_kin")
  table2$coef_pro_sub <- sapply(list2, "[[", "Coef_pro_sub")
  table2$coef_pho_kin <- NA
  table2$pair <- paste(table2$KINASE,":",table2$SUBSTRATE,":",table2$SUB_MOD_RSD,sep="")
  
  # construct table for model3
  KINASE <- sapply(list3, "[[", "KINASE")
  SUBSTRATE <- sapply(list3, "[[", "SUBSTRATE")
  SUB_MOD_RSD  <- sapply(list3, "[[", "SUB_MOD_RSD")
  size <- sapply(list3, "[[", "Size")
  table3 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,size)
  table3$model <- "pho_sub~pro_kin+pro_sub+pho_kin"
  Pvalue <- t(sapply(list3, "[[", "Pvalue"))
  table3$P_pro_kin <- Pvalue[,1]
  table3$P_pro_sub <- Pvalue[,2]
  table3$P_pho_kin <- Pvalue[,3]
  table3$coef_pro_kin <- sapply(list3, "[[", "Coef_pro_kin")
  table3$coef_pro_sub <- sapply(list3, "[[", "Coef_pro_sub")
  table3$coef_pho_kin <- sapply(list3, "[[", "Coef_pho_kin")
  table3$pair <- paste(table3$KINASE,":",table3$SUBSTRATE,":",table3$SUB_MOD_RSD,sep="")
  
  # construct table for model4
  # columns: KINASE, SUBSTRATE, SUB_MOD_RSD, Pvalue, size, pair, self
  KINASE <- sapply(list4, "[[", "KINASE")
  SUBSTRATE <- sapply(list4, "[[", "SUBSTRATE")
  SUB_MOD_RSD  <- sapply(list4, "[[", "SUB_MOD_RSD")
  size <- sapply(list4, "[[", "Size")
  table4 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,size)
  table4$model <- "pho_sub~pho_kin"
  table4$P_pro_kin <- NA
  table4$P_pro_sub <- NA
  table4$P_pho_kin <- sapply(list4, "[[", "Pvalue")
  table4$coef_pro_kin <- NA
  table4$coef_pro_sub <- NA
  table4$coef_pho_kin <- sapply(list4, "[[", "Coef_pho_kin")
  table4$pair <- paste(table4$KINASE,":",table4$SUBSTRATE,":",table4$SUB_MOD_RSD,sep="")
  
  # combine table
  table <- rbind(table1,table2,table3,table4)
  
  # mark cancer and self-regulation
  table$Cancer <- cancer
  table$self <- as.character(table$KINASE) == as.character(table$SUBSTRATE)
  
  
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
table_2can <- rbind(table_BRCA,table_OV)
table_2can$SELF <- "trans"; table_2can$SELF[table_2can$self] <- "cis"
name = c("pro_kin","pro_sub","pho_kin")

## adjust p-values to FDR
for (mod in c("pho_sub~pro_kin", "pho_sub~pro_kin+pro_sub", "pho_sub~pro_kin+pro_sub+pho_kin","pho_sub~pho_kin")) {
  for (cancer in c("BRCA","OV")) {
    for(self in c(TRUE,FALSE)) {
      for(coln in name) {#adjust pvalues for each variable
        row <- (table_2can$self==self) & (table_2can$model==mod) & (table_2can$Cancer==cancer)
        table_2can[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_2can[row,paste("P_",coln,sep = "")],method = "fdr")
      }
    }
  }
}


# write out tables --------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",protein,"_substrate_regression.txt", sep="")
write.table(table_2can, file=tn, quote=F, sep = '\t', row.names = FALSE)

table_sig <- table_2can[table_2can$model=="pho_sub~pro_kin" & table_2can$FDR_pro_kin <= sig,]
table_sig_self <- table_sig[table_sig$self,]
table_sig_other <- table_sig[!table_sig$self,]

