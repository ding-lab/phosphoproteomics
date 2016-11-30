### validate_kinase_sub.R ### 
# Yige Wu @ WashU 2016 Nov
# look at correlations of kinase and downstream substrates phosphorylation status
# sub-P ~ a*kinase + k
# sub-P ~ a*kinase + b*sub + k
setwd("~/proteomics/pan3can_analysis/phospho_network")

### read in the kinase/substrate table/ phosphorylation data ### 
K_S_f ="~/Box Sync/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt"
k_s_table = read.delim(K_S_f)

## choose one between the following two blocks

# BRCA
cancer = "BRCA"
BRCA_pro_f = "~/Box Sync/BRCA/BRCA_PRO_formatted_normalized.txt"
pro_data <- read.delim(BRCA_pro_f)
BRCA_pho_f = "~/Box Sync/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
pho_data = read.delim(BRCA_pho_f)
## read in grouped phosphorylation data!
BRCA_pho_g = "~/Box Sync/BRCA/BRCA_PHO_by_PRO_formatted_normalized.txt"
pho_gdata = read.delim(BRCA_pho_g)

table_2can <- data.frame()

# #OV
# cancer = "OV"
# OV_pro_f = "~/Box Sync/OV/OV_merged_PRO_noOverlap_formatted_normalized.txt"
# pro_data <- read.delim(OV_pro_f)
# OV_pho_f = "~/Box Sync/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
# pho_data = read.delim(OV_pho_f)
# OV_pho_g = "~/Box Sync/OV/OV_PNNL_PHO_by_PRO_formatted_normalized.txt"
# pho_gdata = read.delim(OV_pho_g)

# ordering the columns by sample name
pro_data <- pro_data[,order(names(pro_data))]
pho_data <- pho_data[,order(names(pho_data))]
pho_gdata <- pho_gdata[,order(names(pho_gdata))]#order the grouped phospho data

### loop through each of the kinase
unique_kinase = unique(k_s_table$KINASE) 

#split the SUBSTRATE and SUB_MOD_RSD in the first column
library(stringr)
pho_rsd_split <- str_split_fixed(pho_data$X, ":", 3)
#covert the SUB_MOD_RSD from lowercase to uppercase
pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
colnames(pho_rsd_split) <- c("SUBSTRATE","unknown","SUB_MOD_RSD")

#function for normalize the y for glm regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

##initiating the output
#considering currently we have 3 models:
# sub-P ~ a*kinase + k
# sub-P ~ a*kinase + b*sub + k
# sub-P ~ a*pro_kinase + b*pro_sub + c*pho_kinase_grouped
#we initiate 3 lists: 
#first one containing kinase name, kinase expression level, substrate name, SUB_MOD_RSD, phophorylation level
#second one additionaly contains substrate expression level
#third one additionally contains grouped phsophorylation level for the kinase
#because sometimes we cannot find the expressio level for the substrate

#initiating the 3 lists
list1 <- list()
list2 <- list()
list3 <- list()
count <- 0

#for (kinase in unique_kinase){
for (kinase in "ERBB2"){#test
  # find protein expression level for the kinase
  pro_kinase <- pro_data[pro_data$X == kinase,-78]
  
  if(nrow(pro_kinase) != 0){
    # find grouped phosphorylation level for kinase
    pho_kinase_g <- pho_gdata[pho_gdata$X == kinase,-78]
    
    # find its substrate set
    k_k_s_table = k_s_table[k_s_table$KINASE == kinase,]
    
    #one line at a time, because there're the same substrate with different SUB_MOD_RSD for one kinase
    for (i in 1:nrow(k_k_s_table)){
      count <- count + 1
      
      substrate <- toString(k_k_s_table$SUB_GENE[i])
      sub_mod_rsd <- toString(k_k_s_table$SUB_MOD_RSD[i])
      # find phospher level
      pho_sub <- pho_data[pho_rsd_split[,1] == substrate & pho_rsd_split[,3] == sub_mod_rsd,-78]
      
      if(nrow(pho_sub) != 0){
        #normalize phospho level to 0 to 1
        pho_sub_norm <- range01(unlist(pho_sub),na.rm = T)
        
        #prepare regression data for model1
        data1 <- data.frame(pho_sub_norm,t(pro_kinase))
        
        if( nrow(data1[complete.cases(data1),]) > 2 ){#more than 2 complete dataset
          colnames(data1) <- c("pho_sub_norm","pro_kinase")
          # fit regression model1: sub-P ~ a*kinase + k
          fit1 <- glm(pho_sub_norm ~ pro_kinase,data = data1, family=gaussian())

          # record the kinase name, kinase expression level, substrate name, SUB_MOD_RSD, phophorylation level into the list1
          list1[[length(list1)+1]] <- list(KINASE=kinase,PRO_KINASE=pro_kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,PHO_SUB=pho_sub,FIT=fit1,Pvalue=c(coef(summary(fit1))[2,4]) )

          # find substrate expressio level
          pro_sub <- pro_data[pro_data$X == substrate,-78]
          if(nrow(pro_sub) != 0){
            #prepare regression data for model2
            data2 <- data.frame(data1,t(pro_sub))
            
            if(nrow(data2[complete.cases(data2),]) > 3 ){#because we need pvalue to be not NAN
              colnames(data2) <- c("pho_sub_norm","pro_kinase","pro_sub")
              # fit regressio model2: sub-P ~ a*kinase + b*sub + k
              fit2 <- glm(pho_sub_norm ~ pro_kinase + pro_sub, data = data2, family = gaussian())

              # recording like model1 + substrate expression level
              list2[[length(list2)+1]] <- list(KINASE=kinase,PRO_KINASE=pro_kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,PHO_SUB=pho_sub,PRO_SUB=pro_sub,FIT=fit2,Pvalue=c(coef(summary(fit2))[2,4]) )

              if(nrow(pho_kinase_g) != 0){
                #prepaper regression data for model3
                data3 <- data.frame(data2,t(pho_kinase_g))
                
                if(nrow(data3[complete.cases(data3),]) > 4){
                  colnames(data3) <- c("pho_sub_norm","pro_kinase","pro_sub","pho_kinase_g")
                  # sub-P ~ a*pro_kinase + b*pro_sub + c*pho_kinase_grouped
                  fit3 <- glm(pho_sub_norm ~ pro_kinase + pro_sub + pho_kinase_g, data = data3, family = gaussian())
                  
                  # recording like model1 + substrate expression level
                  list3[[length(list3)+1]] <- list(KINASE=kinase,PRO_KINASE=pro_kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,PHO_SUB=pho_sub,PRO_SUB=pro_sub,PHO_KINASE_G=pho_kinase_g,FIT=fit3,Pvalue=c(coef(summary(fit3))[2,4]) )
                }
              }
            }
          }
        }
      }
    }
  }
}


## count how many pairs of kinase/substrate is significant in regression model1/2
Pvalue1 <- sapply(list1, "[[", "Pvalue")
num_sig1 <- length(Pvalue1[Pvalue1 <= 0.05])
Pvalue2 <- sapply(list2, "[[", "Pvalue")
num_sig2 <- length(Pvalue2[Pvalue2 <= 0.05])

save.image(paste("~/proteomics/pan3can_analysis/phospho_network/",cancer,"_2models.RData",sep=""))
