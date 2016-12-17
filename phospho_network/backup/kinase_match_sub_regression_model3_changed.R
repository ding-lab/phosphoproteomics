# Yige Wu @ WashU 2016 Nov
# look at correlations of kinase and downstream substrates phosphorylation status

## choose one between the following two cancers to process
cancer = "BRCA"
# cancer = "OV"
sig <- 0.05


# input -------------------------------------------------------------------
if (cancer == "BRCA") {
  # BRCA
  cancer = "BRCA"
  BRCA_pro_f = "~/Box Sync/BRCA/BRCA_PRO_formatted_normalized.txt"
  pro_data <- read.delim(BRCA_pro_f)
  BRCA_pho_f = "~/Box Sync/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
  pho_data = read.delim(BRCA_pho_f)
  ## read in grouped phosphorylation data!
  BRCA_pho_g = "~/Box Sync/BRCA/BRCA_PHO_by_PRO_formatted_normalized.txt"
  pho_gdata = read.delim(BRCA_pho_g)
  colx <- 78 # the column of protein name
} else {
  #OV
  cancer = "OV"
  OV_pho_f = "~/Box Sync/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
  pho_data = read.delim(OV_pho_f)
  OV_pho_g = "~/Box Sync/OV/OV_PNNL_PHO_by_PRO_formatted_normalized.txt"
  pho_gdata = read.delim(OV_pho_g)
  OV_pro_f = "~/Box Sync/OV/OV_merged_PRO_noOverlap_formatted_normalized.txt"
  pro_data_merged <- read.delim(OV_pro_f)
  pro_data <- pro_data_merged[,colnames(pho_data)]
  colx <- 1 # the column of protein name
}

# ordering the columns by sample name
pro_data <- pro_data[,order(names(pro_data))]
pho_data <- pho_data[,order(names(pho_data))]
pho_gdata <- pho_gdata[,order(names(pho_gdata))]#order the grouped phospho data

### read in the kinase/substrate table/ phosphorylation data ### 
K_S_f ="~/Box Sync/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt"
k_s_table = read.delim(K_S_f)

#split the SUBSTRATE and SUB_MOD_RSD in the first column
library(stringr)
pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))

#covert the SUB_MOD_RSD from lowercase to uppercase
pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

# initiate ----------------------------------------------------------------

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
list4 <- list()

count <- 0
# looping substrate -----------------------------------------------------------------

#for (kinase in unique_kinase){
for (kinase in "ERBB2"){#test
  # find protein expression level for the kinase
  pro_kinase <- pro_data[pro_data$X == kinase,-colx]
  
  # find grouped phosphorylation level for kinase
  pho_kinase_g <- pho_gdata[pho_gdata$X == kinase,-colx]
  pho_kin_g_norm <- range01(unlist(pho_kinase_g),na.rm = T)
  
  if(nrow(pro_kinase) != 0){
    # find its substrate set
    k_k_s_table = k_s_table[k_s_table$KINASE == kinase,]
    k_sub <- unique(k_k_s_table$SUB_GENE)
    
    for (substrate in k_sub){# for each substrate for one kinase
      # find its phosphosites-row numbers
      s_pho_table <- which(pho_rsd_split[,1]==substrate)
      
      # go through all the phosphosites
      for (i in s_pho_table) {
        count <- count + 1#calculate how many pairs we examined
        
        # find phosphorylation level
        sub_mod_rsd <- pho_rsd_split[i,3]
        pho_sub <- pho_data[i,-colx]

        #normalize phospho level of substrate and protein expression level of the kinase to 0 to 1
        pho_sub_norm <- range01(unlist(pho_sub),na.rm = T)
        pro_kin_norm <- range01(unlist(pro_kinase),na.rm = T)
        
        if(nrow(pho_kinase_g) != 0){
          #prepaper regression data for model3
          data3 <- data.frame(pho_sub_norm,pho_kin_g_norm)
          size <- nrow(data3[complete.cases(data3),])
          
          if( size > 2){
            # sub-P ~ c*pho_kinase_grouped
            fit3 <- glm(pho_sub_norm ~ pho_kin_g_norm, data = data3, family = gaussian())
            
            pvalue3 <- coef(summary(fit3))[-1,4]
            # recording like model1 + substrate expression level
            list3[[length(list3)+1]] <- list(KINASE=kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,DATA=data3,Pvalue= pvalue3,Coef_pho_kin_g=pvalue3,Size=size )
          }
        }

        #prepare regression data for model1
        data1 <- data.frame(pho_sub_norm,pro_kin_norm)
        
        size <- nrow(data1[complete.cases(data1),])
        if( size > 2 ){#more than 2 complete dataset
          # fit regression model1: pho_substrate ~ a*pro_kinase + k
          fit1 <- glm(pho_sub_norm ~ pro_kin_norm,data = data1, family=gaussian())

          # record the kinase name, kinase expression level, substrate name, SUB_MOD_RSD, phophorylation level into the list1
          list1[[length(list1)+1]] <- list(KINASE=kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,DATA=data1,Pvalue=c(coef(summary(fit1))[2,4]),Coef_pro_kin=fit1$coefficients[2],Size=size)

          # find substrate expressio level and normaize
          pro_sub <- pro_data[pro_data$X == substrate,-colx]
          pro_sub_norm <- range01(unlist(pro_sub),na.rm = T)
          
          if(nrow(pro_sub) != 0){
            #prepare regression data for model2
            data2 <- data.frame(data1,pro_sub_norm)
            
            size <- nrow(data2[complete.cases(data2),])
            if(size > 3 ){#because we need pvalue to be not NAN
              # fit regressio model2: sub-P ~ a*kinase + b*sub + k
              fit2 <- glm(pho_sub_norm ~ pro_kin_norm + pro_sub_norm, data = data2, family = gaussian())
              pvalue2 <- coef(summary(fit2))[-1,4]
              if(length(pvalue2) == 1){#this means pro_sub and pro_kinase have linear relation, most likely substrate=kinase
                pvalue2 <- rbind(pvalue2,NA)
              }

              # recording like model1 + substrate expression level
              list2[[length(list2)+1]] <- list(KINASE=kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,DATA=data2,Pvalue= pvalue2,Coef_pro_kin=fit2$coefficients[2], Coef_pro_sub=fit2$coefficients[3],Size=size)
            }
          }
        }
      }
    }
  }
}

# integrate_output --------------------------------------------------------
# print result statistics -------------------------------------------------
# construct table for bubble chart then do analysis over gain for the other cancer dataset -----------------------------------------------------
# bubble chart module -----------------------------------------------------
