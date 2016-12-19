### validate_multiple_kinase_sub.R ### 
# Yige Wu @ WashU 2016 Nov
# look at correlations of kinase and downstream substrates phosphorylation status
#setwd("~/proteomics/pan3can_analysis/phospho_network")
setwd("~/Desktop/Ding Lab/figures")

## choose one between the following two cancers

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

# #OV
# cancer = "OV"
# OV_pho_f = "~/Box Sync/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
# pho_data = read.delim(OV_pho_f)
# OV_pho_g = "~/Box Sync/OV/OV_PNNL_PHO_by_PRO_formatted_normalized.txt"
# pho_gdata = read.delim(OV_pho_g)
# OV_pro_f = "~/Box Sync/OV/OV_merged_PRO_noOverlap_formatted_normalized.txt"
# pro_data_merged <- read.delim(OV_pro_f)
# pro_data <- pro_data_merged[,colnames(pho_data)]
# colx <- 1 # the column of protein name

# ordering the columns by sample name
pro_data <- pro_data[,order(names(pro_data))]
pho_data <- pho_data[,order(names(pho_data))]
pho_gdata <- pho_gdata[,order(names(pho_gdata))]#order the grouped phospho data

### read in the kinase/substrate table/ phosphorylation data ### 
K_S_f ="~/Box Sync/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt"
k_s_table = read.delim(K_S_f)

# fill in the blanks in the SUB_GENE column
# blank <- which(is.na(k_s_table$SUB_GENE_ID))
# substrate <- str_split_fixed(k_s_table$SUBSTRATE[blank]," ",2)[,1]
# for (r in blank){
#   k_s_table$SUB_GENE[r] <- substrate[which(blank==r)]
# }

#split the SUBSTRATE and SUB_MOD_RSD in the first column
library(stringr)
pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))
#covert the SUB_MOD_RSD from lowercase to uppercase
pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")

#function for normalize the y for glm regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

### loop through each of the substrate
unique_substrate = unique(k_s_table$SUB_GENE) 

##initiating the output
# for each model, one list for recording substrate, phosphosite, number of kinases, kinases, P-values, coeffcients, data size
# model 1: pho_sub~pro_kin1+pro_kin2+...
# model 2: pho_sub~pro_kin1+pro_kin2+...+pro_sub
# model 3: pho_sub~pro_kin1+pro_kin2+...+pro_sub+pho_kin1_g+pho_kin2_g+...

#initiating the 3 lists
list1 <- list()
list2 <- list()
list3 <- list()

## initiating 3 vectors alongside k_s_table to record whether a certain kinase:substrate:phosphosite is visited/significant/not-significant in the regression model 1/2/3
# 0 for unvisited, 1 for significant, -1 for not significant
visit_mod1 <- vector(mode = "integer", length = length(k_s_table))
visit_mod2 <- visit_mod1
visit_mod3 <- visit_mod1

for (substrate in unique_substrate){
#for (substrate in "AKT1"){#test
  # find substrate expressio level and normaize
  pro_sub <- pro_data[pro_data$X == substrate,-colx]
  pro_sub_norm <- range01(unlist(pro_sub),na.rm = T)
  
  ### loop through each phosphosite of the substrate
  s_rsd_table <- unique(k_s_table$SUB_MOD_RSD[k_s_table$SUB_GENE==substrate])
  s_rsd_table <- s_rsd_table[!is.na(s_rsd_table)]

  for (sub_mod_rsd in s_rsd_table) {
    
    # find phosphorylation level for the phosphosite
    pho_sub <- pho_data[(pho_rsd_split$SUBSTRATE==substrate) & (pho_rsd_split$SUB_MOD_RSD==sub_mod_rsd),-colx]
    
    # if phosphorylation level can be found
    if (nrow(pho_sub) != 0){
      #normalize phospho level of substrate and protein expression level of the kinase to 0 to 1
      pho_sub_norm <- range01(unlist(pho_sub),na.rm = T)
      
      # initiate the dataset for model
      data1 <- data.frame(pho_sub_norm)
      
      # find kinases for substrate:phosphosite pairs
      k_s_rsd_table <- as.character(k_s_table$GENE[k_s_table$SUB_GENE==substrate & k_s_table$SUB_MOD_RSD==sub_mod_rsd])  
      k_s_rsd_table <- unique(k_s_rsd_table)
      
      # record the number of kinases with protein expression level
      k_pro_table <- c()
      
      for (kinase in k_s_rsd_table) {
        # collect protein expression levels for the kinases
        pro_kinase <- pro_data[pro_data$X == kinase,-colx]
        
        # if protein expression level for kinases can be found
        if (nrow(pro_kinase) != 0){
          k_pro_table[length(k_pro_table)+1] <- kinase
          
          # normalize protein expression level
          pro_kin_norm <- range01(unlist(pro_kinase),na.rm = T)
          
          # add to dataset for model1: pho_sub~pro_kin1+pro_kin2+...
          data1 <- cbind(pro_kin_norm,data1)
          
          # modify the colname
          colnames(data1)[1] <- paste("pro_",kinase,sep="")
        }
      }
      
      size_kin <- length(k_pro_table)
      if(size_kin > 0){
        size_dat <- nrow(data1[complete.cases(data1),])   # record the number of samples having valid data for model1
        
        if( size_dat > size_kin+1 ){# if the above samples size exceed the number of variables in the model1
          # fit regression model1
          fit1 <- glm(pho_sub_norm ~ .,data = data1, family=gaussian())
          
          pvalue <- data.frame(t(coef(summary(fit1))[-1,4]))
          
          # record the regression results
          list1[[length(list1)+1]] <- list(SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,SIZE_kin=size_kin,KINASE=k_pro_table,SIZE_dat=size_dat,
                                           P_pro_kin=pvalue,COEF_pro_kin=data.frame(t(fit1$coefficients[-1])),SIG_pro_kin=length(which(pvalue <= 0.05)) )
          
          # if protein expression level for substrate can be found, proceed to to model2
          if (nrow(pro_sub) != 0) {
            data2 <- cbind(pro_sub_norm,data1)  # add pro_sub into data1 for model2
            
            size_dat <- nrow(data2[complete.cases(data2),])   # record the number of samples having valid data for model2
            if( size_dat > size_kin+2 ){# if the above samples size exceed the number of variables in the model2
              # fit regression model2
              fit2 <- glm(pho_sub_norm ~ .,data = data2, family=gaussian())
              
              pvalue <- data.frame(t(coef(summary(fit2))[-1,4]))
              if (nrow(pvalue) == size_kin ) {# which means there's kinase == substrate
                pvalue[length(pvalue)+1] <- NA
                colnames(pvalue)[length(pvalue)] <- paste("pro_",substrate,sep="")
              }
              
              Coefs <- data.frame(t(fit2$coefficients[-1]))
              col_pro_kin <- paste("pro_",k_pro_table,sep = "")
              # record the regression results
              list2[[length(list2)+1]] <- list(SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,SIZE_kin=size_kin,KINASE=k_pro_table,SIZE_dat=size_dat,
                                               P_pro_kin=pvalue[col_pro_kin],COEF_pro_kin=Coefs[col_pro_kin],SIG_pro_kin=length(which(pvalue[col_pro_kin] <= 0.05)),
                                               P_pro_sub=pvalue[1,"pro_sub_norm"],COEF_pro_sub=Coefs[1,"pro_sub_norm"],SIG_pro_sub=length(which(pvalue[1,"pro_sub_norm"] <= 0.05)))
              
              # record the number of kinases with protein expression level and grouped phosphorylation level
              k_pho_table <- c()
              data3 <- data2[,c("pro_sub_norm","pho_sub_norm")]
              
              for (kinase in k_pro_table) {
                # collect grouped phosphorylation levels for the kinases
                pho_g_kinase <- pho_gdata[pro_data$X == kinase,-colx]
                
                # if phosphorylation level for kinases can be found
                if (nrow(pho_g_kinase) != 0){
                  k_pho_table[length(k_pho_table)+1] <- kinase
                  
                  # normalize phosphorylaiton level
                  pho_g_kin_norm <- range01(unlist(pho_g_kinase),na.rm = T)
                  
                  # add to dataset for model3
                  data3 <- cbind(pho_g_kin_norm,data3)
                  
                  # modify the colname
                  colnames(data3)[1] <- paste("pho_g_",kinase,sep="")
                }
              }
              
              size_kin <- length(k_pho_table)
              if(size_kin > 0){
                col_pho_kin <- paste("pro_",k_pho_table,sep="")
                data3 <- cbind(data2[,c(col_pho_kin)],data3,deparse.level = 1)
                if(size_kin==1){
                  colnames(data3)[1] <- paste("pro_",k_pho_table,sep = "")
                }
                
                size_dat <- nrow(data3[complete.cases(data3),])   # record the number of samples having valid data for model3
                
                if( size_dat > size_kin+3 ){# if the above samples size exceed the number of variables in the model2
                  # fit regression model3
                  fit3 <- glm(pho_sub_norm ~ .,data = data3, family=gaussian())
                  
                  pvalue <- data.frame(t(coef(summary(fit3))[-1,4]))
                  if (nrow(pvalue) == 2*size_kin ) {# which means there's kinase == substrate
                    pvalue[length(pvalue)+1] <- NA
                    colnames(pvalue)[length(pvalue)] <- paste("pro_",substrate,sep="")
                  }
                  
                  Coefs <- data.frame(t(fit3$coefficients[-1]))
                  
                  col_pho_kin <- paste("pho_g_",k_pho_table,sep="")
                  # record the regression results
                  list3[[length(list3)+1]] <- list(SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,SIZE_kin=size_kin,KINASE=k_pho_table,SIZE_dat=size_dat,
                                                   P_pro_kin=pvalue[col_pro_kin],COEF_pro_kin=Coefs[col_pro_kin],SIG_pro_kin=length(which(pvalue[col_pro_kin] <= 0.05)),
                                                   P_pro_sub=pvalue[1,"pro_sub_norm"],COEF_pro_sub=Coefs[1,"pro_sub_norm"],SIG_pro_sub=length(which(pvalue[1,"pro_sub_norm"] <= 0.05)),
                                                   P_pho_kin=pvalue[col_pho_kin],COEF_pho_kin=Coefs[col_pho_kin],SIG_pho_kin=length(which(pvalue[col_pho_kin] <= 0.05)))
                }
              }
            }
          }
        }
      }
    }
  }
}





# construct table for model1
# columns: KINASE, SUBSTRATE, SUB_MOD_RSD, Pvalue, size, pair, self
SUBSTRATE <- sapply(list1, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list1, "[[", "SUB_MOD_RSD")
SIZE_dat <- sapply(list1, "[[", "SIZE_dat")
SIZE_kin <- sapply(list1, "[[", "SIZE_kin")
SIG_pro_kin <- sapply(list1, "[[", "SIG_pro_kin")
table1 <- data.frame(SUBSTRATE,SUB_MOD_RSD,SIZE_dat,SIZE_kin,SIG_pro_kin)


SUBSTRATE <- sapply(list2, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list2, "[[", "SUB_MOD_RSD")
SIZE_dat <- sapply(list2, "[[", "SIZE_dat")
SIZE_kin <- sapply(list2, "[[", "SIZE_kin")
SIG_pro_kin <- sapply(list2, "[[", "SIG_pro_kin")
SIG_pro_sub <- sapply(list2, "[[", "SIG_pro_sub")
P_pro_sub <- sapply(list2, "[[", "P_pro_sub")
table2 <- data.frame(SUBSTRATE,SUB_MOD_RSD,SIZE_dat,SIZE_kin,SIG_pro_kin,SIG_pro_sub,P_pro_sub)

SUBSTRATE <- sapply(list3, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list3, "[[", "SUB_MOD_RSD")
SIZE_dat <- sapply(list3, "[[", "SIZE_dat")
SIZE_kin <- sapply(list3, "[[", "SIZE_kin")
SIG_pro_kin <- sapply(list3, "[[", "SIG_pro_kin")
SIG_pro_sub <- sapply(list3, "[[", "SIG_pro_sub")
SIG_pho_kin <- sapply(list3, "[[", "SIG_pho_kin")
P_pro_sub <- sapply(list3, "[[", "P_pro_sub")
table3 <- data.frame(SUBSTRATE,SUB_MOD_RSD,SIZE_dat,SIZE_kin,SIG_pro_kin,SIG_pro_sub,SIG_pho_kin,P_pro_sub)


table1$model <- "pho_sub~pro_kinase"
table1$P_pro_kin <- sapply(list1, "[[", "Pvalue")
table1$P_pro_sub <- NA
table1$P_pho_kin <- NA
table1$coef_pro_kin <- sapply(list1, "[[", "Coef_pro_kin")
table1$coef_pro_sub <- NA
table1$coef_pho_kin_g <- NA
table1$pair <- paste(table1$KINASE,":",table1$SUBSTRATE,":",table1$SUB_MOD_RSD,sep="")





if(nrow(pro_kinase) != 0){
  # find grouped phosphorylation level for kinase
  pho_kinase_g <- pho_gdata[pho_gdata$X == kinase,-colx]
  pho_kin_g_norm <- range01(unlist(pho_kinase_g),na.rm = T)
  
  # find its substrate set
  k_k_s_table = k_s_table[k_s_table$KINASE == kinase,]
  k_sub <- unique(k_k_s_table$SUB_GENE)
  
  for (substrate in k_sub){# for each substrate for one kinase
    # find its phosphosites-row numbers
    s_pho_table <- which(pho_rsd_split[,1]==substrate)
    
    for (i in s_pho_table) {
      count <- count + 1#calculate how many pairs we examined
      
      
      
      
      
    }
  }
}



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
    
    if(nrow(pho_kinase_g) != 0){
      #prepaper regression data for model3
      data3 <- data.frame(data2,pho_kin_g_norm)
      size <- nrow(data3[complete.cases(data3),])
      
      if( size > 4){
        # sub-P ~ a*pro_kinase + b*pro_sub + c*pho_kinase_grouped
        fit3 <- glm(pho_sub_norm ~ pro_kin_norm + pro_sub_norm + pho_kin_g_norm, data = data3, family = gaussian())
        
        pvalue3 <- coef(summary(fit3))[-1,4]
        if(length(pvalue3)==2){#this means pro_sub and pro_kinase have linear relation, most likely substrate=kinase
          pvalue3 <- rbind(pvalue3[1],NA,pvalue3[2])
        }
        # recording like model1 + substrate expression level
        list3[[length(list3)+1]] <- list(KINASE=kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,DATA=data3,Pvalue= pvalue3,Coef_pro_kin=fit3$coefficients[2], Coef_pro_sub=fit3$coefficients[3],Coef_pho_kin_g=fit3$coefficients[4],Size=size )
      }
    }
  }
}
