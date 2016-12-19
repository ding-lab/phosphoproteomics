### validate_multiple_kinase_sub.R ### 
# Yige Wu @ WashU 2016 Nov
# look at correlations of kinase and downstream substrates phosphorylation status

table_2can <- c()
k_s_2can <- c()
## choose one between the following two cancers to process
cancer = "BRCA"
# cancer = "OV"

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

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

### loop through each of the substrate
unique_substrate = unique(k_s_table$SUB_GENE) 

unique_k_s_rsd <- unique(k_s_table[,c("GENE","SUB_GENE","SUB_MOD_RSD")])

# initiate ----------------------------------------------------------------
##initiating the output
# for each model, one list for recording substrate, phosphosite, number of kinases, kinases, P-values, coeffcients, data size
# model 1: pho_sub~pro_kin1+pro_kin2+...
# model 2: pho_sub~pro_kin1+pro_kin2+...+pro_sub
# model 3: pho_sub~pro_kin1+pro_kin2+...+pro_sub+pho_kin1_g+pho_kin2_g+...

#initiating the 3 lists
list1 <- list()
list2 <- list()
list3 <- list()

## initiating 3 vectors for pvalues, coefficients and FDRs for each kinase:substrate:phosphosite is examined in model1/2/3
k_s_fit <- data.frame(unique_k_s_rsd)   #unique kinase:substrate:sub_mod_rsd
x <- vector(mode = "numeric", length = nrow(unique_k_s_rsd))+NaN
temp <- data.frame(matrix(rep(x,21), ncol=21, byrow=T))
name <- c("pro_kin1","pro_kin2","pro_sub2","pro_kin3","pro_sub3","pho_g_kin3")
colnames(temp) <- c(paste("FDR_",name,sep = ""),paste("P_",name,sep = ""),paste("coef_",name,sep = ""),paste("size_dat",1:3,sep = ""))
k_s_fit <- cbind(unique_k_s_rsd,temp)

# looping substrate -----------------------------------------------------------------
for (substrate in unique_substrate){
#for (substrate in c("LMNA")){#test
  # find substrate expressio level and normaize
  pro_sub <- pro_data[pro_data$X == substrate,-colx]
  pro_sub_norm <- range01(unlist(pro_sub),na.rm = T)

  s_rsd_table <- unique(k_s_fit$SUB_MOD_RSD[k_s_fit$SUB_GENE==substrate])
  s_rsd_table <- s_rsd_table[!is.na(s_rsd_table)]
  
  # looping phosphosite -----------------------------------------------------
  for (sub_mod_rsd in s_rsd_table) {
  #for (sub_mod_rsd in "T19") {#test
      
    # find phosphorylation level for the phosphosite
    pho_sub <- pho_data[(pho_rsd_split$SUBSTRATE==substrate) & (pho_rsd_split$SUB_MOD_RSD==sub_mod_rsd),-colx]
    
    # if phosphorylation level can be found
    if (nrow(pho_sub) != 0){
      #normalize phospho level of substrate and protein expression level of the kinase to 0 to 1
      pho_sub_norm <- range01(unlist(pho_sub[1,]),na.rm = T)
      
      # initiate the dataset for model
      data1 <- data.frame(pho_sub_norm)
      
      # find kinases for substrate:phosphosite pairs
      s_rsd_row <- which(k_s_fit$SUB_GENE==substrate & k_s_fit$SUB_MOD_RSD==sub_mod_rsd) # row number of this phosphosite coupled with this substrate
      k_s_rsd_table <- as.character(k_s_fit$GENE[s_rsd_row])  
      k_s_rsd_table <- unique(k_s_rsd_table)
      
      k_pro_table <- c()      # record the kinases with protein expression level
      size_kin <- 0 # record the number of kinases with protein expression level
      
      k_pro_row <- c()   # record the row number for above kinases coupled with this phosphosite and substrate
      
      for (i in s_rsd_row) {
        kinase <- as.character(k_s_fit$GENE[i])
        
        # collect protein expression levels for the kinases
        pro_kinase <- pro_data[pro_data$X == kinase,-colx]
        
        # if protein expression level for kinases can be found
        if (nrow(pro_kinase) != 0){
          size_kin <- size_kin+1
          
          # record kinase name and row number in k_s_fit
          k_pro_table[size_kin] <- kinase
          k_pro_row[size_kin] <- i
          
          # normalize protein expression level
          pro_kin_norm <- range01(unlist(pro_kinase),na.rm = T)
          
          # add to dataset for model1: pho_sub~pro_kin1+pro_kin2+...
          data1 <- cbind(data1,pro_kin_norm)
          
          # modify the colname
          colnames(data1)[size_kin+1] <- paste("pro_",kinase,sep="")
        }
      }
      
      if(size_kin > 0){
        size_dat <- nrow(data1[complete.cases(data1),])   # record the number of samples having valid data for model1
        # if the above samples size exceed the number of variables in the model1, proceed to model1
        if( size_dat > size_kin+1 ){
          fit1 <- glm(pho_sub_norm ~ .,data = data1, family=gaussian())
          
          pvalues <- data.frame(t(signif(coef(summary(fit1))[-1,4],digits = 3)))
          coefs <- data.frame(t(signif(fit1$coefficients[-1],digits = 3)))
          
          # record the regression results for each substrate:phosohpsite pair for model1
          list1[[length(list1)+1]] <- list(SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,
                                           SIZE_kin=size_kin,KINASE=k_pro_table,SIZE_dat=size_dat,
                                           P_pro_kin=pvalues,COEF_pro_kin=coefs)
          
          # record the regression results for each kinase:substrate:phosohpsite pair for model1
          k_s_fit$P_pro_kin1[k_pro_row] <- pvalues
          k_s_fit$coef_pro_kin1[k_pro_row] <- coefs
          k_s_fit$size_dat1[k_pro_row] <- size_dat
          
          # if protein expression level for substrate can be found, proceed to to model2
          if (nrow(pro_sub) != 0) {
            data2 <- cbind(data1,pro_sub_norm)  # add pro_sub into data1 for model2
            
            size_dat <- nrow(data2[complete.cases(data2),])   # record the number of samples having valid data for model2
            # if the above samples size exceed the number of variables in the model2, proceed to model2
            if( size_dat > size_kin+2 ){
              fit2 <- glm(pho_sub_norm ~ .,data = data2, family=gaussian())
              
              pvalues <- data.frame(t(signif(coef(summary(fit2))[-1,4],digits = 3)))
              # make sure if kinase==substrate, the pvalue will still hold a place for pro_sub
              if (ncol(pvalues) == size_kin ) {# which means there's kinase == substrate
                pvalues[length(pvalues)+1] <- NA
                colnames(pvalues)[length(pvalues)] <- "pro_sub_norm"
                # make sure if kinase==substrate and there's only one kinase, the colnae for pvalue will be right
                if (size_kin==1) {
                  colnames(pvalues)[1] <- paste("pro_",kinase,sep = "")
                }
              }
              coefs <- data.frame(t(signif(fit2$coefficients[-1],digits = 3)))
              col_pro_kin <- paste("pro_",k_pro_table,sep = "") # colnames for pro_kin in the pvalues and coefficients
              
              # record the regression results for each substrate:phosohpsite pair for model2
              list2[[length(list2)+1]] <- list(SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,SIZE_kin=size_kin,KINASE=k_pro_table,SIZE_dat=size_dat,
                                               P_pro_kin=pvalues[col_pro_kin],COEF_pro_kin=coefs[col_pro_kin],
                                               P_pro_sub=pvalues[1,"pro_sub_norm"],COEF_pro_sub=coefs[1,"pro_sub_norm"])
              
              # record the regression results for each kinase:substrate:phosohpsite pair for model2
              k_s_fit$P_pro_kin2[k_pro_row] <- pvalues[col_pro_kin]
              k_s_fit$coef_pro_kin2[k_pro_row] <- coefs[col_pro_kin]
              k_s_fit$P_pro_sub2[k_pro_row] <- pvalues[1,"pro_sub_norm"]
              k_s_fit$coef_pro_sub2[k_pro_row] <- coefs[1,"pro_sub_norm"]
              k_s_fit$size_dat2[k_pro_row] <- size_dat
              
              # model3
              k_pho_table <- c()              # record the kinases with protein expression level and grouped phosphorylation level
              k_pho_row <- c() # record the row number of above pairs
              size_kin <- 0 # record the number of above pairs
              
              data3 <- data2[,c("pho_sub_norm","pro_sub_norm")]
              
              for (i in k_pro_row) {
                kinase <- as.character(k_s_fit$GENE[i])
                # collect grouped phosphorylation levels for the kinases
                pho_g_kinase <- pho_gdata[pro_data$X == kinase,-colx]
                
                # if phosphorylation level for kinases can be found
                if (nrow(pho_g_kinase) != 0){
                  size_kin <- size_kin+1
                  
                  # record kinase name and row number in k_s_fit
                  k_pho_table[size_kin] <- kinase
                  k_pho_row[size_kin] <- i
                  
                  # normalize phosphorylaiton level
                  pho_g_kin_norm <- range01(unlist(pho_g_kinase),na.rm = T)
                  
                  # add to dataset for model3
                  data3 <- cbind(pho_g_kin_norm,data3)
                  
                  # modify the colname
                  colnames(data3)[1] <- paste("pho_g_",kinase,sep="")
                }
              }
              
             
              if(size_kin > 0){
                # assemble dataset for model3
                col_pho_kin <- paste("pro_",k_pho_table,sep="")
                data3 <- cbind(data2[,c(col_pho_kin)],data3,deparse.level = 1)
                # make sure colname is right
                if(size_kin==1){
                  colnames(data3)[1] <- paste("pro_",k_pho_table,sep = "")
                }
                
                size_dat <- nrow(data3[complete.cases(data3),])   # record the number of samples having valid data for model3
                
                if( size_dat > 2*size_kin+2 ){# if the above samples size exceed the number of variables in the model3
                  # fit regression model3
                  fit3 <- glm(pho_sub_norm ~ .,data = data3, family=gaussian())
                  
                  pvalues <- data.frame(t(signif(coef(summary(fit3))[-1,4], digits = 3)))
                  # make sure if kinase==substrate, the pvalue will still hold a place for pro_sub
                  if (ncol(pvalues) == 2*size_kin ) {# which means there's kinase == substrate
                    pvalues[length(pvalues)+1] <- NA
                    colnames(pvalues)[length(pvalues)] <- "pro_sub_norm"
                  }
                  coefs <- data.frame(t(signif(fit3$coefficients[-1], digits = 3)))
                  col_pro_kin <- paste("pro_",k_pho_table,sep="")
                  col_pho_kin <- paste("pho_g_",k_pho_table,sep="")
                  
                  # record the regression results for each substrate:phosohpsite pair for model3
                  list3[[length(list3)+1]] <- list(SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,SIZE_kin=size_kin,KINASE=k_pho_table,SIZE_dat=size_dat,
                                                   P_pro_kin=pvalues[col_pro_kin],COEF_pro_kin=coefs[col_pro_kin],
                                                   P_pro_sub=pvalues[1,"pro_sub_norm"],COEF_pro_sub=coefs[1,"pro_sub_norm"],
                                                   P_pho_kin=pvalues[col_pho_kin],COEF_pho_kin=coefs[col_pho_kin])
                  # record the regression results for each kinase:substrate:phosohpsite pair for model3
                  k_s_fit$P_pro_kin3[k_pho_row] <- pvalues[col_pro_kin]
                  k_s_fit$coef_pro_kin3[k_pho_row] <- coefs[col_pro_kin]
                  k_s_fit$P_pro_sub3[k_pho_row] <- pvalues[1,"pro_sub_norm"]
                  k_s_fit$coef_pro_sub3[k_pho_row] <- coefs[1,"pro_sub_norm"]
                  k_s_fit$P_pho_g_kin3[k_pho_row] <- pvalues[col_pho_kin]
                  k_s_fit$coef_pho_g_kin3[k_pho_row] <- coefs[col_pho_kin]
                  k_s_fit$size_dat3[k_pho_row] <- size_dat
                  
                }
              }
            }
          }
        }
      }
    }
  }
}

# integrate_output --------------------------------------------------------
k_s_fit$self <- as.vector(k_s_fit$GENE)==as.vector(k_s_fit$SUB_GENE) # mark whether kinase is equal to substrate
# adjust pvalues to fdrs according to multiple tests and whether kinase==substrate
for(self in c(TRUE,FALSE)) {
  for(coln in name) {#adjust pvalues for each variable
    self_row <- (k_s_fit$self==self) & (!is.nan(unlist(k_s_fit[,paste("P_",coln,sep = "")]))) # only use pairs have entered in the regression model
    k_s_fit[self_row,paste("FDR_",coln,sep = "")] <-p.adjust(k_s_fit[self_row,paste("P_",coln,sep = "")],method = "fdr")
  }
}


# print result statistics -------------------------------------------------
## input row numbers of k_s_fit, output how many substrates, phosphosites and kinases, etc
KSP_count <- function(x) {
  ksp_count <- length(x)
  sub_count <- length(unique(k_s_fit$SUB_GENE[x]))
  kin_count <- length(unique(k_s_fit$GENE[x]))
  phosite_count <- nrow(unique(k_s_fit[x,c("SUB_GENE","SUB_MOD_RSD")]))
  self_count <- length(which(k_s_fit$self))
  self_ratio <- signif(self_count/ksp_count,digits = 3) 
  other_ratio <- signif((ksp_count-self_count)/ksp_count,digits = 3) 
  all_count <- data.frame(ksp_count,sub_count,kin_count,phosite_count, self_ratio, other_ratio)
  return(all_count)
}

temp <- Sys.time()
sink(file = paste("multi-kinase_match_sub_rsd_log",temp,".txt",sep = ""))

temp <- KSP_count(1:nrow(k_s_fit))
cat(paste("This k_s_table has:", "\n", temp$kin_count," kinases,", temp$sub_count," substrates, ", temp$phosite_count," phosphosites", "\n", temp$ksp_count," substrate:phosphosite:kinase pairs(denote as KSP)", "\n", "kinase == substrate: ", temp$self_ratio," KSPs ; kinase !=substrate: ", temp$other_ratio," KSPs", "\n\n", sep=" "))
cat(paste("For each unique substrate, loop around all phosphosites for the substrate;","\n", sep = " "))
cat(paste("For each phosphosite, I find all kinases known to be phosphorylate the site, then carry out the regression model.","\n\n", sep = " "))

temp <- KSP_count(which(!is.nan(unlist(k_s_fit$P_pro_kin1))))
cat(paste(temp$ksp_count," above KSPs have enough samples providing complete data for phosphorylation level and all corresponding kinase protein expression level, which is larger than number of inputs for the the regression model1: pho_sub~pro_kin1+pro_kin2+....","\n", sep = " "))
cat(paste("(",temp$phosite_count," phosphosites,",temp$sub_count," substrates,",temp$kin_count," kinases)","\n\n", sep = " "))
cat(paste("(kinase == substrate: ",temp$self_ratio,"; kinase !=substrate: ",temp$other_ratio,")","\n", sep = " "))

temp <- KSP_count(which(!is.nan(unlist(k_s_fit$P_pro_kin2))))
cat(paste(temp$ksp_count," above KSPs have enough samples providing complete data for substrate protein expression level in addition, which is larger than number of inputs for the the regression model2: pho_sub~pro_kin1+pro_kin2+....+pro_sub","\n", sep = " "))
cat(paste("(kinase == substrate: ",temp$self_ratio,"; kinase !=substrate: ",temp$other_ratio,")","\n", sep = " "))
cat(paste("(",temp$phosite_count," phosphosites,",temp$sub_count," substrates,",temp$kin_count," kinases)","\n\n", sep = " "))

temp <- KSP_count(which(!is.nan(unlist(k_s_fit$P_pro_kin3))))
cat(paste(temp$ksp_count," above KSPs have enough samples providing complete data for grouped kinase phosphorylation level in addition, which is larger than number of inputs for the the regression model3:pho_sub~pro_kin1+pro_kin2+....+pro_sub+pho_g_kin1+pho_g_kin2+...","\n", sep = " "))
cat(paste("(kinase == substrate: ",temp$self_ratio,"; kinase !=substrate: ",temp$other_ratio,")","\n", sep = " "))
cat(paste("(",temp$phosite_count," phosphosites,",temp$sub_count," substrates,",temp$kin_count," kinases)","\n\n", sep = " "))

cat(paste("\n", sep = " "))
sink()



# construct table for bubble chart then do analysis over gain for the other cancer dataset -----------------------------------------------------
k_s_fit$Cancer <- cancer
k_s_2can <- rbind(k_s_2can,k_s_fit)
table <- k_s_fit[!is.nan(unlist(k_s_fit$P_pro_kin1)),] #only keep those KSPs producing regression result
table_2can <- rbind(table_2can,table)

# bubble chart module -----------------------------------------------------
setwd("~/Desktop/Ding Lab/Figures")
library(ggplot2)
library(ggtern)
require(data.table)
library(reshape)

# keep significant ones, divide self-regulated and others
table_2can$pair <- paste(table_2can$GENE,":",table_2can$SUB_GENE,":",table_2can$SUB_MOD_RSD,sep="")
table_sig <- table_2can[table_2can$FDR_pro_kin1 <= 0.05,]
table_sig_self <- table_sig[table_sig$self,]
table_sig_other <- table_sig[!table_sig$self,]

temp <- c("GENE","SUB_GENE","SUB_MOD_RSD","P_pro_kin1","P_pro_kin2","P_pro_sub2","P_pro_kin3","P_pro_sub3","P_pho_g_kin3","coef_pro_kin1","coef_pro_kin2","coef_pro_sub2","coef_pro_kin3","coef_pro_sub3","coef_pho_g_kin3","self","Cancer","pair")
# choose the data for plotting
t0 <- table_sig_self# choose table_sig_other or table_sig_self

dm1 <- melt(t0[,c(temp,paste("FDR_",name,sep = ""))], id=temp)
colnames(dm1) <- c(temp, "variable", "FDR")
dm1$model <- vector(mode = "numeric", length = nrow(dm1))
dm1$model[dm1$variable=="FDR_pro_kin1"] <- 1
dm1$model[dm1$variable=="FDR_pro_kin2" | dm1$variable=="FDR_pro_sub2"] <- 2
dm1$model[dm1$variable=="FDR_pro_kin3" | dm1$variable=="FDR_pro_sub3" | dm1$variable=="FDR_pho_g_kin3"] <- 3

dm2 <- melt(t0[,c(temp,paste("size_dat",1:3,sep = ""))], id=temp)
colnames(dm2) <- c(temp, "variable2", "sample_size")
dm2$model <- vector(mode = "numeric", length = nrow(dm2))
dm2$model[dm2$variable2=="size_dat1"] <- 1
dm2$model[dm2$variable2=="size_dat2"] <- 2
dm2$model[dm2$variable2=="size_dat3"] <- 3
molten <- merge(dm1, dm2)


## plotting
p = ggplot(molten,aes(x=variable, y=pair))# make this the original ethni
p = p + facet_grid(.~Cancer+model, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
#p = p + facet_grid(GENE~Cancer+model)#, space = "free", scales = "free")
p = p + geom_point(aes(color = sample_size , size =-log(FDR))) + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))#element_text(colour="black", size=14))
p


