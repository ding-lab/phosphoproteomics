# Yige Wu @ WashU 2016 Nov
# do the regression on 4 BRCA subtype groups repectively

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

# choose kinase or phosphotase, significance level, outlier threshold and least sample number-------------------------
sig <- 0.05 # significance level
out_thres <- 1.5 #threshold for outlier
least_samples <- 5# least number of samples with complete data for each model
protein <- "kinase"

# input regardless to cancer type, choose between kinase or phosphotase-------------------------------------------------------------------
if ( protein == "kinase" ) {
  ### read in the kinase/substrate table/ phosphorylation data ###
  K_S_f = paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep="")
  k_s_table = read.delim(K_S_f)
}

if ( protein == "phosphotase" ) {
  ### read in the phosphotase/substrate table/ phosphorylation data ### 
  k_s_table <- read.csv(paste(baseD,"pan3can_shared_data/Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.csv",sep = ""))
  colnames(k_s_table) <- c("Phosphatase_UniProtAC_human","GENE","Substrate_UniProtAC_ref","SUB_GENE","Substrate_Type","DephosphoSite","BioassayType", "PubMed_ID_rev")
}

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

unique_kinase <- unique(k_s_table$GENE)

# looping subtypes,run following sections all at once before plot module -----------------------------------------------------------
cancer = "BRCA"
pro_raw <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized_max10NA.txt",sep=""))
pho_raw = read.delim(paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv",sep=""))
pho_g_raw = read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_by_PRO_formatted_normalized_max10NA.txt",sep=""))
clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))

for (cohort in c("Her2","LumA","LumB","Basal")) {
#for (cohort in c("Her2")) {
# input according to subtype-------------------------------------------------------------------
  # BRCA
  subtype_sample <- na.omit(colnames(clinical)[clinical[1,]==cohort])
  pro_data <- pro_raw[,c(subtype_sample,"X")]
  pho_data <- pho_raw[,c(subtype_sample,"X")]
  pho_gdata <- pho_g_raw[,c(subtype_sample,"X")]

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

# looping for model1 -----------------------------------------------------------------
  for (kinase in unique_kinase){
    #for (kinase in "ERBB2"){#test
    # find protein expression level for the kinase
    pro_kinase <- pro_data[pro_data$X == kinase,subtype_sample]
    
    if(nrow(pro_kinase) > 0 ){
      
      # # find grouped phosphorylation level for kinase
      # pho_kinase_g <- pho_gdata[pho_gdata$X == kinase,-colx]
      # pho_kin_g_norm <- range01(unlist(pho_kinase_g),na.rm = T)
      
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
          pho_sub <- pho_data[i,subtype_sample]
          
          
          #normalize phospho level of substrate and protein expression level of the kinase to 0 to 1
          pho_sub_norm <- range01(unlist(pho_sub),na.rm = T)
          pro_kin_norm <- range01(unlist(pro_kinase),na.rm = T)
          
          #prepare regression data for model1
          data1 <- data.frame(pho_sub_norm,pro_kin_norm)
          
          size <- nrow(data1[complete.cases(data1),])
          if( size > least_samples ){
            # fit regression model1: pho_substrate ~ a*pro_kinase + k
            fit1 <- glm(pho_sub_norm ~ pro_kin_norm,data = data1, family=gaussian())
            
            # record the kinase name, kinase expression level, substrate name, SUB_MOD_RSD, phophorylation level into the list1
            list1[[length(list1)+1]] <- list(KINASE=kinase,SUBSTRATE=substrate,SUB_MOD_RSD=sub_mod_rsd,Pvalue=c(coef(summary(fit1))[2,4]),Coef_pro_kin=fit1$coefficients[2],Size=size)
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
  # table1$P_pro_sub <- NA
  # table1$P_pho_kin <- NA
  table1$coef_pro_kin <- sapply(list1, "[[", "Coef_pro_kin")
  # table1$coef_pro_sub <- NA
  # table1$coef_pho_kin <- NA
  table1$pair <- paste(table1$KINASE,":",table1$SUBSTRATE,":",table1$SUB_MOD_RSD,sep="")
  
  # combine table
  # table <- rbind(table1,table2,table3,table4)

  # mark cancer and self-regulation
  table1$cohort <- cohort
  
  if ( cohort == "Her2" ) {
    table_Her2 <- table1
  } 
  if ( cohort == "LumA" ) {
    table_LumA <- table1
  } 
  if ( cohort == "LumB" ) {
    table_LumB <- table1
  }
  if ( cohort == "Basal" ) {
    table_Basal <- table1
  } 
}

table <- rbind(table_Her2,table_LumA,table_LumB,table_Basal)
table$Cancer <- cancer
table$self <- as.character(table$KINASE) == as.character(table$SUBSTRATE)
table$SELF <- "trans"; table$SELF[table$self] <- "cis"

## adjust p-values to FDR
for (mod in c("pho_sub~pro_kin")) {
  for (cohort in c("Her2","LumA","LumB","Basal")) {
    for(self in c(TRUE,FALSE)) {
      for(coln in c("pro_kin")) {#adjust pvalues for each variable
        row <- (table$self==self) & (table$model==mod) & (table$cohort==cohort)
        table[row,paste("FDR_",coln,sep = "")] <-p.adjust(table[row,paste("P_",coln,sep = "")],method = "fdr")
      }
    }
  }
}

# write out tables --------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",cancer,"_subtypes",protein,"_substrate_regression.txt", sep="")
write.table(table, file=tn, quote=F, sep = '\t', row.names = FALSE)

table_subtype <- table

# volcano plotting module -------------------------------------------------
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
mode <- "pho_sub~pro_kin"

table_subtype$coef_pro_kin_filtered = remove_outliers(table_subtype$coef_pro_kin)
table_subtype_outlier_removed_m = table_subtype[!is.na(table_subtype$coef_pro_kin_filtered) & table_subtype$model == mode,]

if ( protein == "kinase") {
  plot_fdr_scale <- -log10(sig)
}
if ( protein == "phosphotase") {
  plot_fdr_scale <- 2
}
p = ggplot(table_subtype_outlier_removed_m,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin)))
p = p + facet_grid(SELF~cohort,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.05)
p = p + geom_text(aes(label= ifelse(-log10(FDR_pro_kin)>plot_fdr_scale, pair, NA)),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/BRCA_subtype_',protein,'_substrate_volcano.pdf',sep ="")
ggsave(file=fn, height=6, width=12, useDingbats=FALSE)

# bubble chart ---------------------
table_sig <- table_subtype[table_subtype$FDR_pro_kin <= sig,]
# choose the parameters for bubble chart
top <- 50 # choose top n rows for FDR_pro_kin for model1
c <- "LumB" #choose in which cancer extract the top n rows
for (cis in c(TRUE,FALSE)) { # loop around cis and trans
  t0 <- table_sig[table_sig$self==cis,]
  # sort by FDR_pro_kin
  t1 <- t0[t0$cohort==c,]
  t1 <- t1[order(t1$FDR_pro_kin),] 
  
  ## corresponding results for other two models are extracted and ordered
  rows <- c()
  for(i in 1:top){
    r <- unlist(which(t0$pair==t1$pair[i]))
    rows <- c(rows,r)
  }
  table <- t0[rows,]
  #table_2can_brca_top = table_2can[table_2can$pair %in% table$pair,]
  
  ## actual plotting
  lim = max(max(table$coef_pro_kin),min(table$coef_pro_kin))
  if(nrow(table) >0) {
    p = ggplot(table,aes(x=model, y=pair))# make this the original ethni
    p = p + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
    p = p + geom_point(aes(fill=coef_pro_kin, size =-log10(FDR_pro_kin), color=ifelse(sig, "black",NA)),pch=21) 
    p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
    p = p + scale_colour_manual(values=c("black",NA))
    p = p + theme_bw() #+ theme_nogrid()
    p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
    p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
    p
    fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/',c,'_subtypes_cis_',cis,'_fdr',sig,'_top_',top,'.pdf',sep ="")
    ggsave(file=fn, height=12, width=8, useDingbats=FALSE)
  }
}


