# Yige Wu @ WashU 2016 Nov
# try to find samples with potential activated kinases, samples and kinase cooccurence

# choose kinase or phosphotase, outlier threshold and least sample number-------------------------
protein <- "kinase"
cancer <- "BRCA"
center <- 0.5
sig <- 0.05

# library -----------------------------------------------------------------
library(reshape)
library(stringr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# input k_s_table-------------------------------------------------------------------
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

# input according to cancer type-------------------------------------------------------------------
  if (cancer == "BRCA") {
    # BRCA
    BRCA_pro_f = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized_max10NA.txt",sep="")
    pro_data <- read.delim(BRCA_pro_f)
    BRCA_pho_f = paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv",sep="")
    pho_data = read.delim(BRCA_pho_f)
    BRCA_pho_g = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized_max10NA.txt",sep="")
    pho_gdata = read.delim(BRCA_pho_g)
    colx <- 78 # the column of protein name
    clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep = ""))
    somatic <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep = ""))
    somatic_aa <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted_amino_acid.txt",sep = ""))
    cnv_data <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_CNV_formatted_normalized.txt",sep = ""))
    pro_outlier <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/druggable_outlier/2016-03-02/2016-03-02_KH_BRCA druggable proteome normalized_outlier_score_table.txt",sep = ""))
  }
  
  # ordering the columns by sample name
  pro_data <- pro_data[,order(names(pro_data))]
  pho_data <- pho_data[,order(names(pho_data))]
  pho_gdata <- pho_gdata[,order(names(pho_gdata))]
  somatic <- somatic[,order(names(somatic))]
  somatic_aa <- somatic_aa[,order(names(somatic_aa))]
  pro_outlier <- pro_outlier[,order(names(pro_outlier))]; rownames(pro_outlier) <- pro_outlier$X
  clinical <- clinical[,order(names(clinical))]
  cnv_data <- cnv_data[,order(names(cnv_data))]
  
  #split the SUBSTRATE and SUB_MOD_RSD in the first column
  pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))
  
  #covert the SUB_MOD_RSD from lowercase to uppercase
  pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
  colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")

  sample_names <- colnames(pho_data[,-colx])

  table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))
  table_brca_trans_valid <- table_2can[!table_2can$self & table_2can$FDR_pho_kin<=sig & table_2can$coef_pho_kin>0 & table_2can$Cancer=="BRCA",]

## add collapsed phosphorylation level and score------------------------------------------
  rownames(pho_gdata) <- pho_gdata$X
  
  kin_phos_c <- melt(pho_gdata,id="X")
  colnames(kin_phos_c) <- c("kinase","sample","kin_phos_c")
  overlap <- kin_phos_c
  
  pho_cscore <- pho_gdata
  for (kinase in pho_cscore$X) {
    temp <- pho_gdata[kinase,sample_names]
    IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
    pho_cscore[kinase,sample_names] = ( temp - quantile(temp, probs=center, na.rm=T))/IQR
  }
  
  kin_phos_cscore <- melt(pho_cscore, id="X")
  colnames(kin_phos_cscore) <- c("kinase","sample","kin_phos_cscore")
  overlap <- merge(overlap,kin_phos_cscore,all.x = T)

## add averaged upstream kinase collapsed phosphorylation level and score(exclude itself)-------------------------------
  # initiate
  unique_substrate = unique(as.vector(k_s_table$SUB_GENE)) 
  x <- vector(mode = "numeric", length = length(unique_substrate))+NaN
  temp <- data.frame(matrix(rep(x,ncol(pho_data)-1), ncol=ncol(pho_data)-1, byrow=T))
  up_pho_c <- cbind(unique_substrate,temp)
  colnames(up_pho_c) <- c("KINASE", sample_names)
  rownames(up_pho_c) <- up_pho_c$KINASE
  up_pho_cscore <- up_pho_c
  
  # fill in up_pho_c
  n_up_pho <- c()
  for (gene in unique_substrate){
    k_k_s_table <- unique(table_brca_trans_valid$KINASE[table_brca_trans_valid$SUBSTRATE==gene])
    row <- c()
    for (kinase in k_k_s_table) {
      row <- c(row,which(pho_gdata$X==kinase))
    }
    n_up_pho <- c(n_up_pho,length(row))
    for (sam in sample_names) {
      # for each cell, average phosphosites phosphrylation level for all the substrates for this kinase
      up_pho_c[gene, sam] <- mean(pho_gdata[row,sam], na.rm = TRUE)
    }
  }
  n_up_pho <- data.frame(n_up_pho); n_up_pho$kinase <- unique_substrate
  overlap <- merge(overlap,n_up_pho,all.x = T)
  
  # fill in up_pho_cscore
  for (kinase in up_pho_c$KINASE) {
    temp <- up_pho_c[kinase,sample_names]
    IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
    up_pho_cscore[kinase,sample_names] = ( temp - quantile(temp, probs=center, na.rm=T))/IQR
  }
  
  # melt the up_pho_c & up_pho_cscore
  up_phos_cscore <- melt(up_pho_cscore, id="KINASE")
  colnames(up_phos_cscore) <- c("kinase","sample","up_phos_cscore")
  overlap <- merge(overlap,up_phos_cscore,all.x = T)
  
  up_phos_c <- melt(up_pho_c, id="KINASE")
  colnames(up_phos_c) <- c("kinase","sample","up_phos_c")
  overlap <- merge(overlap,up_phos_c,all.x = T)
  


## add avaraged collapsed substrate phosphorylation score excluding itself ------------------
  # initiate
  unique_kinase <- as.vector(unique(k_s_table$GENE) )
  x <- vector(mode = "numeric", length = length(unique_kinase))+NaN
  temp <- data.frame(matrix(rep(x,ncol(pho_data)-1), ncol=ncol(pho_data)-1, byrow=T))
  sub_pho_ctrans <- cbind(unique_kinase,temp)
  colnames(sub_pho_ctrans) <- c("KINASE", sample_names)
  rownames(sub_pho_ctrans) <- sub_pho_ctrans$KINASE
  sub_pho_ctransscore <- sub_pho_ctrans
  
  # fill in sub_pho_ctrans
  n_sub_pho <- c()
  for (kinase in unique_kinase){
    #for (kinase in "EIF2AK1") {
    k_k_s_table <- unique(table_brca_trans_valid$SUBSTRATE[table_brca_trans_valid$KINASE==kinase])
    row <- c()
    for (substrate in k_k_s_table) {
      row <- c(row,which(pho_gdata$X==substrate))
    }
    n_sub_pho <- c(n_sub_pho,length(row))
    for (sam in sample_names) {
      # for each cell, average phosphosites phosphrylation level for all the substrates for this kinase
      sub_pho_ctrans[kinase, sam] <- mean(pho_gdata[row,sam], na.rm = TRUE)
    }
  }
  n_sub_pho <- data.frame(n_sub_pho);n_sub_pho$kinase <- unique_kinase
  overlap <- merge(overlap,n_sub_pho,all.x = T)
  
  for (kinase in unique_kinase) {
    #for (kinase in "ABL1") {
    temp <- sub_pho_ctrans[kinase,sample_names]
    IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
    sub_pho_ctransscore[kinase,sample_names] = (temp - quantile(temp, probs=center, na.rm=T))/IQR
  }
  
  # melt the sub_pho_ctrans & sub_pho_ctransscore
  sub_phos_ctransscore <- melt(sub_pho_ctransscore, id="KINASE")
  colnames(sub_phos_ctransscore) <- c("kinase","sample","sub_phos_ctransscore")
  overlap <- merge(overlap,sub_phos_ctransscore,all.x = T)
  
  sub_phos_ctrans <- melt(sub_pho_ctrans, id="KINASE")
  colnames(sub_phos_ctrans) <- c("kinase","sample","sub_phos_ctrans")
  overlap <- merge(overlap,sub_phos_ctrans,all.x = T)
  
  
  
  
## add mutation ------------------------------------------------------------
rownames(somatic) <- somatic$X
mutation <- melt(somatic,id="X")
colnames(mutation) <- c("kinase","sample","mutation")
overlap <- merge(overlap,mutation,all.x = T)
overlap$mut <- grepl("ins",overlap$mutation) | grepl("del",overlap$mutation) | grepl("sense",overlap$mutation) | grepl("splice",overlap$mutation)

## add somatic aa change ---------------------------------------------------
rownames(somatic_aa) <- somatic_aa$X
aa <- melt(somatic_aa,id="X")
colnames(aa) <- c("kinase","sample","aa")
overlap <- merge(overlap, aa, all.x = T)


## add protein expression level(actual value) ------------------------------
pro_level <- melt(pro_data,id="X")
colnames(pro_level) <- c("kinase","sample","pro_level")
overlap <- merge(overlap, pro_level, all.x = T)

## add protein expression score --------------------------------------------
rownames(pro_data) <- pro_data$X
proscore <- pro_data
for (gene in proscore$X ) {
  temp <- pro_data[gene,sample_names]
  IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
  proscore[gene,sample_names] = (temp - quantile(temp, probs=center, na.rm=T))/IQR
}
pro_score <- melt(proscore,id="X")
colnames(pro_score) <- c("kinase","sample","pro_score")
overlap <- merge(overlap, pro_score, all.x = T)


## add CNV score --------------------------------------------
rownames(cnv_data) <- cnv_data$X
cnvscore <- cnv_data
for (gene in cnvscore$X ) {
  temp <- cnv_data[gene,sample_names]
  IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
  cnvscore[gene,sample_names] = (temp - quantile(temp, probs=center, na.rm=T))/IQR
}
cnv_score <- melt(cnvscore,id="X")
colnames(cnv_score) <- c("kinase","sample","cnv_score")
overlap <- merge(overlap, cnv_score, all.x = T)

cnv_level <- melt(cnv_data, id="X")
colnames(cnv_level) <- c("kinase","sample","cnv_level")
overlap <- merge(overlap, cnv_level, all.x = T)



# write out tables --------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/",cancer,"_ks_score_validated_ks_only_median_scored_",protein,"_table.txt", sep="")
write.table(overlap, file=tn, quote=F, sep = '\t', row.names = FALSE)



