# Yige Wu @ WashU 2016 Nov
# try to find samples with potential activated kinases, samples and kinase cooccurence

# choose kinase or phosphotase, outlier threshold and least sample number-------------------------
out_thres <- 1.5 #threshold for outlier
protein <- "kinase"
cancer <- "BRCA"

# library -----------------------------------------------------------------
library(reshape)
library(stringr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


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

# input according to cancer type-------------------------------------------------------------------
  if (cancer == "BRCA") {
    # BRCA
    BRCA_pro_f = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized_max10NA.txt",sep="")
    pro_data <- read.delim(BRCA_pro_f)
    BRCA_pho_f = paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv",sep="")
    pho_data = read.delim(BRCA_pho_f)
    BRCA_pho_g = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_by_PRO_formatted_normalized_max10NA.txt",sep="")
    pho_gdata = read.delim(BRCA_pho_g)
    colx <- 78 # the column of protein name
    clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep = ""))
    somatic <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep = ""))
    somatic_aa <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted_amino_acid.txt",sep = ""))
    pro_outlier <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/druggable_outlier/2016-03-02/2016-03-02_KH_BRCA druggable proteome normalized_outlier_score_table.txt",sep = ""))
  }
  
  
  # ordering the columns by sample name
  pro_data <- pro_data[,order(names(pro_data))]
  pho_data <- pho_data[,order(names(pho_data))]
  pho_gdata <- pho_gdata[,order(names(pho_gdata))]#order the grouped phospho data
  somatic <- somatic[,order(names(somatic))]
  somatic_aa <- somatic_aa[,order(names(somatic_aa))]
  pro_outlier <- pro_outlier[,order(names(pro_outlier))]; rownames(pro_outlier) <- pro_outlier$X
  clinical <- clinical[,order(names(clinical))]
  
  #split the SUBSTRATE and SUB_MOD_RSD in the first column
  pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))
  
  #covert the SUB_MOD_RSD from lowercase to uppercase
  pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
  colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")


# initiate -------------------------------------------------------------------
# first column == kinase name, column names == sample names
k_s_table_trans <- k_s_table
unique_kinase = unique(as.vector(k_s_table_trans$GENE)) 
x <- vector(mode = "numeric", length = length(unique_kinase))+NaN
temp <- data.frame(matrix(rep(x,ncol(pho_data)-1), ncol=ncol(pho_data)-1, byrow=T))
k_p_table <- cbind(unique_kinase,temp)
sample_names <- colnames(pho_data[-colx])
colnames(k_p_table) <- c("KINASE", sample_names)
rownames(k_p_table) <- k_p_table$KINASE
nsample <- length(sample_names)

# pre-process the k_s_table_trans to remove those rows in which substrate==kinase -------------------------------------------------------------------
start.time <- Sys.time()
for (kinase in unique_kinase){
  k_k_s_table <- unique(as.vector(k_s_table_trans$SUB_GENE[k_s_table_trans$GENE==kinase]))
  row <- c()
  for (substrate in k_k_s_table) {
    row <- c(row,which(pho_rsd_split$SUBSTRATE==substrate))
  }
  for (sam in 1:nsample) {
    # for each cell, average phosphosites phosphrylation level for all the substrates for this kinase
    k_p_table[k_p_table$KINASE==kinase, sam+1] <- mean(pho_data[row,sam], na.rm = TRUE)
    
  }
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

k_p_score <- k_p_table
for (kinase in k_p_table$KINASE) {
  #for (kinase in "ABL1") {
  temp <- k_p_table[kinase,sample_names]
  IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
  k_p_score[kinase,sample_names] = ( temp - quantile(temp, probs=0.75, na.rm=T))/IQR
}

# melt the k_p_table, add variable Her+/- according to sample name --------
k_p_her <- melt(k_p_table, id="KINASE")
Her2 <- vector(mode = "numeric", length = length(k_p_her))+NA
for (sam in sample_names[which(clinical[1,]=="Her2")]) {
  Her2[k_p_her$variable==sam] <- "Her2+"
}
Her2[is.na(Her2)] <- "Her2-"
k_p_her$Her2_status <- Her2
colnames(k_p_her) <- c("kinase","sample","ave_sub_phos","Her_status")

# add pro_outlier into k_p_her --------------------
## pre-process the pro_outlier to have a logical matrix
x <- vector(mode = "logical", length = nrow(pro_outlier))
temp <- data.frame(matrix(rep(x,ncol(pro_outlier)-1), ncol=ncol(pro_outlier)-1, byrow=T))
is.pro_outlier <- cbind(pro_outlier$X,temp)
colnames(is.pro_outlier) <- c("kinase", sample_names)
rownames(is.pro_outlier) <- is.pro_outlier$kinase

## loop around the kinases
for (kinase in pro_outlier$X) {
#for (kinase in "ABL1") {
  temp <- pro_outlier[kinase,sample_names]
  IQR = quantile(temp, probs=0.75, na.rm=T) - quantile(temp, probs=0.25, na.rm=T) 
  is.pro_outlier[kinase,sample_names] = ( pro_outlier[kinase,sample_names] >= quantile(temp, probs=0.75, na.rm=T) + out_thres*IQR)
}


## add pro_outlier into k_p_her
is.pro.outlier <- melt(is.pro_outlier, id="kinase")
colnames(is.pro.outlier) <- c("kinase","sample","pro_outlier_status")
pro.level <- melt(pro_outlier, id="X")
colnames(pro.level) <- c("kinase","sample","pro_level")
pro.level <- merge(pro.level,is.pro.outlier)

score <- melt(k_p_score, id="KINASE")
colnames(score) <- c("kinase","sample","score")
overlap <- merge(k_p_her,score)
overlap <- merge(overlap,pro.level, all.x = T)
overlap$pro_outlier_status[is.na(overlap$pro_outlier_status)] <- FALSE


## add mutatin into it 
rownames(somatic) <- somatic$X
mut <- c()
for (kinase in somatic$X) {
  #for (kinase in "A1BG") {
  mut <- rbind(mut,t(grepl("missense",unlist(somatic[kinase,sample_names]))))
}
mut <- data.frame(mut)
colnames(mut) <- sample_names
mut$kinase <- somatic$X
missense <- melt(mut,id="kinase")
colnames(missense) <- c("kinase","sample","missense")

rownames(somatic_aa) <- somatic_aa$X
aa <- melt(somatic_aa,id="X")
colnames(aa) <- c("kinase","sample","aa")
overlap <- merge(overlap, aa, all.x = T)

overlap <- merge(overlap, missense, all.x = T)
overlap$missense[is.na(overlap$missense)] <- FALSE
overlap$missense.pro_outlier <- paste(overlap$missense,overlap$pro_outlier_status,sep = ",")

# write out tables --------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/",protein,"_",cancer,"_reverse_marking_out_thers_",out_thres,".txt", sep="")
write.table(overlap, file=tn, quote=F, sep = '\t', row.names = FALSE)

