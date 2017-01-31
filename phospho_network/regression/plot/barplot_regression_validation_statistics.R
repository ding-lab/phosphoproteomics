# Yige Wu @ WashU 2016 Dec
# barplot the validation statistics for regression results

# choose kinase/phosphotase, cancer , significance level -----------------------------------------------
protein <- "kinase"
sig <- 0.05
cancer <- "BRCA"

# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)
library(readr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# input k_s_table-------------------------------------------------------------------
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


# input regression processed data -----------------------------------------
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",protein,"_substrate_regression.txt",sep = ""))
## only use data of BRCA and model1
table1 <- table_2can[table_2can$Cancer== cancer & table_2can$model=="pho_sub~pro_kin",]

# different standard of validation for kinase and phosphotase -------------
if ( protein == "kinase" ) {
  valid_type <- "pos_sig"
}
if ( protein == "phosphotase" ) {
  valid_type <- "neg_sig"
}


# initiate ----------------------------------------------------------------
x <- vector(mode = "numeric", length = length(unique_kinase) )
temp <- data.frame(matrix(rep(x,4), ncol=4, byrow=T))
valid_trans <- cbind(unique_kinase,temp)
colnames(valid_trans) <- c("kinase","pos_sig","pos_insig","neg_sig","neg_insig")
rownames(valid_trans) <- unique_kinase
valid_cis <- valid_trans


# make tables for trans pairs ---------------------------------------------
for( kinase in unique_kinase) {
  temp <- table1[table1$KINASE==kinase,]
  valid_trans[kinase,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & !temp$self))
  valid_trans[kinase,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & !temp$self))
  valid_trans[kinase,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & !temp$self))
  valid_trans[kinase,"neg_insig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & !temp$self))
}
if ( protein == "phosphotase") {
  valid_trans <- valid_trans[,c("kinase","neg_sig","neg_insig","pos_sig","pos_insig")]
}
valid_trans$all_count <- rowSums(valid_trans[,c("pos_sig","pos_insig","neg_sig","neg_insig")])
valid_trans$valid_count <- valid_trans[,valid_type]
valid_trans$valid_ratio <- valid_trans$valid_count/valid_trans$all_count
table <- melt(valid_trans[!is.na(valid_trans$valid_ratio),],id=c("kinase","all_count","valid_count","valid_ratio"))
colnames(table) <- c("kinase","all_count","valid_count","valid_ratio","coef_FDR","count")
table$KINASE <- reorder(table$kinase, table$valid_count)


# write out table for trans ---------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",cancer,'_',protein,"_trans_regression_validation_statistics.txt", sep="")
write.table(valid_trans, file=tn, quote=F, sep = '\t', row.names = FALSE)


# barplot ranking validation count for trans----------------------------------------
ggplot() +
  geom_bar(data=table, aes(y = count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab(protein)+ylab("number of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_trans_validation_count_ranking_per_gene.png',sep ="")
ggsave(file=fn, height=12, width=8)

# barplot validation ratio (by count ranking) for trans----------------------------------------
ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab(protein)+ylab("number of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_trans_validation_ratio_per_gene.png',sep ="")
ggsave(file=fn, height=12, width=8)


# scatterplot validation ratio for trans ----------------------------------------
p <- ggplot(valid_trans, aes(all_count, valid_ratio, label = rownames(valid_trans)))
p + geom_text(check_overlap = TRUE)
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_trans_all_count~valid_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)


# make tables for cis pairs -----------------------------------------------
for( kinase in unique_kinase) {
  temp <- table1[table1$KINASE==kinase,]
  valid_cis[kinase,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & temp$self))
  valid_cis[kinase,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & temp$self))
  valid_cis[kinase,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & temp$self))
  valid_cis[kinase,"neg_insig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & temp$self))
}
if ( protein == "phosphotase") {
  valid_cis <- valid_cis[,c("kinase","neg_sig","neg_insig","pos_sig","pos_insig")]
}
valid_cis$all_count <- rowSums(valid_cis[,c("pos_sig","pos_insig","neg_sig","neg_insig")])
valid_cis$valid_count <- valid_cis[,valid_type]
valid_cis$valid_ratio <- valid_cis$valid_count/valid_cis$all_count
table <- melt(valid_cis[!is.na(valid_cis$valid_ratio),],id=c("kinase","all_count","valid_count","valid_ratio"))
colnames(table) <- c("kinase","all_count","valid_count","valid_ratio","coef_FDR","count")
table$KINASE <- reorder(table$kinase, table$valid_count)

# write out table for trans ---------------------------------------------------------
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",cancer,'_',protein,"_cis_regression_validation_statistics.txt", sep="")
write.table(valid_cis, file=tn, quote=F, sep = '\t', row.names = FALSE)

# barplot ranking validation count for cis----------------------------------------
ggplot() +
  geom_bar(data=table, aes(y = count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab(protein)+ylab("number of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_cis_gene_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

# barplot validation ratio (by count ranking) for cis----------------------------------------
ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab(protein)+ylab("number of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_cis_gene_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

# scatterplot validation ratio for trans ----------------------------------
p <- ggplot(valid_cis, aes(all_count, valid_ratio, label = rownames(valid_cis)))
p + geom_text(check_overlap = TRUE)
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_cis_all_count~valid_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)
# trans-residue -----------------------------------------------------------
## trans-residue
## initiate
x <- vector(mode = "numeric", length = length(3) )
temp <- data.frame(matrix(rep(x,4), ncol=4, byrow=T))
valid_rsd_trans <- cbind(c("S","T","Y"),temp)
colnames(valid_rsd_trans) <- c("residue","pos_sig","pos_insig","neg_sig","neg_insig")
rownames(valid_rsd_trans) <- c("S","T","Y")
valid_rsd_cis <- valid_rsd_trans

## trans
for(rsd in c("S","T","Y") ) {
  temp <- table_2can[grepl(rsd,table_2can$SUB_MOD_RSD),]
  valid_rsd_trans[rsd,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & !temp$self))
  valid_rsd_trans[rsd,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & !temp$self))
  valid_rsd_trans[rsd,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & !temp$self))
  valid_rsd_trans[rsd,"neg_insig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & !temp$self))
}
if ( protein == "phosphotase") {
  valid_rsd_trans <- valid_rsd_trans[,c("residue","neg_sig","neg_insig","pos_sig","pos_insig")]
}
valid_rsd_trans$all_count <- rowSums(valid_rsd_trans[,c("pos_sig","pos_insig","neg_sig","neg_insig")])
valid_rsd_trans$valid_count <- valid_rsd_trans[,valid_type]
valid_rsd_trans$valid_ratio <- valid_rsd_trans$valid_count/valid_rsd_trans$all_count
table <- melt(valid_rsd_trans[!is.na(valid_rsd_trans$valid_ratio),],id=c("residue","all_count","valid_count","valid_ratio"))
colnames(table) <- c("residue","all_count","valid_count","valid_ratio","coef_FDR","count")
table$RESIDUE <- reorder(table$residue, table$valid_count)

ggplot() +
  geom_bar(data=table, aes(y = count, x = RESIDUE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("residue")+ylab("number of phosphosite residues with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_trans_residue_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = RESIDUE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("residue")+ylab("ratio of phosphosite residues with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_trans_residue_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

# cis-residue -------------------------------------------------------------
## cis-residue
for(rsd in c("S","T","Y") ) {
  temp <- table_2can[grepl(rsd,table_2can$SUB_MOD_RSD),]
  valid_rsd_cis[rsd,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & temp$self))
  valid_rsd_cis[rsd,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & temp$self))
  valid_rsd_cis[rsd,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & temp$self))
  valid_rsd_cis[rsd,"neg_insig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & temp$self))
}
if ( protein == "phosphotase") {
  valid_rsd_cis <- valid_rsd_cis[,c("residue","neg_sig","neg_insig","pos_sig","pos_insig")]
}
valid_rsd_cis$all_count <- rowSums(valid_rsd_cis[,c("pos_sig","pos_insig","neg_sig","neg_insig")])
valid_rsd_cis$valid_count <- valid_rsd_cis[,valid_type]
valid_rsd_cis$valid_ratio <- valid_rsd_cis$valid_count/valid_rsd_cis$all_count
table <- melt(valid_rsd_cis[!is.na(valid_rsd_cis$valid_ratio),],id=c("residue","all_count","valid_count","valid_ratio"))
colnames(table) <- c("residue","all_count","valid_count","valid_ratio","coef_FDR","count")
table$RESIDUE <- reorder(table$residue, table$valid_count)

ggplot() +
  geom_bar(data=table, aes(y = count, x = RESIDUE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("residue")+ylab("number of phosphosite residues with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_cis_residue_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = RESIDUE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("residue")+ylab("ratio of phosphosite residues with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",cancer,'_',protein,'_cis_residue_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)






