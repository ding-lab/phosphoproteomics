# Yige Wu @ WashU 2016 Nov
# look at correlations of kinase and downstream substrates phosphorylation status


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

# looping cancer,run following sections all at once before plot module -----------------------------------------------------------
for (cancer in c("BRCA","OV")) {
#for (cancer in c("BRCA")) {
  
# input according to cancer type-------------------------------------------------------------------
  if (cancer == "BRCA") {
    # BRCA
    pro_raw <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized_max10NA.txt",sep=""))
    pro_data <- pro_raw[,-1]
    
    pho_raw <- read_delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_RPPA_formatted_wGpos_forR.txt",sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
    pho_raw <- pho_raw[!is.na(pho_raw$genomic_pos) & pho_raw$genomic_pos != ".",]
    pho_data <- pho_raw[,-c(1,2)]
    pho_data <- pho_data[,colnames(pro_data)]
    
    transvar <- read_delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_RPPA_formatted_transvar_output.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
  }
  
  if ( cancer == "OV" ) {
    #OV
    pho_raw <- read_delim(paste(baseD,"pan3can_shared_data/OV/OV_RPPA_formatted_wGpos_forR.txt",sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
    pho_raw <- pho_raw[!is.na(pho_raw$genomic_pos) & pho_raw$genomic_pos != ".",]
    pho_data <- pho_raw[,-c(1,2)]
    
    pro_raw <- read.delim(paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized_max10NA.txt",sep=""))
    pro_data <- pro_raw[,-1]
    
    temp <-  str_split_fixed(colnames(pho_raw[,-c(1,2)]),"-",2)
    temp <- paste("X",temp[,1],".",temp[,2],"_PNNL",sep="")
    colnames(pho_data) <- temp
    temp2 <- intersect(temp, colnames(pro_data))
    
    pro_data <- pro_data[,temp2]
    pho_data <- pho_data[,temp2]
    
    transvar <- read_delim(paste(baseD,"pan3can_shared_data/OV/OV_RPPA_formatted_transvar_output.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
  }

  transvar_split <- data.frame(str_split_fixed(transvar$`coordinates(gDNA/cDNA/protein)`,"/", 3)[,1],str_split_fixed(transvar$input,":",2)[,2])
  colnames(transvar_split) <- c("genomic_pos","SUB_MOD_RSD")
  transvar_split <- transvar_split[transvar_split$genomic_pos != ".",]
  transvar_split <- unique(transvar_split)
  
  #split the SUBSTRATE and SUB_MOD_RSD in the first column
  pho_rsd_split <- data.frame(str_split_fixed(pho_raw$antibody, "[_ ]", 2))
  colnames(pho_rsd_split) <- c("SUBSTRATE","antibody")
  
  pho_rsd_split$genomic_pos <- pho_raw$genomic_pos
  pho_rsd_split <- merge(pho_rsd_split,transvar_split, by=c("genomic_pos"), all.x = TRUE)
  
# initiate ----------------------------------------------------------------
  #initiating the lists to store info for each model
  list1 <- list()
  
# looping for model1 -----------------------------------------------------------------
  for (kinase in unique_kinase){
    # find protein expression level for the kinase
    pro_kinase <- pro_data[pro_raw$X == kinase,]
    
    if(nrow(pro_kinase) != 0){
      # find its substrate set
      k_k_s_table = k_s_table[k_s_table$GENE == kinase,]
      k_sub <- unique(k_k_s_table$SUB_GENE)
      
      for (substrate in k_sub){# for each substrate for one kinase
        # find its phosphosites-row numbers
        s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
        
        # go through all the phosphosites
        for (i in s_pho_table) {
          
          # find phosphorylation level
          sub_mod_rsd <- pho_rsd_split$SUB_MOD_RSD[i]
          pho_sub <- pho_data[i,]
          
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
  
  # combine table
  table <- table1
  
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

## write out tables
tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/RPPA_",protein,"_substrate_regression.txt", sep="")
write.table(table_2can, file=tn, quote=F, sep = '\t', row.names = FALSE)

## only plot result for model1
table_sig <- table_2can[table_2can$model=="pho_sub~pro_kin" & table_2can$FDR_pro_kin <= sig,]
table_sig_self <- table_sig[table_sig$self,]
table_sig_other <- table_sig[!table_sig$self,]

tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/RPPA_",protein,"_substrate_regression_sig_",sig,".txt", sep="")
write.table(table_sig, file=tn, quote=F, sep = '\t', row.names = FALSE)

# barplot for validation statistics ---------------------------------------
## only use data of BRCA and model1
table1 <- table_2can[table_2can$Cancer=="BRCA" & table_2can$model=="pho_sub~pro_kin",]
## different standard of validation for kinase and phosphotase
if ( protein == "kinase" ) {
  valid_type <- "pos_sig"
}
if ( protein == "phosphotase" ) {
  valid_type <- "neg_sig"
}

## initiate
x <- vector(mode = "numeric", length = length(unique_kinase) )
temp <- data.frame(matrix(rep(x,4), ncol=4, byrow=T))
valid_trans <- cbind(unique_kinase,temp)
colnames(valid_trans) <- c("kinase","pos_sig","pos_insig","neg_sig","neg_insig")
rownames(valid_trans) <- unique_kinase
valid_cis <- valid_trans

## make tables
### for trans pairs
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

## for count
ggplot() +
  geom_bar(data=table, aes(y = count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab(protein)+ylab("number of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_trans_gene_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

## for ratio
ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab(protein)+ylab("number of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_trans_gene_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

# scatterplot
p <- ggplot(valid_trans, aes(all_count, valid_ratio, label = rownames(valid_trans)))
p + geom_text(check_overlap = TRUE)
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_trans_all_count~valid_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

## for cis pairs
### for cis pairs
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

## for count
ggplot() +
  geom_bar(data=table, aes(y = count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab(protein)+ylab("number of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_cis_gene_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

## for ratio
ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab(protein)+ylab("number of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_cis_gene_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

# scatterplot
p <- ggplot(valid_cis, aes(all_count, valid_ratio, label = rownames(valid_cis)))
p + geom_text(check_overlap = TRUE)
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_cis_all_count~valid_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)


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
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_trans_residue_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = RESIDUE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("residue")+ylab("ratio of phosphosite residues with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_trans_residue_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

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
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_cis_residue_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = RESIDUE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("residue")+ylab("ratio of phosphosite residues with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,"/",protein,'_cis_residue_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

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

#table <- table_2can[table_2can$Cancer==c & table_2can$model==mode & table_2can$self==cis,]
table_2can$coef_pro_kin_filtered = remove_outliers(table_2can$coef_pro_kin)
table_2can_outlier_removed_m = table_2can[!is.na(table_2can$coef_pro_kin_filtered) & table_2can$model == mode,]

if ( protein == "kinase") {
  plot_fdr_scale <- 5
}
if ( protein == "phosphotase") {
  plot_fdr_scale <- 2
}
p = ggplot(table_2can_outlier_removed_m,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin)))
p = p + facet_grid(SELF~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(alpha=0.05)
p = p + geom_text(aes(label= ifelse(-log10(FDR_pro_kin)>plot_fdr_scale, pair, NA)),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/volcano_plots/',protein,'_substrate_volcano.pdf',sep ="")
ggsave(file=fn, height=6, width=8, useDingbats=FALSE)

# choose the parameters for bubble chart ---------------------
top <- 50 # choose top n rows for FDR_pro_kin for model1
c <- "BRCA" #choose in which cancer extract the top n rows
t0 <- table_sig_other  ## choose table_sig_other or table_sig_self

# extract top XX significant pairs in cis/trans group ---------------------
# sort by FDR_pro_kin
t1 <- t0[t0$Cancer==c,]
t1 <- t1[order(t1$FDR_pro_kin),] 

## corresponding results for other two models are extracted and ordered
rows <- c()
for(i in 1:top){
  r <- unlist(which(t0$pair==t1$pair[i]))
  rows <- c(rows,r)
}
table <- t0[rows,]
#table_2can_brca_top = table_2can[table_2can$pair %in% table$pair,]

# bubble chart plotting module ---------------------------------------------------------
lim = max(max(table$coef_pro_kin),min(table$coef_pro_kin))

## actual plotting
p = ggplot(table,aes(x=model, y=pair))# make this the original ethni
p = p + facet_grid(KINASE~Cancer, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_point(aes(fill=coef_pro_kin, size =-log10(FDR_pro_kin), color=ifelse(sig, "black",NA)),pch=21) 
p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
p = p + scale_colour_manual(values=c("black",NA))
p = p + theme_bw() #+ theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p
fn = paste(pd, 'kinase_substrate_fdr0.05.pdf',sep ="_")
ggsave(file=fn, height=10, width=5, useDingbats=FALSE)

