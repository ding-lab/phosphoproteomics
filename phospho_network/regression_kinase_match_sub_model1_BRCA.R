# Yige Wu @ WashU 2016 Nov
# look at correlations of kinase and downstream substrates phosphorylation status


# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)

# input regardless to cancer type-------------------------------------------------------------------
### read in the kinase/substrate table/ phosphorylation data ### 
K_S_f ="~/Box Sync/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt"
k_s_table = read.delim(K_S_f)

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

unique_kinase <- unique(k_s_table$GENE)

# choose significance level and outlier threshold and least sample number-------------------------
sig <- 0.05 # significance level
out_thres <- 1.5 #threshold for outlier
least_samples <- 5# least number of samples with complete data for each model

# looping cancer,run following sections all at once before plot module -----------------------------------------------------------
# input according to cancer type-------------------------------------------------------------------

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
    pro_kinase <- pro_data[pro_data$X == kinase,-colx]
    
    if(nrow(pro_kinase) != 0){
      # find its substrate set
      k_k_s_table = k_s_table[k_s_table$KINASE == kinase,]
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
  
  table <- table1
  
  # mark cancer and self-regulation
  table$Cancer <- cancer
  table$self <- as.character(table$KINASE) == as.character(table$SUBSTRATE)

# make the table for bubble chart -------------------------------------------------------
## combine table from BRCA and OV
table_2can <- table

## adjust p-values to FDR
for (mod in c("pho_sub~pro_kin", "pho_sub~pro_kin+pro_sub", "pho_sub~pro_kin+pro_sub+pho_kin","pho_sub~pho_kin")) {
  for (cancer in c("BRCA","OV")) {
    for(self in c(TRUE,FALSE)) {
      for(coln in "pro_kin") {#adjust pvalues for each variable
        row <- (table_2can$self==self) & (table_2can$model==mod) & (table_2can$Cancer==cancer)
        table_2can[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_2can[row,paste("P_",coln,sep = "")],method = "fdr")
      }
    }
  }
}


  
## only plot result for model1
table_sig <- table_2can[table_2can$model=="pho_sub~pro_kin" & table_2can$FDR_pro_kin <= sig,]
table_sig_self <- table_sig[table_sig$self,]
table_sig_other <- table_sig[!table_sig$self,]

# choose the parameters for bubble chart ---------------------
top <- 100 # choose top n rows for FDR_pro_kin for model1
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

## reshapte the table 
temp <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","size","model","P_pro_kin","P_pro_sub","P_pho_kin","pair","Cancer","self","order","pairs")
dm1 <- melt(table[,c(temp,paste("FDR_",name,sep = ""))], id=temp)
colnames(dm1) <- c(temp, "variable", "FDR")
dm1$var <- vector(mode = "numeric", length = nrow(dm1))
for (v in name) {
  dm1$var[dm1$variable==paste("FDR_",v,sep = "")] <- v
}


dm2 <- melt(table[,c(temp,paste("coef_",name,sep = ""))], id=temp)
colnames(dm2) <- c(temp, "variable2", "coefficients")
dm2$var <- vector(mode = "numeric", length = nrow(dm2))
for (v in name) {
  dm2$var[dm2$variable2==paste("coef_",v,sep = "")] <- v
}

molten <- merge(dm1, dm2)
molten <- molten[!is.na(molten$FDR),]
molten <- molten[molten$var!="pro_sub" & molten$model!="pho_sub~pro_kinase+pro_sub+pho_kin_grouped",]

# bubble chart plotting module ---------------------------------------------------------
## actual plotting
p = ggplot(table,aes(x=model, y=pair))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point(aes(colour=coef_pro_kin, size =-log10(FDR_pro_kin))) + theme_bw() #+ theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p

# volcano plotting module -------------------------------------------------
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

## choose the cancer, model and cis/trans
c <- "OV"
mode <- "pho_sub~pro_kin"
cis <- FALSE

table <- table_2can[table_2can$Cancer==c & table_2can$model==mode & table_2can$self==cis,]
x1 <- remove_outliers(table$coef_pro_kin)
p <- qplot(coef_pro_kin, -log10(FDR_pro_kin), data = table[table$coef_pro_kin==x1,])
p




# barplot for validation statistics ---------------------------------------
x <- vector(mode = "numeric", length = length(unique_kinase) )
temp <- data.frame(matrix(rep(x,4), ncol=4, byrow=T))
valid_trans <- cbind(unique_kinase,temp)
colnames(valid_trans) <- c("kinase","pos_sig","pos_insig","neg_sig","neg_insig")
rownames(valid_trans) <- unique_kinase
valid_cis <- valid_trans

#for( kinase in "EIF2AK1") {

## for trans pairs
for( kinase in unique_kinase) {
  temp <- table_2can[table_2can$KINASE==kinase,]
  valid_trans[kinase,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & !temp$self))
  valid_trans[kinase,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & !temp$self))
  valid_trans[kinase,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & !temp$self))
  valid_trans[kinase,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & !temp$self))
}
valid_trans$all_count <- rowSums(valid_trans[,c("pos_sig","pos_insig","neg_sig","neg_insig")])
valid_trans$valid_count <- valid_trans$pos_sig
valid_trans$valid_ratio <- valid_trans$pos_sig/valid_trans$all_count
valid_trans$positive <- (valid_trans$pos_sig+valid_trans$pos_insig)/valid_trans$all_count


table <- melt(valid_trans[!is.na(valid_trans$valid_ratio),],id=c("kinase","all_count","valid_count","valid_ratio","positive"))
colnames(table) <- c("kinase","all_count","valid_count","valid_ratio","positive","coef_FDR","count")
table$KINASE <- reorder(table$kinase, table$valid_count)

ggplot() +
  geom_bar(data=table, aes(y = count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("kinase")+ylab("number of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))

ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("kinase")+ylab("ratio of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))

