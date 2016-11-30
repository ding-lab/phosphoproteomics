### make_table_for_bubble.R ### 
# Yige Wu @ WashU 2016 Nov
# construct the table for bubbleplot
# load data from running validate_kinase_sub.R first

setwd("~/proteomics/pan3can_analysis/phospho_network")

# extract the sample size of complete data for model1
size <- c()
for(i in 1:length(list1)){
  data_t <- list1[[i]]$DATA
  size[1+length(size)] <- nrow(data_t[complete.cases(data_t),])
}

# construct table for model1
# columns: KINASE, SUBSTRATE, SUB_MOD_RSD, Pvalue, size, pair, self
KINASE <- sapply(list1, "[[", "KINASE")
SUBSTRATE <- sapply(list1, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list1, "[[", "SUB_MOD_RSD")
P_pro_kin <- sapply(list1, "[[", "Pvalue")
table1 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,P_pro_kin,size)
table1$pair <- paste(table1$KINASE,":",table1$SUBSTRATE,":",table1$SUB_MOD_RSD,sep="")
table1$model <- "pho_sub~pro_kinase"
table1$P_pro_sub <- NA
table1$P_pho_kin <- NA


# extract the sample size of complete data for model2
size <- c()
for(i in 1:length(list2)){
  data_t <- list2[[i]]$DATA
  size[1+length(size)] <- nrow(data_t[complete.cases(data_t),])
}

# construct table for model2
KINASE <- sapply(list2, "[[", "KINASE")
SUBSTRATE <- sapply(list2, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list2, "[[", "SUB_MOD_RSD")
Pvalue <- t(sapply(list2, "[[", "Pvalue")); P_pro_kin <- Pvalue[,1]
table2 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,P_pro_kin,size)
table2$pair <- paste(table2$KINASE,":",table2$SUBSTRATE,":",table2$SUB_MOD_RSD,sep="")
table2$model <- "pho_sub~pro_kinase+pro_sub"
table2$P_pro_sub <- Pvalue[,2]
table2$P_pho_kin <- NA

# extract the sample size of complete data for model3
size <- c()
for(i in 1:length(list3)){
  data_t <- list3[[i]]$DATA
  size[1+length(size)] <- nrow(data_t[complete.cases(data_t),])
}

# construct table for model3
KINASE <- sapply(list3, "[[", "KINASE")
SUBSTRATE <- sapply(list3, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list3, "[[", "SUB_MOD_RSD")
Pvalue <- t(sapply(list3, "[[", "Pvalue")); P_pro_kin <- Pvalue[,1]
table3 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,P_pro_kin,size)
table3$pair <- paste(table3$KINASE,":",table3$SUBSTRATE,":",table3$SUB_MOD_RSD,sep="")
table3$model <- "pho_sub~pro_kinase+pro_sub+pho_kin_grouped"
table3$P_pro_sub <- Pvalue[,2]
table3$P_pho_kin <- Pvalue[,3]

# combine table1,2,3
table <- rbind(table1,table2,table3)

# mark cancer and self-regulation
table$Cancer <- cancer
table$self <- as.character(table$KINASE) == as.character(table$SUBSTRATE)

#choose one command
table_BRCA <- table
# tble_OV <- table
