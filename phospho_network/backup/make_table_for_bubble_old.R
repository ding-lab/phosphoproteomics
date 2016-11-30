### make_table_for_bubble.R ### 
# Yige Wu @ WashU 2016 Nov
# construct the table for bubbleplot
# load data from running validate_kinase_sub.R first

setwd("~/proteomics/pan3can_analysis/phospho_network")

# extract the sample size of complete data for model1
size <- c()
for(i in 1:length(list1)){
  data_t <- list1[[i]]$FIT$data
  size[1+length(size)] <- nrow(data_t[complete.cases(data_t),])
}

# construct table for model1
# columns: KINASE, SUBSTRATE, SUB_MOD_RSD, Pvalue, size, pair, self
KINASE <- sapply(list1, "[[", "KINASE")
SUBSTRATE <- sapply(list1, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list1, "[[", "SUB_MOD_RSD")
Pvalue <- sapply(list1, "[[", "Pvalue")
table1 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,Pvalue,size)
table1$pair <- paste(table1$KINASE,":",table1$SUBSTRATE,":",table1$SUB_MOD_RSD,sep="")
table1$model <- "pho_sub~pro_kinase"

# extract the sample size of complete data for model2
size <- c()
for(i in 1:length(list2)){
  data_t <- list2[[i]]$FIT$data
  size[1+length(size)] <- nrow(data_t[complete.cases(data_t),])
}

# construct table for model2
KINASE <- sapply(list2, "[[", "KINASE")
SUBSTRATE <- sapply(list2, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list2, "[[", "SUB_MOD_RSD")
Pvalue <- sapply(list2, "[[", "Pvalue")
table2 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,Pvalue,size)
table2$pair <- paste(table2$KINASE,":",table2$SUBSTRATE,":",table2$SUB_MOD_RSD,sep="")
table2$model <- "pho_sub~pro_kinase+pro_sub"

# combine table1$2
table <- rbind(table1,table2)

# mark cancer and self-regulation
table1$Cancer <- cancer
table$self <- as.character(table$KINASE) == as.character(table$SUBSTRATE)

# combine the table for two cancers
table_2can <- rbind(table_2can,table)
