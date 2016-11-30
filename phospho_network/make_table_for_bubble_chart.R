### make_table_for_bubble.R ### 
# Yige Wu @ WashU 2016 Nov
# construct the table for bubbleplot
# load data from running validate_kinase_sub.R first

# construct table for model1
# columns: KINASE, SUBSTRATE, SUB_MOD_RSD, Pvalue, size, pair, self
KINASE <- sapply(list1, "[[", "KINASE")
SUBSTRATE <- sapply(list1, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list1, "[[", "SUB_MOD_RSD")
size <- sapply(list1, "[[", "Size")
table1 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,size)
table1$model <- "pho_sub~pro_kinase"
table1$P_pro_kin <- sapply(list1, "[[", "Pvalue")
table1$P_pro_sub <- NA
table1$P_pho_kin <- NA
table1$coef_pro_kin <- sapply(list1, "[[", "Coef_pro_kin")
table1$coef_pro_sub <- NA
table1$coef_pho_kin_g <- NA
table1$pair <- paste(table1$KINASE,":",table1$SUBSTRATE,":",table1$SUB_MOD_RSD,sep="")

# construct table for model2
KINASE <- sapply(list2, "[[", "KINASE")
SUBSTRATE <- sapply(list2, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list2, "[[", "SUB_MOD_RSD")
size <- sapply(list2, "[[", "Size")
table2 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,size)
table2$model <- "pho_sub~pro_kinase+pro_sub"
Pvalue <- t(sapply(list2, "[[", "Pvalue"))
table2$P_pro_kin <- Pvalue[,1]
table2$P_pro_sub <- Pvalue[,2]
table2$P_pho_kin <- NA
table2$coef_pro_kin <- sapply(list2, "[[", "Coef_pro_kin")
table2$coef_pro_sub <- sapply(list2, "[[", "Coef_pro_sub")
table2$coef_pho_kin_g <- NA
table2$pair <- paste(table2$KINASE,":",table2$SUBSTRATE,":",table2$SUB_MOD_RSD,sep="")

# construct table for model3
KINASE <- sapply(list3, "[[", "KINASE")
SUBSTRATE <- sapply(list3, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list3, "[[", "SUB_MOD_RSD")
size <- sapply(list3, "[[", "Size")
table3 <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,size)
table3$model <- "pho_sub~pro_kinase+pro_sub+pho_kin_grouped"
Pvalue <- t(sapply(list3, "[[", "Pvalue"))
table3$P_pro_kin <- Pvalue[,1]
table3$P_pro_sub <- Pvalue[,2]
table3$P_pho_kin <- Pvalue[,3]
table3$coef_pro_kin <- sapply(list3, "[[", "Coef_pro_kin")
table3$coef_pro_sub <- sapply(list3, "[[", "Coef_pro_sub")
table3$coef_pho_kin_g <- sapply(list3, "[[", "Coef_pho_kin_g")
table3$pair <- paste(table3$KINASE,":",table3$SUBSTRATE,":",table3$SUB_MOD_RSD,sep="")

# combine table1,2,3
table <- rbind(table1,table2,table3)

# mark cancer and self-regulation
table$Cancer <- cancer
table$self <- as.character(table$KINASE) == as.character(table$SUBSTRATE)

#choose one command
#table_BRCA <- table
table_OV <- table
