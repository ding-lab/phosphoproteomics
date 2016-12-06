c <- "BRCA" #choose the cancer you want to sort

t1 <- data.table(t0[t0$Cancer==c,], key="FDR_pro_kin1")# choose sort by which variable

## extract top XX significant pairs in self-regulated or other group and ordered
top <- 50 # choose top n rows for FDR_pro_kin1
rows <- c()
for(i in 1:top){
  r <- unlist(which(t0$pair==t1$pair[i]))
  rows <- c(rows,r)
}
table <- t0[rows,]
table$order <- nrow(table):1
table$pairs <- reorder(table$pair,table$order)
molten = melt(table, id = c("KINASE","SUBSTRATE","SUB_MOD_RSD","size","model","P_pro_kin","P_pro_sub","P_pho_kin","pair","Cancer","self","P_pro_kin_scale","order","pairs"))



sig1 <- vector(mode = "numeric", length = nrow(table1))
for (i in 1:nrow(table1)){
  substrate <- as.character(table1$SUBSTRATE[i])
  sub_mod_rsd <- as.character(table1$SUB_MOD_RSD[i])
  sig1[i] <- length(which(k_s_fit$FDR_pro_kin1[k_s_fit$SUB_GENE==substrate & k_s_fit$SUB_MOD_RSD==sub_mod_rsd] <=0.05))
}
table1$sig1 <- sig1

k1 <- c(); k2 <- c()
for( i in which(table1$SIZE_kin==2)) {
  substrate <- as.character(table1$SUBSTRATE[i])
  sub_mod_rsd <- as.character(table1$SUB_MOD_RSD[i])
  kinase1 <- list1[[i]]$KINASE[1]
  kinase2 <- list1[[i]]$KINASE[2]
  k1[length(k1)+1] <- k_s_fit$FDR_pro_kin1[k_s_fit$SUB_GENE==substrate & k_s_fit$SUB_MOD_RSD==sub_mod_rsd & k_s_fit$GENE==kinase1]
  k2[length(k2)+1] <- k_s_fit$FDR_pro_kin1[k_s_fit$SUB_GENE==substrate & k_s_fit$SUB_MOD_RSD==sub_mod_rsd & k_s_fit$GENE==kinase2]
}



k_s_fit[k_s_fit$FDR_pro_kin1,]



num_NA= nrow(m2)
test1 <- k_s_fit[k_s_fit$SUB_GENE=="AKT1",]
test2 <- k_s_fit[k_s_fit$SUB_GENE=="AKT1" & k_s_fit$SUB_MOD_RSD=="S473",]

test <- k_s_fit[!is.nan(unlist(k_s_fit$P_pro_kin1)),]

# construct table for model1
# columns: KINASE, SUBSTRATE, SUB_MOD_RSD, Pvalue, size, pair, self
SUBSTRATE <- sapply(list1, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list1, "[[", "SUB_MOD_RSD")
SIZE_dat <- sapply(list1, "[[", "SIZE_dat")
SIZE_kin <- sapply(list1, "[[", "SIZE_kin")
table1 <- data.frame(SUBSTRATE,SUB_MOD_RSD,SIZE_dat,SIZE_kin)

SUBSTRATE <- sapply(list2, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list2, "[[", "SUB_MOD_RSD")
SIZE_dat <- sapply(list2, "[[", "SIZE_dat")
SIZE_kin <- sapply(list2, "[[", "SIZE_kin")
SIG_pro_kin <- sapply(list2, "[[", "SIG_pro_kin")
SIG_pro_sub <- sapply(list2, "[[", "SIG_pro_sub")
P_pro_sub <- sapply(list2, "[[", "P_pro_sub")
table2 <- data.frame(SUBSTRATE,SUB_MOD_RSD,SIZE_dat,SIZE_kin,SIG_pro_kin,SIG_pro_sub,P_pro_sub)

SUBSTRATE <- sapply(list3, "[[", "SUBSTRATE")
SUB_MOD_RSD  <- sapply(list3, "[[", "SUB_MOD_RSD")
SIZE_dat <- sapply(list3, "[[", "SIZE_dat")
SIZE_kin <- sapply(list3, "[[", "SIZE_kin")
SIG_pro_kin <- sapply(list3, "[[", "SIG_pro_kin")
SIG_pro_sub <- sapply(list3, "[[", "SIG_pro_sub")
SIG_pho_kin <- sapply(list3, "[[", "SIG_pho_kin")
P_pro_sub <- sapply(list3, "[[", "P_pro_sub")
table3 <- data.frame(SUBSTRATE,SUB_MOD_RSD,SIZE_dat,SIZE_kin,SIG_pro_kin,SIG_pro_sub,SIG_pho_kin,P_pro_sub)



cat(paste("\n", sep = " "))
cat(paste("\n", sep = " "))
cat(paste("\n", sep = " "))
cat(paste("\n", sep = " "))
cat(paste("\n", sep = " "))
cat(paste("\n", sep = " "))
cat(paste("\n", sep = " "))
cat(paste("\n", sep = " "))
cat(paste("\n", sep = " "))
