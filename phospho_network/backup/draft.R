erbb2 <- k_p_her[k_p_her$kinase=="ERBB2",]
rownames(erbb2) <- erbb2$sample

IQR = quantile(erbb2$ave_sub_phos, probs=0.75, na.rm=T) - quantile(erbb2$ave_sub_phos, probs=0.25, na.rm=T) 
erbb2$phos_outlier = (erbb2$ave_sub_phos >= quantile(erbb2$ave_sub_phos, probs=0.75, na.rm=T) + 1.5*IQR)

erbb2$pro_level <- vector(mode = "numeric", length = nrow(erbb2))+NaN
erbb2$pro_outlier <- vector(mode = "logical", length = nrow(erbb2))+NaN
for (sam in sample_names) {
  erbb2[sam,"pro_level"] <- pro_outlier["ERBB2",sam]
}
IQR = quantile(erbb2$pro_level, probs=0.75, na.rm=T) - quantile(erbb2$pro_level, probs=0.25, na.rm=T) 
erbb2$pro_outlier = (erbb2$pro_level >= quantile(erbb2$pro_level, probs=0.75, na.rm=T) + 1.5*IQR)




# scatterplot regression plotting module ----------------------------------
lm_eqn = function(df){
  m = glm(pho_sub_norm ~ pro_kinase, data = df);
  eq <- substitute(italic(pho_sub_norm) == a + b %.% italic(pro_kinase),
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2)));
  as.character(as.expression(eq));
}
table <- table_sig[table_sig$Cancer=="BRCA",]

# choose kinase
kinase = "MAPK3"

# extract the data for this kinase in model1 with FDR<0.05
data_k <- table[table$KINASE == kinase,]

eq <- ddply(data_k,.(SUBSTRATE,SUB_MOD_RSD),lm_eqn)

p <- ggplot(data = data_k, aes(x = pro_kinase, y = pho_sub_norm)) +
  geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ggtitle(paste("kinase = ",kinase,sep=""))+
  ylim(c(0,1))

p1 = p + facet_grid(SUBSTRATE+SUB_MOD_RSD~.)
p2 = p1 + geom_text(data=eq,aes(x = 0.5*max(data_k$pro_kinase), y = 0.9,label=V1, size=5), parse = TRUE, inherit.aes=FALSE)
p2
#BRCA_cor_pro_kinase_phosphosite_

## one kinase, one substrate

kinase = "CDK2";
sub = "APC";
rows <- which(table$KINASE == kinase & table$SUBSTRATE == sub & table$P_pro_kin <= 0.05 & table$model == "pho_sub~pro_kinase")

data_k <- data.frame()
for(i in rows ){
  # extract the data for each pair
  data <- list1[[i]]$DATA
  data$SUBSTRATE <- list1[[i]]$SUBSTRATE
  data$SUB_MOD_RSD <- list1[[i]]$SUB_MOD_RSD
  data_k <- rbind(data_k,data)
}
data_k <- data_k[complete.cases(data_k),]

eq <- ddply(data_k,.(SUBSTRATE,SUB_MOD_RSD),lm_eqn)

p <- ggplot(data = data_k, aes(x = pro_kinase, y = pho_sub_norm)) +
  geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ggtitle(paste("kinase = ",kinase,sep=""))+
  ylim(c(0,1))

p1 = p + facet_grid(SUBSTRATE+SUB_MOD_RSD~.)
p2 = p1 + geom_text(data=eq,aes(x = 0.5*max(data_k$pro_kinase), y = 0.9,label=V1, size=5), parse = TRUE, inherit.aes=FALSE)
p2





for (i1 in 1:nkinase) {
  for (i2 in 1:i1) {
    P_co_kin[i1,i2] <- NaN
  }
}
test <- as.vector(t(P_co_kin_all))


library(dplyr)
somatic %>% 
  sapply(levels)



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
