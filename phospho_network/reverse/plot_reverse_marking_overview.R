# choose the kinases to plot --------------------------------------------------------
outlier_num <- data.frame(is.pro_outlier$kinase,rowSums(is.pro_outlier[,sample_names]))
colnames(outlier_num) <- c("kinase","num_pro_outlier")
outlier_num$num_pro_outlier[is.na(outlier_num$num_pro_outlier)] <- 0

missense_num <- data.frame(mut$kinase,rowSums(mut[,sample_names]))
colnames(missense_num) <- c("kinase","missense_num")
missense_num$missense_num[is.na(missense_num$missense_num)] <- 0

outlier_num <- merge(outlier_num,missense_num, all.y = TRUE)
outlier_num$num_pro_outlier[is.na(outlier_num$num_pro_outlier)] <- 0

range <- outlier_num$kinase[outlier_num$missense_num==1]


# Plot kinase score and outlier_status---------------------------------------
row <- c()
order <- c()
for (kinase in range) {
  temp <- which(overlap$kinase==kinase & !is.na(overlap$score))
  row <- c(row,temp)
  if (length(temp) > 0) {
    order <- c(order, rep(outlier_num$missense_num[outlier_num$kinase==kinase], length(temp)))
  }
}
table <- overlap[row,]
table$kinases <- reorder(table$kinase,order)
pd <- position_dodge(.65)
ggplot(data = table,aes(x= kinases, y = score, colour = missense.pro_outlier, shape = missense.pro_outlier)) +    
  geom_point(position = pd, size = 4) +
  scale_colour_manual(name = "missense.pro_outlier",
                      labels = c("wt,pro_normal", "wt,pro_outlier", "missense,pro_normal", "missense,pro_outlier"),
                      values = c("blue", "red", "blue", "red")) +   
  scale_shape_manual(name = "missense.pro_outlier",
                     labels = c("wt,pro_normal", "wt,pro_outlier", "missense,pro_normal", "missense,pro_outlier"),
                     values = c(19, 19, 17, 17))+
  coord_flip()

dev.off()

# calculate substrate-phospho-outlier sample number for each kinase -------------------------
## construct a table, 1==outlier, 0==not outlier, NA==no data
x <- vector(mode = "logical", length = length(unique_kinase))
temp <- data.frame(matrix(rep(x,ncol(pho_data)-1), ncol=ncol(pho_data)-1, byrow=T))
outlierk <- cbind(unique_kinase,temp)
sample_names <- colnames(pho_data[-colx])
colnames(outlierk) <- c("KINASE", sample_names)

for (i in 1:length(unique_kinase)) {
  IQR = quantile(k_p_table[i,-1], probs=0.75, na.rm=T) - quantile(k_p_table[i,-1], probs=0.25, na.rm=T) 
  outlierk[i,-1] = (k_p_table[i,-1] >= quantile(k_p_table[i,-1], probs=0.75, na.rm=T) + out_thres*IQR)
}
outlierk$sum <- rowSums(outlierk[,-1], na.rm = TRUE)
row.names(outlierk) <- outlierk$KINASE

# barplot for number of phospho-outlier sample for each kinase -------------------------------------------------------------
outlierk_order <- outlierk[outlierk$sum>0,]
outlierk_order$kinase <- reorder(outlierk_order$KINASE,outlierk_order$sum)

ggplot(data=outlierk_order, aes(x=kinase, y=sum)) +
  geom_bar(stat="identity") +
  coord_flip()

# calculate outlier kinase number for each sample -------------------------
kinase_sum <- t(vector(mode = "numeric", length = length(sample_names))+NaN)
kinase_sum <- colSums(outlierk[,sample_names], na.rm = TRUE)

outliers <- data.frame(kinase_sum,sample_names)
outliers$sample <- reorder(outliers$sample_names,outliers$kinase_sum)

# barplot for number of outlier kinases for each sample -------------------
ggplot(data=outliers, aes(x=sample, y=kinase_sum)) +
  geom_bar(stat="identity") +
  coord_flip()

# Calculate per kinase mutation distribution in phospho-outlier/no --------
## cols: kinase, mutation type, outlier_type, sample number
## rows: No.kinase(346)*No.mutation_type*outlier/normal(2)
KINASE <- c();  outlier_type <- c(); mutation_type <- c();sample_num <- c();noutlier <- c()
len <- 0
for (kinase in outlierk_order$kinase) {
  #for (kinase in "EIF2AK1"){
  mut_table <- levels(as.factor(as.vector(t(somatic[somatic$X ==kinase, -colx ]))))
  nmut <- length(mut_table)
  order <- outlierk$sum[outlierk$KINASE==kinase]
  for (out in c(0,1)) {
    if (nmut == 0){
      len <- len+1
      KINASE[len] <- kinase
      outlier_type[len] <- out
      mutation_type[len] <- "unknown"
      sample_num[len] <- length(which(outlierk_order[outlierk_order$kinase==kinase,sample_names]==out))
      noutlier[len] <- order
    } else {
      for (mut in mut_table) {
        len <- len+1
        KINASE[len] <- kinase
        outlier_type[len] <- out
        mutation_type[len] <- mut
        sample_num[len] <- length(which(somatic[somatic$X==kinase,sample_names]==mut & outlierk_order[outlierk_order$kinase==kinase,sample_names]==out))
        noutlier[len] <- order
      }
    }
  }
}
mut_pho_table <- data.frame(KINASE,outlier_type,mutation_type,sample_num,noutlier)
table <- mut_pho_table[mut_pho_table$outlier_type==1,]
table$kinase <- reorder(table$KINASE,table$noutlier)

# Barplot for somatic mutation & phosphor-outlier distribution ------------
ggplot() +
  geom_bar(data=table, aes(y = sample_num, x = kinase, fill = mutation_type), stat="identity",
           position='stack') +
  theme_bw() + 
  #facet_grid( ~ outlier_type)+ 
  xlab("kinase")+ylab("number of phospho-outlier samples")+
  coord_flip()

# Calculate per kinase protein outlier distribution in phospho-outlier/no --------
KINASE <- c();  outlier_type <- c(); pro_type <- c();sample_num <- c();noutlier <- c()
len <- 0
for (kinase in outlierk_order$kinase) {
  #for (kinase in "JAK3"){
  temp <- length(which(pro_outlier$X==kinase))
  order <- outlierk$sum[outlierk$KINASE==kinase]
  for (out in c(0,1)) {
    if (temp == 0){
      len <- len+1
      KINASE[len] <- kinase
      outlier_type[len] <- out
      noutlier[len] <- order
      pro_type[len] <- "normal_pro_gene"
      sample_num[len] <- length(which(outlierk_order[outlierk_order$kinase==kinase,sample_names]==out))
    } else {
      len <- len+1
      KINASE[len] <- kinase
      outlier_type[len] <- out
      pro_type[len] <- "outlier_pro_gene_pro_outlier"
      sample_num[len] <- length(which(pro_outlier[pro_outlier$X==kinase,sample_names] >= 1 & outlierk_order[outlierk_order$kinase==kinase,sample_names]==out))
      noutlier[len] <- order
      
      len <- len+1
      KINASE[len] <- kinase
      outlier_type[len] <- out
      pro_type[len] <- "outlier_pro_gene_pro_normal"
      sample_num[len] <- length(which(pro_outlier[pro_outlier$X==kinase,sample_names] < 1 & outlierk_order[outlierk_order$kinase==kinase,sample_names]==out))
      noutlier[len] <- order
      
      len <- len+1
      KINASE[len] <- kinase
      outlier_type[len] <- out
      pro_type[len] <- "outlier_pro_gene_pro_NA"
      sample_num[len] <- length(which(is.na(pro_outlier[pro_outlier$X==kinase,sample_names]) & outlierk_order[outlierk_order$kinase==kinase,sample_names]==out))
      noutlier[len] <- order
    }
  }
}
pro_pho_table <- data.frame(KINASE,outlier_type,pro_type,sample_num,noutlier)
table <- pro_pho_table[pro_pho_table$outlier_type==1,]
table$kinase <- reorder(table$KINASE,table$noutlier)

# Barplot for protein & phosphor-outlier distribution ------------
ggplot() +
  geom_bar(data=table, aes(y = sample_num, x = kinase, fill = pro_type), stat="identity",
           position='stack') +
  theme_bw() + 
  #facet_grid( ~ outlier_type)+ 
  xlab("kinase")+ylab("number of phospho-outlier samples")+
  coord_flip()

# Find kinase activation coocurrence --------------------------------------
## construct kinase*kinase table, each cell is the P-value for fisher-test
nkinase <- length(unique_kinase)
x <- vector(mode = "numeric", length = nkinase)+NaN
P_co_kin <- data.frame(matrix(rep(x,nkinase) ,ncol=nkinase, byrow=T))
colnames(P_co_kin) <- unique_kinase
rownames(P_co_kin) <- unique_kinase

## initiate a 2*2 table for fisher test
fisher_table <- matrix(data = NA, nrow = 2, ncol = 2)
unique_kinase <- as.vector(unique_kinase)
## 2 loops around each kinase
for (i1 in 1:(nkinase-1)) {
  #for (k1 in "EIF2AK1") {
  k1 <- unique_kinase[i1]
  out1 <- outlierk[k1,sample_names]
  for (i2 in (i1+1):nkinase) {
    k2 <- unique_kinase[i2]
    #for (k2 in "PRKCD") {
    out2 <- outlierk[k2,sample_names]
    for (on1 in c(0,1)) {#### kinase1 activation status: 0-off, 1-on
      for (on2 in c(0,1)) {#### kinase2 activation status: 0-off, 1-on
        ### collect the numbers for fisher test
        fisher_table[on1+1,on2+1] <- length(which(out1==on1 & out2==on2))
      }
    }
    #print(fisher_table)
    ### fisher test
    temp <- fisher.test(fisher_table)
    P_co_kin[k1,k2] <- temp$p.value
  }
}

## adjust p-values to FDR
### reshapte P-value table to be a long vector
temp <- as.vector(t(P_co_kin))

### adjust to FDR
p2fdr.na.omit <- function(x) {
  y <- x
  pos <- which(!is.na(y))
  pvalues <- y[pos]
  fdrs <- p.adjust(pvalues, method = "fdr")
  y[pos] <- fdrs
  return(y)
}#actually there's no need to use this function, p.adjust will ignore NAs
FDR_co_kin <- p2fdr.na.omit(temp)

### reshapte to make a FDR table
FDR_co_kin <- matrix(data = FDR_co_kin, nrow = nkinase, ncol = nkinase, byrow = TRUE)
colnames(FDR_co_kin) <- unique_kinase
rownames(FDR_co_kin) <- unique_kinase

## histogram of number of cooccurence
nco <- colSums(FDR_co_kin <= sig, na.rm = TRUE)
co_sig <- which(FDR_co_kin <= sig, arr.ind = TRUE)
for (i in 1:nrow(co_sig)) {
  print(paste(unique_kinase[co_sig[i,1]],unique_kinase[co_sig[i,2]],FDR_co_kin[co_sig[i,1],co_sig[i,2]], sep = ":"))
}

## record fisher test result for those with significant FDRs
fisher_sig <- list()
for (i in 1:nrow(co_sig)) {
  k1 <- unique_kinase[co_sig[i,1]]
  out1 <- outlierk[k1,sample_names]
  k2 <- unique_kinase[co_sig[i,2]]
  out2 <- outlierk[k2,sample_names]
  for (on1 in c(0,1)) {#### kinase1 activation status: 0-off, 1-on
    for (on2 in c(0,1)) {#### kinase2 activation status: 0-off, 1-on
      ### collect the numbers for fisher test
      fisher_table[on1+1,on2+1] <- length(which(out1==on1 & out2==on2))
    }
  }
  temp <- fisher.test(fisher_table)
  fisher_sig[[i]] <- list(kinase1 = k1, kinase2 = k2, FDR = FDR_co_kin[co_sig[i,1],co_sig[i,2]], fisher_test = temp, estimate = temp$estimate)  
}

barplot(nco[nco>0])

# Check the scale of tests to be performed --------------------------------
## check how big will be the table and how many test will be performed
unique_kin_pho <- intersect(unique_kinase, pho_rsd_split$SUBSTRATE)

### if only test substrates phosphosites in the k_s_table
n_sub_phos <- c()
n_kin_phos <- c()

for (kinase in unique_kin_pho) {
  
  k_k_s_table <- k_s_table[k_s_table$GENE==kinase,]
  
  nrsdk <- c()
  ns <- c()
  nrsds <- c()
  
  for (i in 1:nrow(k_k_s_table)) {
    substrate <- as.character(k_k_s_table$SUB_GENE[i])
    sub_mod_rsd <- as.character(k_k_s_table$SUB_MOD_RSD[i])
    rows <- pho_rsd_split$SUBSTRATE==substrate & pho_rsd_split$SUB_MOD_RSD==sub_mod_rsd
    if ( length(which(rows)) >0 ) {
      ns <- c(ns,substrate)
      nrsds <- c(nrsds,sub_mod_rsd)
    }
  }
  
  n_sub_phos <- rbind(n_sub_phos, cbind(ns, nrsds))
  n_kin_phos <- rbind(n_kin_phos, cbind(kinase, nrsdk))
}

### if test all substrates phosphosites in the pho_data
# n_sub_phos <- 0
# n_kin_phos <- 0
# n_sub_kin_phos <- 0
# for (kinase in unique_kin_pho) {
#   k_k_s_table <- as.vector(unique(k_s_table$SUB_GENE[k_s_table$GENE==kinase]))
#   nk <- 0
#   ns <- 0
#   for (substrate in k_k_s_table) {
#       ns <- ns + length(which(pho_rsd_split$SUBSTRATE==substrate))
#   }
#   
#   nk <- length(which(pho_rsd_split$SUBSTRATE==kinase))
#   n_sub_phos <- n_sub_phos + ns
#   n_kin_phos <- n_kin_phos + nk
#   n_sub_kin_phos <- n_sub_kin_phos + ns*nk
# }

# Test for phosphorylation level correlation between kinase phosphosites and substrate phosphosites(only match substrate phosphosites in k_s_table) --------
## initiate a phosphosites(kinase)*phosphosites(substrate) table

## record all the substrate phosphosite in k_s_table that has phosphorylation data
phos_all <- paste(pho_rsd_split$SUBSTRATE,pho_rsd_split$SUB_MOD_RSD, sep = ":")
phos_sub <- paste(k_s_table$SUB_GENE,k_s_table$SUB_MOD_RSD,sep = ":")
rows<- match(phos_sub, phos_all) ### a vector of row number in pho_data for each substrate phosphosite in k_s_table
phos_sub <- intersect(phos_sub, phos_all)
phos_sub <- data.frame(str_split_fixed(phos_sub, ":", 2))
colnames(phos_sub) <- c("substrate","rsd")

## record all the kinase phosphosite in k_s_table that has phosphorylation data
phos_kin <- match(pho_rsd_split$SUBSTRATE,unique_kinase)
rowk <- which(!is.na(phos_kin)) ### a vector of row number in pho_data for each kinase phosphosite
phos_kin <- pho_rsd_split[rowk,c("SUBSTRATE","SUB_MOD_RSD")]
colnames(phos_kin) <- c("kinase","rsd")

## 3 tables for p-values, coefficients and FDR, rows for kinase, cols for substrates
x <- vector(mode = "numeric", length = nrow(phos_kin))+NaN
P_corr_pho <- data.frame(matrix(rep(x,nrow(phos_sub)) ,ncol=nrow(phos_sub), byrow=T))
coef_corr_pho <- P_corr_pho
FDR_corr_pho <- P_corr_pho

## preprocess the phosphorylation level for phosphosites?
list_rsdk <- list()
for (kinase in as.character(phos_kin$kinase)) {
  list_rsdk[[length(list_rsdk)+1]] <- list(KINASE=kinase, rsdk=as.character(phos_kin$rsd[phos_kin$kinase==kinase]))
}
kins_pho <- as.character(phos_kin$kinase)
subs <- as.character(k_s_table$SUB_GENE)
rsdss <- as.character(k_s_table$SUB_MOD_RSD)

start.time <- Sys.time()
for (j in 1:length(kins_pho)) {
  kinase <- kins_pho[j]
  
  for (rsdk in list_rsdk[[j]]$rsdk) { ## loop around all the substrate phosphosites
    k_pho <- pho_data[pho_rsd_split$SUBSTRATE==kinase & pho_rsd_split$SUB_MOD_RSD==rsdk,-colx]
    
    row <- which(phos_kin$kinase==kinase & phos_kin$rsd==rsdk)
    for (i in which(k_s_table$KINASE==kinase & !is.na(rows))) { ## loop around all the substrate phosphosites
      
      s_pho <- pho_data[rows[i],-colx]
      table <- data.frame(t(rbind(k_pho, s_pho)))
      if ( length(which(complete.cases(table))) >= 5) {
        colnames(table) <- c("k_pho","s_pho")
        fit <- glm(s_pho~k_pho, data = table )
        
        substrate <- subs[i]
        rsds <- rsdss[i]
        col <- which(phos_sub$substrate==substrate & phos_sub$rsd==rsds)
        
        P_corr_pho[row,col] <- signif(coef(summary(fit))[-1,4],digits = 3)
        coef_corr_pho[row,col] <- signif(fit$coefficients[-1],digits = 3)
      }
    }
  }
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
# Time difference of 19.1747 mins

## adjust to FDR
temp <- as.vector(t(P_corr_pho))
FDR_corr_pho <- p.adjust(temp, method = "fdr")

## reshapte to make a FDR table
FDR_corr_pho <- matrix(data = FDR_corr_pho, nrow = nrow(phos_kin), ncol = nrow(phos_sub), byrow = TRUE)
corr_sig <- which(FDR_corr_pho <= sig, arr.ind = TRUE)
FDR <- FDR_corr_pho[corr_sig]
corr_sig_phos <- cbind(phos_kin[corr_sig[,1],], phos_sub[corr_sig[,2],], FDR )
colnames(corr_sig_phos) <- c("kinase","rsdk","substrate","rsds","FDR")
same <- which(corr_sig_phos$kinase==as.character(corr_sig_phos$substrate) & corr_sig_phos$rsdk==as.character(corr_sig_phos$rsds))
corr_sig_phos <- corr_sig_phos[setdiff(1:nrow(corr_sig_phos),same),]  
write.table(corr_sig_phos, file = "phosphosites-correlation-fdr.txt", row.names = F)


# scale by kinase/row -----------------------------------------------------
k_p_scale <- as.matrix(k_p_table[,-1])
k_p_scale <- data.frame(t(scale(t(k_p_scale))))
k_p_scale <- cbind(k_p_table$KINASE,k_p_scale)
