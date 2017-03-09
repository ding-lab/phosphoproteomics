# barplot for validation statistics ---------------------------------------
## only use data of BRCA and model1
baseD <-"/Users/yigewu/Box Sync/"
load(paste(baseD,"pan3can_shared_data/analysis_results/regression/Rdata/kinase.RData",sep = ""))
table1 <- table_2can[table_2can$Cancer=="BRCA" & table_2can$model=="pho_sub~pro_kin",]

## RTK gene list
RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
rtk <- as.vector(t(RTK_file))
protein <- "kinase"

# common kinase of RTKs ---------------------------------------------------
row_regression <- c()
for (gene in rtk) {
  temp <- which(table1$SUBSTRATE==gene)
  row_regression <- c(row_regression,temp)
}
table(c(as.vector(table1$KINASE))[row_regression])
# ABL1   AKT1   CDK1   CDK5  CHEK1    CSK   EGFR  ERBB2  ERBB4  FGFR4   GRK6  GSK3B  IGF1R    KIT  MAPK1  MAPK3 
# 11     13     13      6      2      6     13     22      3      1      1      5      1      2     13     13 
# PDGFRA PDGFRB   PDK1   PKN1  PRKCA   PTK2    RET    SRC    SYK 
# 4     10      5     13     44      2      2     42      3 
#for(kinase in unique(c(as.vector(table1$KINASE))[row_regression])) {
for(kinase in "PRKCA"){
  row_rtk <- c()
  for (gene in rtk) {
    temp <- which(table1$SUBSTRATE==gene & table1$KINASE==kinase)
    row_rtk <- c(row_rtk,temp)
  }
  if (length(which(table1$FDR_pro_kin[row_rtk] <= sig)) > 0) {
    print(kinase)
  }
}

# [1] "ERBB2"
# [1] "PRKCA" EGFR T693
# [1] "SRC" 
# [1] "EGFR"
# [1] "IGF1R"
# [1] "KIT"
# [1] "PDGFRB"
# [1] "SYK"
# [1] "ERBB4"


# barplot -----------------------------------------------------------------
test <- table1[row_rtk,]

valid_type <- "pos_sig"

## initiate
x <- vector(mode = "numeric", length = length(rtk) )
temp <- data.frame(matrix(rep(x,4), ncol=4, byrow=T))
valid_trans <- cbind(rtk,temp)
colnames(valid_trans) <- c("kinase","pos_sig","pos_insig","neg_sig","neg_insig")
rownames(valid_trans) <- rtk
valid_cis <- valid_trans

## for cis pairs
### for cis pairs
table_rtk <- c()
for( kinase in rtk) {
  temp <- table1[table1$KINASE==kinase,]
  if ( nrow(temp) > 0 ) {
    table_rtk <- rbind(table_rtk,temp)
    valid_cis[kinase,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & temp$self))
    valid_cis[kinase,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & temp$self))
    valid_cis[kinase,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & temp$self))
    valid_cis[kinase,"neg_insig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & temp$self))
  }
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
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/RTKs/RTKs_cis_gene_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

# scatterplot
p <- ggplot(valid_cis, aes(all_count, valid_ratio, label = rownames(valid_cis)))
p + geom_text(check_overlap = TRUE)
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/RTKs/RTKs_cis_all_count~valid_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

table_rtk_trans <- c()
for( kinase in rtk) {
  temp <- table1[table1$KINASE==kinase,]
  if ( nrow(temp) > 0 ) {
    table_rtk_trans <- rbind(table_rtk_trans,temp)
    valid_trans[kinase,"pos_sig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin <= sig & !temp$self))
    valid_trans[kinase,"pos_insig"] <- length(which(temp$coef_pro_kin>0 & temp$FDR_pro_kin > sig & !temp$self))
    valid_trans[kinase,"neg_sig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin <= sig & !temp$self))
    valid_trans[kinase,"neg_insig"] <- length(which(temp$coef_pro_kin<0 & temp$FDR_pro_kin > sig & !temp$self))
  }
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
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/RTKs/RTKs_trans_gene_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

## for ratio
ggplot() +
  geom_bar(data=table, aes(y = count/all_count, x = KINASE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab(protein)+ylab("ratio of phosphosites with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/RTKs/RTKs_trans_gene_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

# scatterplot
p <- ggplot(valid_trans, aes(all_count, valid_ratio, label = rownames(valid_trans)))
p + geom_text(check_overlap = TRUE)
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/RTKs/RTKs_trans_all_count~valid_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

## residue
## initiate
x <- vector(mode = "numeric", length = length(3) )
temp <- data.frame(matrix(rep(x,4), ncol=4, byrow=T))
valid_rsd_trans <- cbind(c("S","T","Y"),temp)
colnames(valid_rsd_trans) <- c("residue","pos_sig","pos_insig","neg_sig","neg_insig")
rownames(valid_rsd_trans) <- c("S","T","Y")
valid_rsd_cis <- valid_rsd_trans

## trans-residue
for(rsd in c("S","T","Y") ) {
  temp <- table_rtk_trans[grepl(rsd,table_rtk_trans$SUB_MOD_RSD),]
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
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/RTKs/RTKs_trans_residue_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

ggplot() +
  geom_bar(data=table, aes(y = count/all_count , x = RESIDUE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("residue")+ylab("ratio of phosphosite residues with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/RTKs/RTKs_trans_residue_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)


## cis-residue
for(rsd in c("S","T","Y") ) {
  temp <- table_rtk[grepl(rsd,table_rtk$SUB_MOD_RSD),]
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
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/RTKs/RTKs_cis_residue_count.png',sep ="")
ggsave(file=fn, height=6, width=8)

ggplot() +
  geom_bar(data=table, aes(y = count/all_count , x = RESIDUE, fill = coef_FDR), stat="identity",
           position='stack') +
  theme_bw() + 
  xlab("residue")+ylab("ratio of phosphosite residues with different regression result")+
  coord_flip()+ theme(axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
fn = paste(baseD,'pan3can_shared_data/analysis_results/regression/overview/',protein,'/RTKs/RTKs_cis_residue_ratio.png',sep ="")
ggsave(file=fn, height=6, width=8)

# make table for cytoscape ------------------------------------------------
table <- table_rtk
self <- TRUE
sig <- 0.1
table <- table[table$FDR_pro_kin <= sig & table$self==self,]
table$scale_FDR <- -log10(table$FDR_pro_kin)

tn = paste(baseD,"pan3can_shared_data/analysis_results/regression/table/RTKs_self_",self,"_fdr_",sig,".txt", sep="")
write.table(table, file=tn, quote=F, sep = '\t', row.names = FALSE)

