### plot_bubble_chart.R ### 
# Yige Wu @ WashU 2016 Nov
# plot correlation p-value and sample size between substrate phosphorylation level and kinase protein expression level, etc
# use make_table_for_bubble.R to prepapre the table

setwd("~/proteomics/pan3can_analysis/phospho_network")
library(ggplot2)
library(ggtern)
require(data.table)
library(reshape)

load("~/proteomics/pan3can_analysis/phospho_network/tables_for_bubble_chart.RData")

# merge table for BRCA and OV, keep significant ones, divide self-regulated and others
table_2can <- rbind(table_BRCA,table_OV)
table_2can$P_pro_kin_scale <- round(-log10(table_2can$P_pro_kin),digit=0)

# adjust pvalues to fdrs according to multiple tests and whether kinase==substrate
x <- vector(mode = "numeric", length = nrow(table_2can))
temp <- data.frame(matrix(rep(x,3), ncol=3, byrow=T))
name <- c("pro_kin","pro_sub","pho_kin")
colnames(temp) <- c(paste("FDR_",name,sep = ""))
table_2can <- cbind(table_2can,temp)

for (mod in c("pho_sub~pro_kinase", "pho_sub~pro_kinase+pro_sub", "pho_sub~pro_kinase+pro_sub+pho_kin_grouped")) {
  for (cancer in c("BRCA","OV")) {
    for(self in c(TRUE,FALSE)) {
      for(coln in name) {#adjust pvalues for each variable
        row <- (table_2can$self==self) & (table_2can$model==mod) & (table_2can$Cancer==cancer)
        table_2can[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_2can[row,paste("P_",coln,sep = "")],method = "fdr")
      }
    }
  }
}

table_sig <- table_2can[table_2can$FDR_pro_kin <= 0.1,]
table_sig_self <- table_sig[table_sig$self,]
table_sig_other <- table_sig[!table_sig$self,]

## extract top XX significant pairs in self-regulated or other group
top <- 10 # choose top n rows for P_pro_kin for model1
cancer <- "BRCA" #choose the cancer you want to sort
t0 <- table_sig_other# choose table_sig_other or table_sig_self
t1 <- data.table(t0[t0$model=="pho_sub~pro_kinase"& t0$Cancer==cancer,], key="FDR_pro_kin")# choose sort by which variable

## corresponding results for other two models are extracted and ordered
rows <- c()
for(i in 1:top){
  r <- unlist(which(t0$pair==t1$pair[i]))
  rows <- c(rows,r)
}
table <- t0[rows,]
table$order <- nrow(table):1
table$pairs <- reorder(table$pair,table$order)

temp <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","size","model","P_pro_kin","P_pro_sub","P_pho_kin","pair","Cancer","self","P_pro_kin_scale","order","pairs")
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

## actual plotting
p = ggplot(molten,aes(x=var, y=pairs))# make this the original ethni
p = p + facet_grid(.~Cancer+model, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point(aes(alpha =  abs(coefficients), colour=size, size =-log10(FDR))) + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p
