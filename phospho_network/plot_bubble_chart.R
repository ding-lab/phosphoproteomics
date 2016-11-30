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
table_sig <- table_2can[table_2can$P_pro_kin <= 0.05,]
table_sig_self <- table_sig[table_sig$self,]
table_sig_other <- table_sig[!table_sig$self,]

## extract top XX significant pairs in self-regulated or other group
top <- 50 # choose top n rows for P_pro_kin for model1
cancer <- "BRCA" #choose the cancer you want to sort
t0 <- table_sig_other# choose table_sig_other or table_sig_self
t1 <- data.table(t0[t0$model=="pho_sub~pro_kinase"& t0$Cancer==cancer,], key="P_pro_kin")# choose sort by which variable

## corresponding results for other two models are extracted and ordered
rows <- c()
for(i in 1:top){
  r <- unlist(which(t0$pair==t1$pair[i]))
  rows <- c(rows,r)
}
table <- t0[rows,]
table$order <- nrow(table):1
table$pairs <- reorder(table$pair,table$order)
molten = melt(table, id = c("KINASE","SUBSTRATE","SUB_MOD_RSD","size","model","P_pro_kin","P_pro_sub","P_pho_kin","pair","Cancer","self","P_pro_kin_scale","order","pairs"))

## actual plotting
p = ggplot(molten,aes(x=variable, y=pairs))# make this the original ethni
p = p + facet_grid(.~Cancer+model, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point(aes(alpha = size, colour=log10(P_pro_kin), size = abs(value))) + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=5))#element_text(colour="black", size=14))
p
