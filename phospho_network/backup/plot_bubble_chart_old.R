### bubble_plot_phosphosite.R ### 
# Yige Wu @ WashU 2016 Nov
# plot correlation p-value and sample size between substrate phosphorylation level and kinase protein expression level, etc
# use make_table_for_bubble.R to prepapre the table

setwd("~/proteomics/pan3can_analysis/phospho_network")
load("~/proteomics/pan3can_analysis/phospho_network/BRCA_3models.RData")
load("~/proteomics/pan3can_analysis/phospho_network/OV_3models.RData")
library(ggplot2)
library(ggtern)
require(data.table)
library(reshape)


# merge table for BRCA and OV, keep significant ones
table_2can <- rbind(table_BRCA,table_OV)
table_2can$P_pro_kin_scale <- round(-log10(table_2can$P_pro_kin),digit=0)
table_sig <- table_2can[table_2can$P_pro_kin <= 0.05,]
table_sig_self <- table_sig[table_sig$self,]
table_sig_other <- table_sig[!table_sig$self,]

# plot top significant ones for self-regualted
top <- 50 # choose top n rows for P_pro_kin for model1
cancer <- "BRCA" #choose the cancer you want to sort
# sort table_sig first by P_pro_kin
t0 <- table_sig_other
t1 <- data.table(t0[t0$model=="pho_sub~pro_kinase"& t0$Cancer==cancer,], key="P_pro_kin")

rows <- c()
for(i in 1:top){
  r <- unlist(which(t0$pair==t1$pair[i]))
  rows <- c(rows,r)
}

table <- t0[rows,]
table$order <- nrow(table):1
table$pairs <- reorder(table$pair,table$order)
p = ggplot(table,aes(x=model, y=pairs, colour=-log10(P_pro_kin), size =size))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point() + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p



# plot self-regulated in significant ones
#width=421, height = 800 for one cancer
p = ggplot(table_sig[table_sig$self,],aes(x=model, y=pair, colour=-log(P_pro_kin), size =size))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point() + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p

# plot non self-regulated separately in significant ones
#width=500, height = 3000 for one cancer
p = ggplot(table_sig[!table_sig$self,],aes(x=model, y=pair, colour=-log(P_pro_kin), size =size))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point() + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p


# plot all singnificant pairs
p = ggplot(table_sig,aes(x=model, y=pair, colour=-log(P_pro_kin), size =size))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point() + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
