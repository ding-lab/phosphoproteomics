### bubble_plot_phosphosite.R ### 
# Yige Wu @ WashU 2016 Nov
# plot correlation p-value and sample size between substrate phosphorylation level and kinase protein expression level, etc
# use make_table_for_bubble.R to prepapre the table

setwd("~/proteomics/pan3can_analysis/phospho_network")
library(ggplot2)
library(ggtern)

# keep significant ones
table_sig <- table[table_2can$Pvalue <= 0.05,]

# plot all singnificant pairs
p = ggplot(table_sig,aes(x=model, y=pair, colour=-log(Pvalue), size =size))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point() + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p

# plot self-regulated in significant ones
#width=421, height = 800 for one cancer
p = ggplot(table_sig[table_sig$self,],aes(x=model, y=pair, colour=-log(Pvalue), size =size))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point() + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p

# plot non self-regulated separately in significant ones
#width=500, height = 3000 for one cancer
p = ggplot(table_sig[!table_sig$self,],aes(x=model, y=pair, colour=-log(Pvalue), size =size))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point() + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p



