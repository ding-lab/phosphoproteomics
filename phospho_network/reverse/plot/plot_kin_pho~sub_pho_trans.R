# Yige Wu @ WashU 2017 Feb
# For RTKs,plot protein expression ~ collapsed phosphorylation score, marked with mutation, protein_outlier, upstream kinase phosphorylation level(not include rtk itself)

# library -----------------------------------------------------------------
library(ggplot2)
library(readr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# input -------------------------------------------------------------------
# overlap <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/kinase_BRCA_reverse_marking_out_thers_1.5.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
overlap <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/kinase_BRCA_reverse_marking_median_out_thers_1.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)

RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
rtk <- as.vector(t(RTK_file))

# pull out result for rtks ------------------------------------------------
row <- c()
for (kinase in rtk) {
  temp <- which(overlap$kinase==kinase & !is.na(overlap$kin_phos_score) & !is.na(overlap$sub_phos_transscore) )
  row <- c(row,temp)
}
table <- overlap[row,]


# score~pro_kin for RTKs --------------------------------------------------
lim = max(abs(max(table$up_phos_score, na.rm = T)),abs(min(table$up_phos_score, na.rm = T)))
# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
getPalette = colorRampPalette(c("#ffffcc","#fd8d3c","#800026"))
outlier.colors=c("NA", "#000000")

p = ggplot(data = table,aes(x= sub_phos_transscore, y = kin_phos_score, shape = missense.pro_outlier, colour = up_phos_score , alpha = 1))
p = p + scale_colour_gradientn(name= "up_phos_score", na.value="grey40", colours=getPalette(100)) # , limit=c(-lim,lim)
p = p + geom_point(size = 2, alpha=0.8, stroke = 0) 
p = p + geom_text(data=subset(table, aa != "wt" ),
                  aes(x= pro_level, y = score,label=aa), hjust = 0, nudge_x = -1.2 ,nudge_y = 0.2, size = 2) +
  facet_wrap(~ kinase, nrow = 4)
p = p + theme_bw()
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/BRCA_RTK_kin_phos_score~sub_phos_transscore_.pdf',sep ="")
ggsave(file=fn, height=10, width=10)

# ggplot(data = table,aes(x= pro_level, y = score, colour = kinase , shape = missense.pro_outlier, alpha = (transparency+1)/3)) +    
#   geom_point(size = 2) + 
#   geom_text(data=subset(table, pro_outlier_status | missense ),
#             aes(x= pro_level, y = score,label=kinase),hjust = 0, nudge_x = -0.3 ,nudge_y = 0.1)
# fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/RTK_BRCA_score~pro_level.pdf',sep ="")
# ggsave(file=fn, height=8, width=10)


