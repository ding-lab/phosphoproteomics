# Yige Wu @ WashU 2017 Jan
# For RTKs,plot kinase protein expression ~ substrate average phosphorylation score, marked with protein_outlier, colored with kinase phosphorylation level

# library -----------------------------------------------------------------
library(ggplot2)
library(readr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# score~pro_kin for RTKs --------------------------------------------------
overlap <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/kinase_BRCA_reverse_marking_out_thers_1.5.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
rtk <- as.vector(t(RTK_file))

row <- c()
for (kinase in rtk) {
  temp <- which(overlap$kinase==kinase & !is.na(overlap$score) & !is.na(overlap$pro_level) )
  row <- c(row,temp)
}
table <- overlap[row,]
table$transparency <- vector(mode = "numeric", length = nrow(table))
table$transparency[table$pro_outlier_status] <- table$transparency[table$pro_outlier_status] + 1
table$transparency[table$missense] <- table$transparency[table$missense] + 1
table$transparency <- (table$transparency+3)/5

lim = max(abs(max(table$kin_phos_score, na.rm = T)),abs(min(table$kin_phos_score, na.rm = T)))
# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
getPalette = colorRampPalette(c("#ffffcc","#fd8d3c","#800026"))

outlier.colors=c("NA", "#000000")
p = ggplot(data = table,aes(x= pro_level, y = score, shape = pro_outlier_status, colour = kin_phos_score , alpha = 1))
p = p + scale_colour_gradientn(name= "kin_phos_score", na.value="grey40", colours=getPalette(100)) # , limit=c(-lim,lim)
p = p + geom_point(size = 2, alpha=0.8, stroke = 0) 
p = p + geom_text(data=subset(table, aa != "wt" ),
            aes(x= pro_level, y = score,label=aa), hjust = 0, nudge_x = -1.2 ,nudge_y = 0.2, size = 2) +
  facet_wrap(~ kinase, nrow = 4)
p = p + theme_bw()
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/BRCA_RTK_score~pro_level_aa_change_kin_phos_score.pdf',sep ="")
ggsave(file=fn, height=10, width=13)

p = ggplot(data = table,aes(x= pro_level, y = score, colour = kinase , shape = missense.pro_outlier)) +
  geom_point(alpha=0.2) +
  geom_text(data=subset(table, pro_outlier_status | missense ),
            aes(x= pro_level, y = score,label=kinase),hjust = 0, nudge_x = -0.3 ,nudge_y = 0.1, size=2, alpha=0.2)
p = p + theme_bw()
p = p + labs(x="Protein expression", y = "Substrate Phosphorylation Score")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/RTK_BRCA_score~pro_level.pdf',sep ="")
ggsave(file=fn, height=4, width=5)


