# HotPho.R
# Kuan Huang @ WashU 2017 March
# extract association data to find hot "subnetworks" of each kinase

# choose kinase/phosphotase, cancer , significance level, model -----------------------------------------------
protein <- "kinase"
sig = 0.05
strict_sig = 0.01

# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)

# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"

setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# input regression processed data -----------------------------------------
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/manuscript/supplementary_tables/kinase_substrate_regression_trans_edited.txt",sep = ""))

# trans
table_trans_pos = table_2can[!table_2can$self & table_2can$coef_pho_kin > 0,]

# recursive function
# to conduct statistics testing, we can likely permute the table and do this recursively to observe the size
kinase_subnet_recursive = function(kinase, subnetwork, heat_ini, thres){
  k_s = table_trans_pos[table_trans_pos$KINASE==kinase,]
  substrates = unique(k_s$SUBSTRATE)
  for (substrate in substrates){
    edge_score = 0
    FDR_pho_kin = min(k_s[k_s$SUBSTRATE==substrate,]$FDR_pho_kin, na.rm = T)
    edge_score = (1 - FDR_pho_kin)
    heat = heat_ini*edge_score
    if ( heat > thres ){
      subnetwork = c(subnetwork,substrate)
      # recursively find downstream subnetwork
      subnetwork = kinase_subnet_recursive(kinase = substrate, subnetwork = subnetwork, heat_ini = heat, thres = thres)
    }
  }
  return(subnetwork)
}

# track the heat dissipates by the kinase and find subnetworks for each kinase
all_kinases = unique(table_trans_pos$KINASE)
kinase_substrates = vector("list")
kinase_substrates_count = vector("list")
kinase_subnets = vector("list")
kinase_subnets_count = vector("list")
for (kinase in all_kinases){
  subnetwork = c(kinase)
  thres = 0.75
  k_substrates = unique(table_trans_pos$SUBSTRATE[table_trans_pos$KINASE==kinase & table_trans_pos$FDR_pho_kin < (1 - thres)])
  k_substrates_str = paste(k_substrates, collapse = ";")
  kinase_substrates[[kinase]] = k_substrates_str
  kinase_substrates_count[[kinase]] = length(k_substrates)
  
  k_subnetwork = kinase_subnet_recursive(kinase, subnetwork, heat_ini = 1, thres = thres)
  k_subnetwork_str = paste(k_subnetwork, collapse = ";")
  kinase_subnets[[kinase]] = k_subnetwork_str
  kinase_subnets_count[[kinase]] = length(k_subnetwork)
}
kinase_substrates_m = do.call(rbind,kinase_substrates)
kinase_substrates_count_m = do.call(rbind,kinase_substrates_count)

kinase_subnets_m = do.call(rbind,kinase_subnets)
kinase_subnets_count_m = do.call(rbind,kinase_subnets_count)

# save the resulting network
#save(kinase_subnets)

# some plotting of the subnetwork size
colnames(kinase_subnets_count_m) = "Subnetwork_size"
colnames(kinase_substrates_count_m) = "Substrate_count"
k_counts = merge(kinase_substrates_count_m, kinase_subnets_count_m, by="row.names")
colnames(k_counts)[1] = "Kinase"
k_counts$Subnetwork_nonsubstrate = k_counts$Subnetwork_size - k_counts$Substrate_count -1 
k_counts$ratio = k_counts$Subnetwork_nonsubstrate/k_counts$Subnetwork_size
k_counts_m = melt(k_counts,id.var=c("Kinase","Subnetwork_size","ratio"))

k_counts_m_sele = k_counts_m[k_counts_m$ratio>0.7,]
k_counts_m_sele$Kinase = factor(k_counts_m_sele$Kinase, 
                      levels = k_counts_m_sele$Kinase[order(k_counts_m_sele$Subnetwork_size, decreasing = T)])

p = ggplot(data=k_counts_m_sele)
p = p + geom_bar(aes(x=Kinase, y=value, fill=variable),alpha=0.8, stat="identity") #+ geom_point(aes(x=Status, y=value)) 
p = p + theme_bw() #+ guides(fill=FALSE) 
p = p + theme(axis.text.x=element_text(colour="black", size=12,angle=90,vjust=0.5))
p = p + scale_y_log10()
p = p + labs(y="Subnetwork member counts")
p
fn = paste(pd, "phospho_subnetwork.pdf", sep="_")
ggsave(file=fn, useDingbats=FALSE)

# kinase_subnets_count_m8 = data.frame(kinase_subnets_count_m[kinase_subnets_count_m>7,,drop=F])
# colnames(kinase_subnets_count_m8) = "Subnetwork_size"
# kinase_subnets_count_m8$kinase = row.names(kinase_subnets_count_m8)
# kinase_subnets_count_m8$kinase = factor(kinase_subnets_count_m8$kinase, 
#                                         levels = kinase_subnets_count_m8$kinase[order(kinase_subnets_count_m8$Subnetwork_size, decreasing = T)])
# 
# p = ggplot(data=kinase_subnets_count_m8)
# p = p + geom_bar(aes(x=kinase, y=Subnetwork_size, fill=as.factor(Subnetwork_size)),alpha=0.8, stat="identity") #+ geom_point(aes(x=Status, y=value)) 
# p = p + theme_bw() + guides(fill=FALSE) 
# p = p + theme(axis.text.x=element_text(colour="black", size=12,angle=90,vjust=0.5))
# p
# fn = paste(pd, "genelist_CNV_vs_pho_pro_ERBB2.pdf", sep="_")
# ggsave(file=fn, useDingbats=FALSE)