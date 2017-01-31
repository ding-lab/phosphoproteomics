# Yige Wu @ WashU 2017 Jan
# Pathway/subtype relations 
# Incorporate subtype information. Color the node as the difference of the phosphoprotein expression between subtypes. 
#? I use the grouped phosphorylation level for now

# library -----------------------------------------------------------------
library(ggplot2)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# choose one between the following two cancers to process -----------------
cancer = "BRCA"
# cancer = "OV"

sig <- 0.1
out_thres <- 1.5

# input according to cancer type-------------------------------------------------------------------
if (cancer == "BRCA") {
  # BRCA
  cancer = "BRCA"
  # BRCA_pro_f = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized_max10NA.txt",sep="")
  # pro_data <- read.delim(BRCA_pro_f)
  BRCA_pho_f = paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA.tsv",sep="")
  pho_data = read.delim(BRCA_pho_f)
  ## read in grouped phosphorylation data
  BRCA_pho_g = paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_by_PRO_formatted_normalized_max10NA.txt",sep="")
  pho_gdata = read.delim(BRCA_pho_g)
  clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))
  colx <- 78 # the column of protein name
}

if ( cancer == "OV" ) {
  #OV
  cancer = "OV"
  OV_pho_f = paste(baseD,"pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA.tsv",sep="")
  pho_data = read.delim(OV_pho_f)
  OV_pho_g = paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_by_PRO_formatted_normalized_max10NA.txt",sep="")
  pho_gdata = read.delim(OV_pho_g)
  OV_pro_f = paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized_max10NA.txt",sep="")
  pro_data <- read.delim(OV_pro_f)
  colnames(pro_data) <- str_split_fixed(colnames(pro_data),"_",2)[,1]
  pro_data <- pro_data[,colnames(pho_data)]
  colx <- 1 # the column of protein name
}

# ordering the columns by sample name
# pro_data <- pro_data[,order(names(pro_data))]
pho_data <- pho_data[,order(names(pho_data))]
pho_gdata <- pho_gdata[,order(names(pho_gdata))]#order the grouped phospho data

# #split the SUBSTRATE and SUB_MOD_RSD in the first column
# pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))
# 
# #covert the SUB_MOD_RSD from lowercase to uppercase
# pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
# colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")

# initiate----------------------------------------------------------------
x <- vector(mode = "numeric", length = nrow(pho_gdata))
pho_diff <- data.frame(matrix(rep(x,2), ncol=2, byrow=T))
rownames(pho_diff) <- pho_gdata$X
colnames(pho_diff) <- c("lum_mean","basal_mean")

# calculate the differential phospho level between subtypes ---------------
## get LumA&LumB-basal for all phosphorylation level
sample_names <- colnames(pho_gdata[,-colx])
pho_diff$lum_mean <- rowMeans(pho_gdata[,na.omit(sample_names[clinical[1,-1]=="LumA" | clinical[1,-1]=="LumB"])], na.rm = TRUE)
pho_diff$basal_mean <- rowMeans(pho_gdata[,na.omit(sample_names[clinical[1,-1]=="Basal"])], na.rm = TRUE)
pho_diff$lum_basal <- pho_diff$lum_mean - pho_diff$basal_mean

# input pathway info ----------------------------------------------
load("~/Box Sync/pan3can_shared_data/analysis_results/2015-08-01_Gene_Set.RData")


# output results for certain pathways -------------------------------------
path <- KEGG$`hsa04010	MAPK signaling pathway`
temp <- pho_diff[path,]
temp <- temp[!is.na(temp[,1]),]
path_pho_diff <- data.frame(temp$lum_basal)
rownames(path_pho_diff) <- row.names(temp)
tn = paste(baseD,"pan3can_shared_data/analysis_results/pathway/MAPK_lum_basal_phospho_diff.txt", sep="")
write.table(path_pho_diff, file=tn, quote=F, sep = '\t', row.names = TRUE)
