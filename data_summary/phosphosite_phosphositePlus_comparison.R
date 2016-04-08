##### phosphosite_summary.R #####
# Kuan-lin Huang @ WashU 2016 Feb
# cross-correlation between different phosphosites

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/data_summary/")
source("/Users/khuang/bin/LIB_exp.R")
library(biomaRt)

get_sites = function(site_list){
  sites_covered_type = c()
  for (site in site_list){
    type = substr(unlist(strsplit(site,"\\."))[2],1,1)#[1]
    sites_covered_type = c(sites_covered_type,type)
  }
  table(sites_covered_type)
}

# read in 
pps_cancer = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/PhosphositePlus/data/Disease-associated_sites_term=cancer.oma.leuk.txt")
colnames(pps_cancer) <- c(colnames(pps_cancer)[-1],"x")
pps_cancer$x <- NULL
# biomart conversion: modularize this to a function sometime
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
mapTab = getBM(attributes = c("hgnc_symbol", "uniprot_swissprot"), filters = "uniprot_swissprot", values = pps_cancer$ACC_ID, mart = ensembl, uniqueRows=FALSE)
dupRows = union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
if (length(dupRows)==0){mapTabNoDup = mapTab} else {mapTabNoDup = mapTab[-dupRows, ]}
colnames(pps_cancer)[4]="uniprot_swissprot"
pps_cancer_wgene = merge(pps_cancer, mapTabNoDup, by="uniprot_swissprot")
pps_cancer_wgene$geneSite = paste(pps_cancer_wgene$hgnc_symbol,tolower(gsub("-p","",pps_cancer_wgene$MOD_RSD)),sep=".")
pps_cancer_wgene = pps_cancer_wgene[-grep("p",pps_cancer_wgene$geneSite),]
pps_cancer_wgene = pps_cancer_wgene[-grep("k",pps_cancer_wgene$geneSite),]

brca_file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
ov_file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/data/201602_pancan_phospho/TCGA_Breast_BI_Phosphoproteome/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"
brca_pho = read.table(row.names=1, header=TRUE, sep="\t", file= brca_file)
ov_pho = read.table(row.names=1, header=TRUE, sep="\t", file= ov_file)

# genes = gsub("(.*?)\\.(.).*","\\1",row.names(brca_pho))
# sort(table(genes))
# BRCA
brca_pho_10NA = brca_pho[rowSums(!is.na(brca_pho))>=10,]
brca_siteSD = rowSds(as.matrix(brca_pho_10NA), na.rm=T)
brca_site_table = cbind(row.names(brca_pho_10NA),brca_siteSD,row.names(brca_pho_10NA) %in% pps_cancer_wgene$geneSite, "BRCA")
colnames(brca_site_table) = c("Phosphosite","SiteSD","Disease_site", "Cancer")
# OV
ov_pho_10NA = ov_pho[rowSums(!is.na(ov_pho))>=10,]
ov_siteSD = rowSds(as.matrix(ov_pho_10NA), na.rm=T)
ov_site_table = cbind(row.names(ov_pho_10NA),ov_siteSD,row.names(ov_pho_10NA) %in% pps_cancer_wgene$geneSite, "OV")
colnames(ov_site_table) = c("Phosphosite","SiteSD","Disease_site", "Cancer")

site_table = as.data.frame(rbind(brca_site_table,ov_site_table))
site_table$SiteSD = as.numeric(as.character(site_table$SiteSD))

sites_known = length(unique(pps_cancer_wgene$geneSite))
sites_covered = pps_cancer_wgene$geneSite[pps_cancer_wgene$geneSite %in% c(row.names(brca_pho),row.names(ov_pho))]
cat("Number of known cancer sites covered: ", length(sites_covered), "out of ", sites_known)

get_sites(unique(pps_cancer_wgene$geneSite))
get_sites(sites_covered)

# how many of these sites are covered by RPPA data
RPPA = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/collaborations/TCPA_RPPA_2015-10-30/TCGA-PANCAN16-RBN-Trans-Gene_sampleProbe_matrix.tsv")
RPPA_sites = row.names(RPPA)[grep ("_p", row.names(RPPA))] # get only the phosphosites
RPPA_sites_m = c()
for (i in 1:length(RPPA_sites)){
  site = tolower(sub("p","",unlist(strsplit(RPPA_sites[i], "_"))[3]))
  gene = unlist(strsplit(RPPA_sites[i], "_"))[1]
  gene_s = unlist(strsplit(gene, " "))
  for (gene in gene_s){
    RPPA_sites_m =c(RPPA_sites_m, paste(gene, site, sep=".")) 
  }
}
#sum(sites_covered %in% RPPA_sites_m)
cat("Sites replicated by RPPA out of the covered sites: ", sum(sites_covered %in% RPPA_sites_m), "out of ", length(sites_covered))

cat("Known disease site vs. non-disease site variance in BRCA")
wilcox.test(site_table[site_table$Disease_site == TRUE & site_table$Cancer == "BRCA",]$SiteSD,site_table[site_table$Disease_site ==FALSE & site_table$Cancer == "BRCA",]$SiteSD)
cat("Known disease site vs. non-disease site variance in OV")
wilcox.test(site_table[site_table$Disease_site == TRUE & site_table$Cancer == "OV",]$SiteSD,site_table[site_table$Disease_site ==FALSE & site_table$Cancer == "OV",]$SiteSD)

fn = paste(pd, "phosphosite_SD_by_disease_violin.pdf", sep="_")
p = ggplot(data=site_table)
p = p + facet_grid(~Cancer)
p = p + geom_violin(aes(x=Disease_site, y=SiteSD, fill=Disease_site),alpha=0.5) + guides(fill=FALSE) 
#p = p + geom_point(aes(x=Mutation_Status, y=Expression, color=Mutation_Type, position = position_jitter(w = 0.1, h = 0)))
p = p + geom_jitter(aes(x=Disease_site, y=SiteSD, fill=Disease_site), alpha=0.05) #+ geom_point(aes(x=Status, y=value)) 
#p = p + geom_dotplot(aes(binaxis = "y",x=Mutation_Status, y=as.numeric(Expression), color=Mutation_Type, 
#                         stackdir = "center",binwidth = 5, stackgroups = TRUE,method = "histodot")) #+ geom_point(aes(x=Status, y=value)) 
p = p + labs(x = "Known cancer phosphosites", y = "Standard Deviation")
#p = p + ylim(c(min_bound,max_bound))
p = p + theme_bw() + theme_nogrid()
p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
              axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
p
ggsave(file=fn, limitsize=FALSE, useDingbats=FALSE)

