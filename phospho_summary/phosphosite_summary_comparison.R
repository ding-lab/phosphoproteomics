##### phosphosite_summary_comparison.R #####
# Kuan-lin Huang @ WashU 2016 Feb
# analysis on Broad phosphosite data

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/phospho_summary/")
#source("/Users/khuang/bin/LIB_exp.R")
library(biomaRt)
library(UpSetR)
library(matrixStats)
library(ggplot2)
source("../pan3can_aes.R")

get_sites = function(site_list){
  sites_cpdxered_type = c()
  for (site in site_list){
    type = substr(unlist(strsplit(site,"\\."))[2],1,1)#[1]
    sites_cpdxered_type = c(sites_cpdxered_type,type)
  }
  table(sites_cpdxered_type)
}

# read in data
pdb_phosphosites = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/HotPho/data_runs/PDB_phosphosite_wGpos.maf")
# pps_cancer = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Disease-associated_sites_term=cancer.oma.leuk.txt")
# pps_cancer = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Disease-associated_sites_term_cancer.oma.leuk_wGPos.txt")
# colnames(pps_cancer) <- c(colnames(pps_cancer)[-1],"x")
# pps_cancer$x <- NULL
# # biomart conversion: modularize this to a function sometime
# ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
# mapTab = getBM(attributes = c("hgnc_symbol", "uniprot_swissprot"), filters = "uniprot_swissprot", values = pps_cancer$ACC_ID, mart = ensembl, uniqueRows=FALSE)
# dupRows = union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
# if (length(dupRows)==0){mapTabNoDup = mapTab} else {mapTabNoDup = mapTab[-dupRows, ]}
# colnames(pps_cancer)[5]="uniprot_swissprot"
# pps_cancer_wgene = merge(pps_cancer, mapTabNoDup, by="uniprot_swissprot")
# pps_cancer_wgene$geneSite = paste(pps_cancer_wgene$hgnc_symbol,tolower(gsub("-p","",pps_cancer_wgene$MOD_RSD)),sep=".")
# pps_cancer_wgene = pps_cancer_wgene[-grep("p",pps_cancer_wgene$geneSite),]
# pps_cancer_wgene = pps_cancer_wgene[-grep("k",pps_cancer_wgene$geneSite),]
pps_cancer_wgene = read.table(header=T, quote="", sep = '\t', file="known_cancer_sites_cpdxerage_breakdown.tsv")

brca_file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/BRCA77_unimodal_phosphoproteome-ratio-norm_wGpos_cleaned.txt"
pdx_file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM_WUPCC_shared/WHIM\ Proteomic\ Data/WHIM_phosphoproteome/phosphoproteome-ratio-norm_exp_filtered.txt"
brca_pho = read.table(stringsAsFactors = F,row.names=NULL, header=TRUE, sep="\t", file= brca_file)
pdx_pho = read.table(stringsAsFactors = F,row.names=NULL, header=TRUE, sep="\t", file= pdx_file)

num_unique_sites = length(unique(c(brca_pho$genomic_pos,pdx_pho$genomic_pos)))
# genes = gsub("(.*?)\\.(.).*","\\1",row.names(brca_pho))
# sort(table(genes))
# BRCA
brca_siteSD = rowSds(as.matrix(brca_pho[,-c(1:3)]), na.rm=T)
brca_site_table = cbind(brca_pho$Gene.site,brca_siteSD,brca_pho$genomic_pos %in% pps_cancer_wgene$genomic_pos & brca_pho$genomic_pos != ".", "Human")
colnames(brca_site_table) = c("Phosphosite","SiteSD","Disease_site", "Cancer")
# pdx
pdx_pho = pdx_pho#[rowSums(!is.na(pdx_pho))>=10,]
pdx_siteSD = rowSds(as.matrix(pdx_pho[,-c(1:3)]), na.rm=T)
pdx_site_table = cbind(pdx_pho$gene.site,pdx_siteSD,pdx_pho$genomic_pos %in% pps_cancer_wgene$genomic_pos & pdx_pho$genomic_pos != ".", "PDX")
colnames(pdx_site_table) = c("Phosphosite","SiteSD","Disease_site", "Cancer")

site_table = as.data.frame(rbind(brca_site_table,pdx_site_table))
site_table$SiteSD = as.numeric(as.character(site_table$SiteSD))

#rule out the ones without genomic coordinate first
pps_cancer_wgene = pps_cancer_wgene[pps_cancer_wgene$genomic_pos != ".",]

sites_known = length(unique(pps_cancer_wgene$genomic_pos))
sites_covered = pps_cancer_wgene$genomic_pos[pps_cancer_wgene$genomic_pos %in% c(brca_pho$genomic_pos,pdx_pho$genomic_pos) & pps_cancer_wgene$genomic_pos != "."]
cat("Number of known cancer sites covered: ", length(unique(sites_covered)), "out of ", sites_known,"\n")

cat("Distribution of cancer sites: ","\n")
get_sites(unique(pps_cancer_wgene$geneSite))
cat("Distribution of cancer sites cpdxered: ","\n")
get_sites(unique(pps_cancer_wgene$geneSite[pps_cancer_wgene$genomic_pos %in% c(brca_pho$genomic_pos,pdx_pho$genomic_pos)& pps_cancer_wgene$genomic_pos != "."]))
# length(unique(pps_cancer_wgene$genomic_pos[pps_cancer_wgene$genomic_pos %in% c(brca_pho$genomic_pos,pdx_pho$genomic_pos)& pps_cancer_wgene$genomic_pos != "."]))

# # output the disease site table with each cohort status
# pdb_phosphosites_gpos = as.character(pdb_phosphosites$genomic_pos[!is.na(pdb_phosphosites$genomic_pos)])
# pps_cancer_wgene$In_CPTAC_BRCA = pps_cancer_wgene$genomic_pos %in% brca_pho$genomic_pos
# pps_cancer_wgene$In_CPTAC_pdx = pps_cancer_wgene$genomic_pos %in% pdx_pho$genomic_pos
# pps_cancer_wgene$In_PDB_phospho_site = pps_cancer_wgene$genomic_pos %in% pdb_phosphosites_gpos
# write.table(pps_cancer_wgene, col.names=NA, quote=F, sep = '\t', file="known_cancer_sites_cpdxerage_breakdown.tsv")

# # how many of these sites are cpdxered by RPPA data
# RPPA = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/collaborations/TCPA_RPPA_2015-10-30/TCGA-PANCAN16-RBN-Trans-Gene_sampleProbe_matrix.tsv")
# RPPA_sites = row.names(RPPA)[grep ("_p", row.names(RPPA))] # get only the phosphosites
# RPPA_sites_m = c()
# for (i in 1:length(RPPA_sites)){
#   site = tolower(sub("p","",unlist(strsplit(RPPA_sites[i], "_"))[3]))
#   gene = unlist(strsplit(RPPA_sites[i], "_"))[1]
#   gene_s = unlist(strsplit(gene, " "))
#   for (gene in gene_s){
#     RPPA_sites_m =c(RPPA_sites_m, paste(gene, site, sep=".")) 
#   }
# }
# #sum(sites_cpdxered %in% RPPA_sites_m)
# cat("Sites replicated by RPPA out of the cpdxered sites: ", sum(sites_cpdxered %in% RPPA_sites_m), "out of ", length(sites_cpdxered))

cat("Known disease site vs. non-disease site variance in BRCA")
wilcox.test(site_table[site_table$Disease_site == TRUE & site_table$Cancer == "Human",]$SiteSD,site_table[site_table$Disease_site ==FALSE & site_table$Cancer == "Human",]$SiteSD)
cat("Known disease site vs. non-disease site variance in pdx")
wilcox.test(site_table[site_table$Disease_site == TRUE & site_table$Cancer == "PDX",]$SiteSD,site_table[site_table$Disease_site ==FALSE & site_table$Cancer == "PDX",]$SiteSD)

site_table_merge = merge(brca_site_table,pdx_site_table,by="Phosphosite")
colnames(site_table_merge) = c("Phosphosite","SiteSD.Human", "Disease_site.Human", "Cancer.x", "SiteSD.PDX", "Disease_site.PDX", "Cancer.y")

site_table_merge_cancer = site_table_merge[site_table_merge$Disease_site.Human == TRUE,]
site_table_merge_cancer$SiteSD.Human = as.numeric(as.character(site_table_merge_cancer$SiteSD.Human))
site_table_merge_cancer$SiteSD.PDX = as.numeric(as.character(site_table_merge_cancer$SiteSD.PDX))
site_table_merge_cancer$Gene_site = gsub(":NP.*:"," ",site_table_merge_cancer$Phosphosite)

p = ggplot(data=site_table,aes(x=Disease_site, y=SiteSD))
p = p + facet_grid(.~Cancer)
p = p + geom_violin(aes(x=Disease_site, y=SiteSD, fill=Disease_site),alpha=0.5) + guides(fill=FALSE)
p = p + geom_jitter(aes(x=Disease_site, y=SiteSD, color=Disease_site),alpha = 0.2)
p = p + scale_colour_manual(name = "Disease Site", values = setNames(c('black',NA),c(T, F)))
p = p + geom_text(aes(x=Disease_site, y=SiteSD, label = ifelse(Disease_site & SiteSD > 1.5,Gene_site,NA)),size=2)
# p = p + labs(x = "Known cancer phosphosites", y = "Standard Deviation")
p = p + stat_summary(fun.y="median")
p = p + theme_bw() + theme_nogrid()
p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
              axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8),
              legend.position="none")
p = p + labs(y="Phosphosite standard deviation", x="Known cancer site")
p
fn = paste(pd,"BRCA_PDX_siteSD_comparison.pdf",sep="_")
ggsave(file=fn, height=5, useDingbats=FALSE)

# ### upsetR between all the universe ###
# all_sites = data.frame(unique(as.character(c(brca_pho$genomic_pos,pdx_pho$genomic_pos,pps_cancer_wgene$genomic_pos,pdb_phosphosites_gpos,pho_elm_phosphosites_gpos))))
# colnames(all_sites)[1] = "Phosphosite"
# 
# #all_sites$Phosphosite = as.character(all_sites$Phosphosite)
# #all_sites = all_sites[all_sites$Phosphosite != ".",]
# all_sites$CPTAC_HUMAN = all_sites$Phosphosite %in% brca_pho$genomic_pos
# all_sites$CPTAC_PDX = all_sites$Phosphosite %in% pdx_pho$genomic_pos
# all_sites$Known_cancer_site = all_sites$Phosphosite %in% pps_cancer_wgene$genomic_pos
# all_sites$PDB_phospho_site = all_sites$Phosphosite %in% pdb_phosphosites_gpos
# all_sites$Phospho_ELM_site = all_sites$Phosphosite %in% pho_elm_phosphosites_gpos
# all_sites = all_sites[all_sites$Phosphosite != ".",]
# 
# all_sites$CPTAC_HUMAN[all_sites$CPTAC_HUMAN]=1
# all_sites$CPTAC_HUMAN[!all_sites$CPTAC_HUMAN]=0
# all_sites$CPTAC_PDX[all_sites$CPTAC_PDX]=1
# all_sites$CPTAC_PDX[!all_sites$CPTAC_PDX]=0
# all_sites$Known_cancer_site[all_sites$Known_cancer_site]=1
# all_sites$Known_cancer_site[!all_sites$Known_cancer_site]=0
# all_sites$PDB_phospho_site[all_sites$PDB_phospho_site]=1
# all_sites$PDB_phospho_site[!all_sites$PDB_phospho_site]=0
# all_sites$Phospho_ELM_site[all_sites$Phospho_ELM_site]=1
# all_sites$Phospho_ELM_site[!all_sites$Phospho_ELM_site]=0
# fn = paste(pd,"phosphosite_upset.pdf",sep="_")
# pdf(fn, height=4, useDingbats = F)
# upset(all_sites,sets = c("CPTAC_HUMAN","CPTAC_PDX","Known_cancer_site","PDB_phospho_site","Phospho_ELM_site" ),order.by = "freq")
# dev.off()
# 
# # PDB
# pdb_phosphosites_non_unique = gsub("p.","",pdb_phosphosites$amino_acid_change)
# pdb_phosphosites_residue = substring(pdb_phosphosites_non_unique, 1, 1)
# pdb_phosphosites_residue = pdb_phosphosites_residue[pdb_phosphosites_residue!="H"]
# pdb_phosphosites_r_count = table(pdb_phosphosites_residue)
# pdb_phosphosites_r_count
# 
# # phospho.ELM
# pho_elm_phosphosites_sites = gsub(".*:","",pho_elm_phosphosites$input)
# pho_elm_phosphosites_residue = substring(pho_elm_phosphosites_sites, 1, 1)
# pho_elm_phosphosites_residue = pho_elm_phosphosites_residue[pho_elm_phosphosites_residue!="H"]
# pho_elm_phosphosites_r_count = table(pho_elm_phosphosites_residue)
# pho_elm_phosphosites_r_count
# 
# # BRCA
# non_unique = gsub(".*:","",row.names(brca_pho))
# brca_residue = substring(non_unique, 1, 1)
# BRCA = table(brca_residue)
# BRCA
# 
# # OV
# non_unique = gsub(".*:","",row.names(ov_pho))
# OV_residue = substring(non_unique, 1, 1)
# OV = table(OV_residue)
# OV
# 
# # both
# residue_table = rbind(BRCA,OV,pdb_phosphosites_r_count,pho_elm_phosphosites_r_count)
# residue_table_m = melt(residue_table)
# colnames(residue_table_m) = c("Cancer","Residue","Number_of_Phosphosites")
# residue_table_m$Cancer = c("CPTAC BRCA","CPTAC OV", "PDB","Phospho.ELM", 
#                            "CPTAC BRCA","CPTAC OV", "PDB","Phospho.ELM",
#                            "CPTAC BRCA","CPTAC OV", "PDB","Phospho.ELM")
# 
# p = ggplot(residue_table_m,aes(x=Cancer, y=Number_of_Phosphosites, fill=Residue))
# #p = p + facet_grid(Cancer~.,drop=T,space="free",scale="free")
# p = p + geom_bar(stat="identity") + theme_bw() #+ theme_nogrid()
# p = p + labs(x = "Data", y="# of phosphosites")
# p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
# p = p + theme(legend.position="bottom")
# p = p + coord_flip()
# p
# fn = paste(pd, 'phosphosite_brca.ov_sty_summary.pdf',sep ="_")
# ggsave(file=fn, height = 4, width=4, useDingbats=FALSE)
# 
# ### distribution among kinase subgroups and families
