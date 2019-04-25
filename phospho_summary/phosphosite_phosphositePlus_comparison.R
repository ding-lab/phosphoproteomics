##### phosphosite_summary.R #####
# Kuan-lin Huang @ WashU 2016 Feb
# cross-correlation between different phosphosites

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/phospho_summary/")
source("/Users/khuang/bin/LIB_exp.R")
library(biomaRt)
library(UpSetR)

get_sites = function(site_list){
  sites_covered_type = c()
  for (site in site_list){
    type = substr(unlist(strsplit(site,"\\."))[2],1,1)#[1]
    sites_covered_type = c(sites_covered_type,type)
  }
  table(sites_covered_type)
}

# read in data
psp_phosphosites = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Phosphorylation_site_dataset_transvar_in_uniq.txt")

pdb_phosphosites = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/HotPho/data_runs/PDB_phosphosite_wGpos.maf")
pps_cancer = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Disease-associated_sites_term=cancer.oma.leuk.txt")
pps_cancer = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Disease-associated_sites_term_cancer.oma.leuk_wGPos.txt")
colnames(pps_cancer) <- c(colnames(pps_cancer)[-1],"x")
pps_cancer$x <- NULL
# biomart conversion: modularize this to a function sometime
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hmale_gene_ensembl",host="www.ensembl.org")
#ensembl = useMart(biomart = "ENSEMBL_MART_HUMAN")
mapTab = getBM(attributes = c("hgnc_symbol", "uniprotswissprot"), filters = "uniprotswissprot", values = pps_cancer$ACC_ID, mart = ensembl, uniqueRows=FALSE)
dupRows = union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
if (length(dupRows)==0){mapTabNoDup = mapTab} else {mapTabNoDup = mapTab[-dupRows, ]}
colnames(pps_cancer)[5]="uniprot_swissprot"
pps_cancer_wgene = merge(pps_cancer, mapTabNoDup, by="uniprot_swissprot")
pps_cancer_wgene$geneSite = paste(pps_cancer_wgene$hgnc_symbol,tolower(gsub("-p","",pps_cancer_wgene$MOD_RSD)),sep=".")
pps_cancer_wgene = pps_cancer_wgene[-grep("p",pps_cancer_wgene$geneSite),]
pps_cancer_wgene = pps_cancer_wgene[-grep("k",pps_cancer_wgene$geneSite),]

brca_file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA_wGpos.tsv"
ov_file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA_wGpos.tsv"
brca_pho = read.table(row.names=2, header=TRUE, sep="\t", file= brca_file)
ov_pho = read.table(row.names=2, header=TRUE, sep="\t", file= ov_file)

brca_pho_genome = brca_pho$row.names
ov_pho_genome = ov_pho$row.names
brca_pho = brca_pho[,-1]
ov_pho = ov_pho[,-1]
brca_pho_genome_10NA = brca_pho_genome#[rowSums(!is.na(brca_pho))>=10]
ov_pho_genome_10NA = ov_pho_genome#[rowSums(!is.na(ov_pho))>=10]

num_unique_sites = length(unique(c(brca_pho_genome_10NA,ov_pho_genome_10NA)))
# genes = gsub("(.*?)\\.(.).*","\\1",row.names(brca_pho))
# sort(table(genes))
# BRCA
brca_pho_10NA = brca_pho#[rowSums(!is.na(brca_pho))>=10,]
brca_siteSD = rowSds(as.matrix(brca_pho_10NA), na.rm=T)
brca_site_table = cbind(row.names(brca_pho_10NA),brca_siteSD,brca_pho_genome_10NA %in% pps_cancer_wgene$genomic_pos & brca_pho_genome_10NA != ".", "BRCA")
colnames(brca_site_table) = c("Phosphosite","SiteSD","Disease_site", "Cancer")
# OV
ov_pho_10NA = ov_pho#[rowSums(!is.na(ov_pho))>=10,]
ov_siteSD = rowSds(as.matrix(ov_pho_10NA), na.rm=T)
ov_site_table = cbind(row.names(ov_pho_10NA),ov_siteSD,ov_pho_genome_10NA %in% pps_cancer_wgene$genomic_pos & ov_pho_genome_10NA != ".", "OV")
colnames(ov_site_table) = c("Phosphosite","SiteSD","Disease_site", "Cancer")

site_table = as.data.frame(rbind(brca_site_table,ov_site_table))
site_table$SiteSD = as.numeric(as.character(site_table$SiteSD))

# Count the ones in PSP
colnames(psp_phosphosites) = "sites"
psp_phosphosites = psp_phosphosites[-grep("^:",psp_phosphosites$sites),]
cat("PSP phosphosites:\n")
length(psp_phosphosites)
brca_pho_10NA_sites = toupper(gsub(":NP.*:",":",row.names(brca_pho_10NA)))
cat("BRCA phosphosites in PSP:\n")
table(brca_pho_10NA_sites %in% gsub("p.","",psp_phosphosites))

#rule out the ones without genomic coordinate first
pps_cancer_wgene = pps_cancer_wgene[pps_cancer_wgene$genomic_pos != ".",]

sites_known = length(unique(pps_cancer_wgene$genomic_pos))
sites_covered = pps_cancer_wgene$genomic_pos[pps_cancer_wgene$genomic_pos %in% c(brca_pho_genome,ov_pho_genome) & pps_cancer_wgene$genomic_pos != "."]
cat("Number of known cancer sites covered: ", length(unique(sites_covered)), "out of ", sites_known,"\n")

cat("Distribution of cancer sites: ","\n")
get_sites(unique(pps_cancer_wgene$geneSite))
cat("Distribution of cancer sites covered: ","\n")
get_sites(unique(pps_cancer_wgene$geneSite[pps_cancer_wgene$genomic_pos %in% c(brca_pho_genome,ov_pho_genome)& pps_cancer_wgene$genomic_pos != "."]))
# length(unique(pps_cancer_wgene$genomic_pos[pps_cancer_wgene$genomic_pos %in% c(brca_pho_genome,ov_pho_genome)& pps_cancer_wgene$genomic_pos != "."]))

# output the disease site table with each cohort status
pdb_phosphosites_gpos = as.character(pdb_phosphosites$genomic_pos[!is.na(pdb_phosphosites$genomic_pos)])
pps_cancer_wgene$In_CPTAC_BRCA = pps_cancer_wgene$genomic_pos %in% brca_pho_genome_10NA
pps_cancer_wgene$In_CPTAC_OV = pps_cancer_wgene$genomic_pos %in% ov_pho_genome_10NA
pps_cancer_wgene$In_PDB_phospho_site = pps_cancer_wgene$genomic_pos %in% pdb_phosphosites_gpos
write.table(pps_cancer_wgene, col.names=NA, quote=F, sep = '\t', file="known_cancer_sites_coverage_breakdown.tsv")

# # how many of these sites are covered by RPPA data
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
# #sum(sites_covered %in% RPPA_sites_m)
# cat("Sites replicated by RPPA out of the covered sites: ", sum(sites_covered %in% RPPA_sites_m), "out of ", length(sites_covered))

cat("Known disease site vs. non-disease site variance in BRCA")
wilcox.test(site_table[site_table$Disease_site == TRUE & site_table$Cancer == "BRCA",]$SiteSD,site_table[site_table$Disease_site ==FALSE & site_table$Cancer == "BRCA",]$SiteSD)
cat("Known disease site vs. non-disease site variance in OV")
wilcox.test(site_table[site_table$Disease_site == TRUE & site_table$Cancer == "OV",]$SiteSD,site_table[site_table$Disease_site ==FALSE & site_table$Cancer == "OV",]$SiteSD)

# fn = paste(pd, "phosphosite_SD_by_disease_violin.pdf", sep="_")
# p = ggplot(data=site_table)
# #p = p + facet_grid(~Cancer)
# p = p + geom_violin(aes(x=Disease_site, y=SiteSD, fill=Disease_site),alpha=0.5) + guides(fill=FALSE)
# p = p + geom_jitter(aes(x=Disease_site, y=SiteSD, color=ifelse(as.logical(site_table$Disease_site),"black",NA)), alpha=0.05)
# p = p + scale_colour_gradientn(na.value=NA,colors=c("black"))
# p = p + labs(x = "Known cancer phosphosites", y = "Standard Deviation")
# p = p + theme_bw() + theme_nogrid()
# p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
#               axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
# p
# ggsave(file=fn, limitsize=FALSE, useDingbats=FALSE)

site_table_merge = merge(brca_site_table,ov_site_table,by="Phosphosite")
colnames(site_table_merge) = c("Phosphosite","SiteSD.BRCA", "Disease_site.BRCA", "Cancer.x", "SiteSD.OV", "Disease_site.OV", "Cancer.y")

site_table_merge_cancer = site_table_merge[site_table_merge$Disease_site.BRCA == TRUE,]
site_table_merge_cancer$SiteSD.BRCA = as.numeric(as.character(site_table_merge_cancer$SiteSD.BRCA))
site_table_merge_cancer$SiteSD.OV = as.numeric(as.character(site_table_merge_cancer$SiteSD.OV))
site_table_merge_cancer$Gene_site = gsub(":NP.*:","",site_table_merge_cancer$Phosphosite)

p = ggplot(data=site_table_merge_cancer)
p = p + geom_point(aes(x=SiteSD.BRCA, y=SiteSD.OV, fill=Disease_site.BRCA),alpha=0.5) + guides(fill=FALSE)
#p = p + geom_jitter(aes(x=Disease_site, y=SiteSD, color=ifelse(as.logical(site_table$Disease_site),"black",NA)), alpha=0.05)
# p = p + scale_colour_gradientn(na.value=NA,colors=c("black"))
# p = p + labs(x = "Known cancer phosphosites", y = "Standard Deviation")
p = p + geom_text(aes(x=SiteSD.BRCA, y=SiteSD.OV, label = Gene_site),vjust=-0.2,alpha=0.8)
p = p + theme_bw() + theme_nogrid()
p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
              axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
p = p + coord_equal() + xlim(0,2) + ylim(0,2)
p
fn = paste(pd,"BRCA_OV_cancerOnly_siteSD.pdf",sep="_")
ggsave(file=fn, useDingbats=FALSE)

site_table$Disease_site = as.logical(site_table$Disease_site)
site_table$Gene_site = gsub(":NP.*:"," ",site_table$Phosphosite)

p = ggplot(data=site_table)
p = p + facet_grid(.~Cancer)
p = p + geom_violin(aes(x=Disease_site, y=SiteSD, fill=Disease_site),alpha=0.5) + guides(fill=FALSE)
p = p + geom_jitter(aes(x=Disease_site, y=SiteSD, color=Disease_site),alpha = 0.2)
p = p + scale_colour_manual(name = "Disease Site", values = setNames(c('black',NA),c(T, F)))
p = p + geom_text(aes(x=Disease_site, y=SiteSD, label = ifelse(Disease_site & SiteSD > 1.5,Gene_site,NA)),size=2)
# p = p + labs(x = "Known cancer phosphosites", y = "Standard Deviation")
p = p + theme_bw() + theme_nogrid()
p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
              axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8),
              legend.position="none")
p = p + labs(y="Phosphosite standard deviation", x="Known cancer site")
p
fn = paste(pd,"BRCA_OV_siteSD_comparison.pdf",sep="_")
ggsave(file=fn, height=5, useDingbats=FALSE)

### upsetR between all the universe ###
all_sites = data.frame(unique(as.character(c(brca_pho_genome_10NA,ov_pho_genome_10NA,pps_cancer_wgene$genomic_pos,pdb_phosphosites_gpos))))
colnames(all_sites)[1] = "Phosphosite"

#all_sites$Phosphosite = as.character(all_sites$Phosphosite)
#all_sites = all_sites[all_sites$Phosphosite != ".",]
all_sites$CPTAC_BRCA = all_sites$Phosphosite %in% brca_pho_genome_10NA
all_sites$CPTAC_OV = all_sites$Phosphosite %in% ov_pho_genome_10NA
all_sites$Known_cancer_site = all_sites$Phosphosite %in% pps_cancer_wgene$genomic_pos
all_sites$PDB_phospho_site = all_sites$Phosphosite %in% pdb_phosphosites_gpos
all_sites$PSP_phospho_site = all_sites$Phosphosite %in% XXX
all_sites = all_sites[all_sites$Phosphosite != ".",]

all_sites$CPTAC_BRCA[all_sites$CPTAC_BRCA]=1
all_sites$CPTAC_BRCA[!all_sites$CPTAC_BRCA]=0
all_sites$CPTAC_OV[all_sites$CPTAC_OV]=1
all_sites$CPTAC_OV[!all_sites$CPTAC_OV]=0
all_sites$Known_cancer_site[all_sites$Known_cancer_site]=1
all_sites$Known_cancer_site[!all_sites$Known_cancer_site]=0
all_sites$PDB_phospho_site[all_sites$PDB_phospho_site]=1
all_sites$PDB_phospho_site[!all_sites$PDB_phospho_site]=0
all_sites$PSP_phospho_site[all_sites$PSP_phospho_site]=1
all_sites$PSP_phospho_site[!all_sites$PSP_phospho_site]=0
fn = paste(pd,"phosphosite_upset.pdf",sep="_")
pdf(fn, height=4, useDingbats = F)
upset(all_sites,sets = c("CPTAC_BRCA","CPTAC_OV","Known_cancer_site","PDB_phospho_site","PSP_phosphosite" ),order.by = "freq")
dev.off()
