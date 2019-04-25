##### phosphosite_summary.R #####
# Kuan-lin Huang @ WashU 2016-2017 Feb
# cross-correlation between different phosphosites

baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/"

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/phospho_summary/")
#source("/Users/khuang/bin/LIB_exp.R")
source("../pan3can_aes.R")
library(biomaRt)
library(UpSetR)
library(VennDiagram)
library(ggrepel)

render = F

get_phosphosites = function(site_genomic_pos){unique(site_genomic_pos[!is.na(site_genomic_pos)])}

cancer_gene_fn = "/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/cgenes_and_druggable.list"
cancer_genes_f = read.table(file=cancer_gene_fn, header=FALSE, stringsAsFactors = F)
cancer_genes = as.vector(t(cancer_genes_f))

# databases
psp_phosphosites = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/HotPho/data_runs/PSP_phosphosite_")
pho_elm_phosphosites = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/Phospho.ELM/phosphoELM_all_2015-04.dump_human_transvar_out.txt")
pdb_phosphosites = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/HotPho/data_runs/PDB_phosphosite_wGpos.maf")
hgnc_fn = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Gene_family/HGNC_protein-coding_gene.txt"
hgnc = read.table(header=TRUE, quote = "", sep="\t", fill = T, file= hgnc_fn)
hgnc_fam = hgnc[,c("symbol","gene_family")]
remove(hgnc)
pps_cancer_wgene = read.table(header=TRUE, quote = "", sep="\t",file="known_cancer_sites_coverage_breakdown.tsv")
# manning_kinome = read.table(header=TRUE, quote = "", sep="\t",file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Gene_family/knome_render_peerj-01-126-s001.txt")
manning_kinome_wgene = read.table(header =T, quote = "",sep = '\t', row.names = NULL,file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Gene_family/knome_render_peerj-01-126-s001_wgene.txt")
#manning_kinome_wgene = manning_kinome_wgene[,-c(12:14)]

pdb_phosphosites = pdb_phosphosites[!is.na(pdb_phosphosites$genomic_pos),]
pho_elm_phosphosites = pho_elm_phosphosites[!is.na(pho_elm_phosphosites) & !duplicated(pho_elm_phosphosites[,5]) & pho_elm_phosphosites[,5] != "././.",]

pho_elm_phosphosites_gpos = as.character(unique(gsub("\\/c.*","",pho_elm_phosphosites[,5])))

cat("# number of unique, mapped PSP sites: ", length(get_phosphosites(XX$genomic_pos)),"\n")
cat("# number of unique, mapped PDB sites: ", length(get_phosphosites(pdb_phosphosites$genomic_pos)),"\n")
cat("# number of unique, mapped Phospho.ELM sites: ", length(get_phosphosites(pho_elm_phosphosites_gpos)),"\n")

# # biomart conversion to gene name: modularize this to a function
# ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
# mapTab = getBM(attributes = c("hgnc_symbol", "uniprotswissprot"), filters = "uniprotswissprot", values = manning_kinome$UniprotID, mart = ensembl, uniqueRows=FALSE)
# dupRows = union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
# if (length(dupRows)==0){mapTabNoDup = mapTab} else {mapTabNoDup = mapTab[-dupRows, ]}
# colnames(manning_kinome)[4]="uniprot_swissprot"
# manning_kinome_wgene = merge(manning_kinome, mapTabNoDup, by="uniprot_swissprot",all.x=T)
# manning_kinome_wgene$hgnc_symbol[is.na(manning_kinome_wgene$hgnc_symbol)] = as.character(manning_kinome_wgene$Name[is.na(manning_kinome_wgene$hgnc_symbol)])
# write.table(manning_kinome_wgene,quote=F, sep = '\t', row.names = FALSE,file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Gene_family/knome_render_peerj-01-126-s001_wgene.txt")

# all sites
HUMAN_file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/BRCA77_unimodal_phosphoproteome-ratio-norm_wGpos_cleaned.txt"
PDX_file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM_WUPCC_shared/WHIM\ Proteomic\ Data/WHIM_phosphoproteome/phosphoproteome-ratio-norm_exp_filtered.txt"
HUMAN_pho = read.table(stringsAsFactors = F,row.names=NULL, header=TRUE, sep="\t", file= HUMAN_file)
PDX_pho = read.table(stringsAsFactors = F,row.names=NULL, header=TRUE, sep="\t", file= PDX_file)
num_unique_sites = length(unique(c(HUMAN_pho$genomic_pos,PDX_pho$genomic_pos)))

cat("# number of unique, mapped BRCA sites: ", length(get_phosphosites(HUMAN_pho$genomic_pos)),"\n")
cat("# number of unique, mapped BRCA_PDX sites: ", length(get_phosphosites(PDX_pho$genomic_pos)),"\n")

# gene count
HUMAN_genes_count = as.data.frame(table(HUMAN_pho$Gene))
PDX_genes_count = as.data.frame(table(PDX_pho$gene))
colnames(HUMAN_genes_count) = c("symbol","cohort_count")
colnames(PDX_genes_count) = c("symbol","cohort_count")
HUMAN_genes_count$data = "HUMAN"
PDX_genes_count$data = "PDX"
pdb_phosphosites_count = as.data.frame(table(pdb_phosphosites$Hugo_Symbol))
colnames(pdb_phosphosites_count) = c("symbol","database_count")
pdb_phosphosites_count$database = "PDB"
phospho_ELM_phosphosites_count = as.data.frame(table(pho_elm_phosphosites$gene))
colnames(phospho_ELM_phosphosites_count) = c("symbol","database_count")
phospho_ELM_phosphosites_count$database = "Phospho.ELM"

merge_HUMAN_count_PDB = merge(HUMAN_genes_count,pdb_phosphosites_count,by="symbol",all=T)
merge_HUMAN_count_ELM = merge(HUMAN_genes_count,phospho_ELM_phosphosites_count,by="symbol",all=T)
merge_PDX_count_PDB = merge(PDX_genes_count,pdb_phosphosites_count,by="symbol",all=T)
merge_PDX_count_ELM = merge(PDX_genes_count,phospho_ELM_phosphosites_count,by="symbol",all=T)

both_count = rbind(merge_HUMAN_count_PDB,merge_PDX_count_PDB,merge_HUMAN_count_ELM, merge_PDX_count_ELM)
both_count = both_count[!is.na(both_count$data) & !is.na(both_count$database),]
both_count = both_count[both_count$database_count < 300,] # remPDXe outlier SRRM2
# > both_count[both_count$database_count> 300,]
# symbol cohort_count data database_count    database
# 16772  SRRM2          328 HUMAN            381 Phospho.ELM
# 21511  SRRM2           79   PDX            381 Phospho.ELM

p = ggplot(data=both_count, aes(x=cohort_count,y=database_count))
p = p + facet_grid(database~data, drop=T)#, space ="free")
p = p + geom_text(aes(label = ifelse((data=="HUMAN" & (cohort_count/database_count > 2) & (cohort_count > 20))
                                     ,as.character(symbol),NA)),size=2, alpha=0.8,hjust=-0.15)            
p = p + geom_point(alpha = 0.1,size=0.5, shape=16)
p = p + geom_abline(slope=1, alpha =0.5)
p = p + theme_bw()
p = p + labs(x="Cohort phosphosite count", x="Database phosphosite count")
p
fn = paste(pd, 'phosphosite_HUMAN.PDX.pdb.count_each_gene_summary.pdf',sep ="_")
ggsave(file=fn, height = 5, width=6, useDingbats=FALSE)

# PDB
pdb_phosphosites_non_unique = gsub("p.","",pdb_phosphosites$amino_acid_change)
pdb_phosphosites_residue = substring(pdb_phosphosites_non_unique, 1, 1)
pdb_phosphosites_residue = pdb_phosphosites_residue[pdb_phosphosites_residue!="H"]
pdb_phosphosites_r_count = table(pdb_phosphosites_residue)
pdb_phosphosites_r_count

# phospho.ELM
pho_elm_phosphosites_sites = gsub(".*:","",pho_elm_phosphosites$input)
pho_elm_phosphosites_residue = substring(pho_elm_phosphosites_sites, 1, 1)
pho_elm_phosphosites_residue = pho_elm_phosphosites_residue[pho_elm_phosphosites_residue!="H"]
pho_elm_phosphosites_r_count = table(pho_elm_phosphosites_residue)
pho_elm_phosphosites_r_count

# HUMAN
non_unique = gsub(".*:","",HUMAN_pho$Gene.site)
HUMAN_residue = substring(non_unique, 1, 1)
HUMAN = table(HUMAN_residue)
HUMAN

# PDX
non_unique = gsub(".*:","",PDX_pho$gene.site)
PDX_residue = substring(non_unique, 1, 1)
PDX = table(PDX_residue)
PDX

# both
residue_table = rbind(HUMAN,PDX,pdb_phosphosites_r_count,pho_elm_phosphosites_r_count)
residue_table_m = melt(residue_table)
colnames(residue_table_m) = c("Cancer","Residue","Number_of_Phosphosites")
residue_table_m$Cancer = c("BRCA","BRCA PDX", "PDB","Phospho.ELM", 
                           "BRCA","BRCA PDX", "PDB","Phospho.ELM",
                           "BRCA","BRCA PDX", "PDB","Phospho.ELM")

p = ggplot(residue_table_m,aes(x=Cancer, y=Number_of_Phosphosites, fill=Residue))
#p = p + facet_grid(Cancer~.,drop=T,space="free",scale="free")
p = p + geom_bar(stat="identity") + theme_bw() #+ theme_nogrid()
p = p + labs(x = "Dataset", y="# of phosphosites")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
#p = p + theme(legend.position="bottom")
#p = p + coord_flip()
#p = p + geom_text(aes(label=Number_of_Phosphosites), vjust=-0.25)
p
fn = paste(pd, 'phosphosite_HUMAN.PDX_sty_summary.pdf',sep ="_")
ggsave(file=fn, height = 4, width=3, useDingbats=FALSE)

# ##### upsetR between all the universe #####
# all_sites = data.frame(stringsAsFactors =F, unique(as.character(c(HUMAN_pho$genomic_pos,PDX_pho$genomic_pos,pps_cancer_wgene$genomic_pos,pdb_phosphosites_gpos,pho_elm_phosphosites_gpos))))
# colnames(all_sites)[1] = "Phosphosite"
# 
# all_sites$CPTAC_HUMAN = all_sites$Phosphosite %in% HUMAN_pho$genomic_pos
# all_sites$CPTAC_PDX = all_sites$Phosphosite %in% PDX_pho$genomic_pos
# all_sites$Known_cancer_site = all_sites$Phosphosite %in% pps_cancer_wgene$genomic_pos
# all_sites$PDB_phospho_site = all_sites$Phosphosite %in% pdb_phosphosites_gpos
# all_sites$Phospho_ELM_site = all_sites$Phosphosite %in% pho_elm_phosphosites_gpos
# all_sites = all_sites[!is.na(all_sites$Phosphosite) & all_sites$Phosphosite != ".",]
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
# upset(all_sites,sets = c("BRCA","BRCA PDX","Known_cancer_site","PDB_phospho_site","Phospho_ELM_site" ),order.by = "freq")
# dev.off()

# venn diagram

phosphosite_list = list(BRCA = get_phosphosites(HUMAN_pho$genomic_pos), #cancer_gene = get_phosphosites(pps_cancer_wgene$genomic_pos),
                  PDB = get_phosphosites(pdb_phosphosites$genomic_pos),Phospho.ELM = get_phosphosites(pho_elm_phosphosites_gpos))
venn.plot = venn.diagram(
  x = phosphosite_list, 
  filename = NULL,
  lwd = 4,
  #fill = c("royalblue1", "orange1"),
  alpha = 0.2,
  scaled =T,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  #cat.col = c("royalblue1", "orange1"),
  cat.cex = 1,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  #cat.dist = c(0.16, 0.16),
  #cat.pos = c(240, 120), 
  cat.default.pos = c("text"), 
  #height=2200, width=2200, 
  margin=0.2
)

fn = paste(pd, '_phosphosite_overlap.pdf', sep = "")
pdf(file=fn)
grid.draw(venn.plot)
dev.off()

##### analyze in the context of manning kinome #####
# BRCA_sites_data = cbind(BRCA_sites,BRCA_genes)
HUMAN_pho_count = HUMAN_pho[,c("genomic_pos","Gene","Gene.site")]
HUMAN_pho_count$Status = "Novel"
HUMAN_pho_count$Status[HUMAN_pho_count$genomic_pos %in% 
                  c(pdb_phosphosites$genomic_pos, pho_elm_phosphosites_gpos)] = "In_Databases"

BRCA_genes_count = as.data.frame(table(HUMAN_pho_count$Gene,HUMAN_pho_count$Status))
colnames(BRCA_genes_count) = c("hgnc_symbol","status","unique_phosphosite_count")
BRCA_genes_count_c = dcast(data=BRCA_genes_count,hgnc_symbol~status)
BRCA_genes_count_c$unique_phosphosite_count = BRCA_genes_count_c$In_Databases + BRCA_genes_count_c$Novel
# add novel sites counts
manning_kinome_wgene_wcount = merge(manning_kinome_wgene, BRCA_genes_count_c, by="hgnc_symbol", all.x=T)
tn = "manning_kinome_wgene_wphositecount.tsv"
write.table(manning_kinome_wgene_wcount, row.names = F, quote=F, sep = '\t', file = tn)

cat ("# Phosphosites covering ",sum(!is.na(manning_kinome_wgene_wcount$unique_phosphosite_count))," out of ",nrow(manning_kinome_wgene_wcount)," kinases (manning)\n")

manning_kinome_wgene_wcount = manning_kinome_wgene_wcount[!is.na(manning_kinome_wgene_wcount$unique_phosphosite_count),]
p = ggplot(manning_kinome_wgene_wcount,aes(x=GroupName, y=unique_phosphosite_count, color=GroupName))
#p = p + facet_grid(.~GroupName, drop=T,space="free",scale="free")
p = p + geom_jitter(data = subset(manning_kinome_wgene_wcount,unique_phosphosite_count<=10),alpha=0.2,stroke=0) 
#p = p + geom_label(aes(label=ifelse(unique_phosphosite_count >= 10,as.character(hgnc_symbol),NA)),alpha=0.6, position=position_jitter(width=0.5, height=0))
p = p + geom_point(data = subset(manning_kinome_wgene_wcount,unique_phosphosite_count>10),alpha=0.2,stroke=0) + theme_bw() #+ theme_nogrid()
p = p + geom_text_repel(aes(label=ifelse(unique_phosphosite_count > 10,as.character(hgnc_symbol),NA)),alpha=0.8)
p = p + labs(x = "Kinase group", y="# of phosphosites") #+ scale_y_log10(breaks = c(0, 2, 10, 20, 40, 80))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + theme(legend.position="none")
p
fn = paste(pd, 'phosphosite_kinome_group_gene_summary.pdf',sep ="_")
ggsave(file=fn, width = 8, height=5, useDingbats=FALSE)

manning_kinome_wgene_wcount_10 = manning_kinome_wgene_wcount[!is.na(manning_kinome_wgene_wcount$Novel) & manning_kinome_wgene_wcount$Novel>=10,]
manning_kinome_wgene_wcount_10_m = melt(manning_kinome_wgene_wcount_10[,c("hgnc_symbol","GroupName","In_Databases","Novel")], id.var=c("hgnc_symbol","GroupName"))
manning_kinome_wgene_wcount_10_m$variable = factor(manning_kinome_wgene_wcount_10_m$variable, levels=c("Novel","In_Databases"))

# analyze by family
p = ggplot(manning_kinome_wgene_wcount_10_m,aes(x=hgnc_symbol, y=value, fill=variable))
p = p + facet_grid(.~GroupName, drop=T,space="free",scale="free")
p = p + geom_bar(stat="identity") + theme_bw() #+ theme_nogrid()
p = p + labs(x = "Kinase", y="Number of phosphosites") #+ scale_y_log10(breaks = c(0, 2, 10, 20, 40, 80))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position="bottom")
p
fn = paste(pd, 'phosphosite_kinome_fam_novel_w10_summary.pdf',sep ="_")
ggsave(file=fn, width = 8, height=5,useDingbats=FALSE)

# create file for kinome render
if (render){
  cat("color 0.6 0 0\n\n") #RGB for text/circle color
  for (i in 1:nrow(manning_kinome_wgene_wcount)){
    name = as.character(manning_kinome_wgene_wcount[i,]$Name)
    hgnc_name = as.character(manning_kinome_wgene_wcount[i,]$hgnc_symbol)
    phosite_count = 0
    if (!is.na(manning_kinome_wgene_wcount[i,]$unique_phosphosite_count)) { 
      phosite_count = manning_kinome_wgene_wcount[i,]$unique_phosphosite_count
      cat("at ",name,"\n", sep="")
      cat("scale 10\n", sep="")
      cat("circle-lined\n\n", sep="")
      scale = phosite_count
      
      if (scale > 10){
        if (scale > 30) {scale=30}
        cat("at ",name,"\n", sep="")
        cat("scale ",scale,"\n", sep="")
        cat("text  ",hgnc_name,":",phosite_count,"\n\n", sep="")
      }
      
    }
  }
}


##### pro-pho correlation #####
# brca_pho_collapsed_f = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt"
# brca_pho_collapsed = read.table(row.names=1, header=TRUE, sep="\t", file= brca_pho_collapsed_f)
# 
# brca_pro_f = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt"
# brca_pro = read.table(row.names=1, header=TRUE, sep="\t", file=brca_pro_f)

# # correlation between collapsed and non-collapsed
# brca_pho$gene_site = row.names(brca_pho)
# brca_pho_m = melt(brca_pho,id="gene_site")
# brca_pho_collapsed$gene = row.names(brca_pho_collapsed)
# brca_pho_collapsed_m = melt(brca_pho_collapsed,id="gene")
# brca_pho_m$gene = gsub(":.*","",brca_pho_m$gene_site)
# 
# brca_pro$gene = row.names(brca_pro)
# brca_pro_m = melt(brca_pro,id="gene")
# colnames(brca_pro_m) = c("gene" ,    "sample", "protein_level")
# 
# colnames(brca_pho_m) = c("gene_site" ,"sample" , "phosphosite_level"  ,   "gene")
# colnames(brca_pho_collapsed_m) = c("gene" ,"sample" , "phosphoprotein_level")
# brca_pho_collapsed_m_merge = merge(brca_pho_m,brca_pho_collapsed_m,by=c("sample","gene"))
# brca_pro_pho_collapsed_m_merge = merge(brca_pro_m,brca_pho_collapsed_m_merge,by=c("sample","gene"))
# ##### investigate expression vs. variance #####
# RTK_file = read.table("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/reference_files/RTKs_list.txt")
# RTK = as.vector(t(RTK_file))
# 
# BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))
# gene = row.names(BRCA_RNA)
# gene_exp = rowMeans(BRCA_RNA,na.rm=T)
# RNA_gene_exp = cbind(gene, gene_exp)
# 
# site = row.names(brca_pho_10NA)
# gene = gsub(":.*","",site)
# site_SD = rowSds(as.matrix(brca_pho_10NA),na.rm=T)
# pho_site_gene_SD = cbind(gene, site, site_SD)
# 
# RNA_gene_exp_pho_site_gene_SD = merge(pho_site_gene_SD, RNA_gene_exp, by="gene", all.x=T)
# 
# RNA_gene_exp_pho_site_gene_SD$gene_exp = as.numeric(as.character(RNA_gene_exp_pho_site_gene_SD$gene_exp))
# RNA_gene_exp_pho_site_gene_SD$site_SD = as.numeric(as.character(RNA_gene_exp_pho_site_gene_SD$site_SD))
# 
# RNA_gene_exp_pho_site_gene_SD_RTK = RNA_gene_exp_pho_site_gene_SD[RNA_gene_exp_pho_site_gene_SD$gene %in% RTK,] 
# 
# p = ggplot(RNA_gene_exp_pho_site_gene_SD_RTK,aes(x=gene_exp, y=site_SD))
# #p = p + facet_grid(self~Cancer,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
# #p = p + geom_point(alpha=0.5)
# #p = p + geom_text(aes(label= ifelse(gene_exp > 10 & site_SD > 2.6, as.character(site), NA)),size=2,alpha=0.5)
# p = p + geom_text(aes(label= ifelse(gene %in% RTK, as.character(site), NA), color = gene),size=2,alpha=0.5)
# p = p + theme_bw() #+ theme_nogrid()
# #p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
# p = p + theme(legend.position="none")
# p = p + labs(x = "Average gene expression (log2)", y="Phosphosite expression SD")
# p = p + expand_limits(x = 0, y = 0) + xlim(0,15)
# p
# fn = paste(pd, 'phosphositeSD_vs_geneExp_RTK.pdf',sep ="_")
# ggsave(file=fn, height=5, width=5, useDingbats=FALSE)
# 
# ## violin plot of expression
# BRCA_Clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", file=paste(baseD,"BRCA/BRCA_clinical_summary.txt",sep=""))
# brca_pho_10NA_RTK = brca_pho_10NA[gene %in% RTK,]
# brca_pho_10NA_RTK_t = t(brca_pho_10NA_RTK)
# brca_pho_10NA_RTK_t_m = melt(brca_pho_10NA_RTK_t)
# colnames(brca_pho_10NA_RTK_t_m) = c("Sample","Site","Expression")
# BRCA_Clin_t = as.data.frame(t(BRCA_Clin))
# BRCA_Clin_t$Sample = row.names(BRCA_Clin_t)
# brca_pho_10NA_RTK_clin = merge(brca_pho_10NA_RTK_t_m, BRCA_Clin_t, by="Sample")
# brca_pho_10NA_RTK_clin$Gene = gsub(":.*","",brca_pho_10NA_RTK_clin$Site)
# 
# # plotting
# get.clinical.scale = function() {
#   # Set1 colors
#   colors = c(NA, "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
#   # use Perou's intrinsic subtype colors instead
#   colors = c(NA, "#101010", NA, "#636363", "red", "pink", "#ffeda0", "#deebf7", "#3182bd", "#FFFF33", "#A65628", "#F781BF", "#999999") #positive is dark grey       
#   
#   color.names = c("wt","mut","negative", "positive", "Basal", "Her2", "CLDN low", "LumB", "LumA") # add lumA and check color
#   #color.names = c("negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB")
#   names(colors) = color.names
#   clinical.color.scale = scale_color_manual(name="Intrinsic Subtype", values=colors)
#   
#   return(clinical.color.scale)
# }
# 
# p = ggplot(data=brca_pho_10NA_RTK_clin)
# #p = p + facet_wrap(Gene~Site, nrow=5)
# p = p + facet_grid(.~Gene,scale="free", space="free")
# #p = p + geom_boxplot(aes(x=Site, y=Expression, fill=NULL),alpha=0.1, outlier.shape = NA) 
# p = p + geom_jitter(aes(x=Site, y=Expression, color=pam50), size=1, alpha=0.4) #+ geom_point(aes(x=Status, y=value)) 
# #p = p + geom_text(aes(x=Species, y=Pho, label = ifelse(outlier,as.character(Sample),NA)),size=2)
# p = p + theme_nogrid() + guides(fill=FALSE) 
# p = p + labs(x = "", y = "Phosphosite expression")
# p = p + get.clinical.scale()
# p = p + theme(text = element_text(colour="black", size=16), axis.ticks.x = element_blank(),
#               axis.text.x = element_text(colour="black", size=5,angle=90,vjust=0.5),#element_blank(), axis.title.x = element_blank(),  
#               axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
# p = p + theme(legend.position="bottom") + ylim(-7.2,7.2)
# p
# fn = paste(pd, "BRCA_phosite_RTK_by_subtype.pdf", sep="_")
# ggsave(file=fn, width=12, height=6,useDingbats=FALSE)