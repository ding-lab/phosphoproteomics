##### find_plot_pQTL.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run pQTL analysis for 3 cancer types and plot the result

##### dependencies #####
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

setwd(paste(baseD,"pan3can_analysis/mutation_impact/", sep=""))
source("/Users/khuang/bin/LIB_exp.R")
source("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/mutation_impact/mutation_impact.R")

pd = paste(pd,"all", sep="_")
# system("mkdir logs")
# logFile = paste("logs/", date, "_mutation_impact_analysis.log", sep="")
# sink(file=logFile)
#sink(file=NULL)
kinaseList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/2014-12-05_CPTAC_Kinase.MATRIX.v3b5_sheet1_genes.list', header=FALSE, stringsAsFactors = F)
drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv_list.txt_list', header=FALSE, stringsAsFactors = F)
cancer_genes127 = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/Kandoth_127SMG_list_brca_refseq_ncbi.txt.list', header=FALSE, stringsAsFactors = F)

kinome = as.vector(t(kinaseList))
druggable = as.vector(t(drugList))
druggable = as.vector(t(cancer_genes127))
SMGs = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C","TP53", "NF1", "KRAS", "BRCA1", "BRCA2", "CDK12",
         "TP53","KRAS","APC","PIK3CA","SMAD4")

##### BRCA #####
### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep=""))
BRCA_mut_g = BRCA_mut[row.names(BRCA_mut) %in% SMGs,]

### CNV ###
BRCA_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_CNV_formatted_normalized.txt",sep=""))
BRCA_CNV_c = BRCA_CNV[row.names(BRCA_CNV) %in% SMGs,]
BRCA_CNV_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_CNV_c,name="BRCA CNV")

### RNA ###
BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))
BRCA_RNA_c = BRCA_RNA[row.names(BRCA_RNA) %in% SMGs,]
BRCA_RNA_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_RNA_c,name="BRCA RNA")

### Proteome ###
BRCA_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
BRCA_Pro_c = BRCA_Pro[row.names(BRCA_Pro) %in% SMGs,]
BRCA_Pro_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_Pro_c,name="BRCA Proteome")

### Phosphoproteome ###
BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
BRCA_Pho_c = BRCA_Pho[row.names(BRCA_Pho) %in% SMGs,]
BRCA_Pho_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_Pho_c,name="BRCA Phosphoproteome")

### all levels ###
# do this later, may be more benefitial to do the all outlier table and summarize that instead

##### OV #####
### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_SOMATIC_formatted.txt",sep=""))
OV_mut_g = OV_mut[row.names(OV_mut) %in% SMGs,]

### CNV ###
OV_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_CNV_formatted_normalized.txt",sep=""))
OV_CNV_c = OV_CNV[row.names(OV_CNV) %in% SMGs,]
OV_CNV_diff_exp = find_diff_exp(OV_mut_g,OV_CNV_c,name="OV CNV")

### RNA ###
OV_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_mRNA_formatted_normalized.txt",sep=""))
OV_RNA_c = OV_RNA[row.names(OV_RNA) %in% SMGs,]
OV_RNA_diff_exp = find_diff_exp(OV_mut_g,OV_RNA_c,name="OV RNA")

### Proteome ###
OV_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))
OV_Pro_c = OV_Pro[row.names(OV_Pro) %in% SMGs,]
OV_Pro_diff_exp = find_diff_exp(OV_mut_g,OV_Pro_c,name="OV Proteome")

### Phosphoproteome ###
OV_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))
OV_Pho_c = OV_Pho[row.names(OV_Pho) %in% SMGs,]
OV_Pho_diff_exp = find_diff_exp(OV_mut_g,OV_Pho_c,name="OV Phosphoproteome")

##### CRC #####
### Mutation matrix ###
CRC_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_SOMATIC_formatted.txt",sep=""))
CRC_mut_g = CRC_mut[row.names(CRC_mut) %in% SMGs,]

### CNV ###
CRC_CNV = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_CNV_formatted_normalized.txt",sep=""))
CRC_CNV_c = CRC_CNV[row.names(CRC_CNV) %in% SMGs,]
CRC_CNV_diff_exp = find_diff_exp(CRC_mut_g,CRC_CNV_c,name="CRC CNV")

### RNA ###
CRC_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_mRNA_formatted_normalized.txt",sep=""))
CRC_RNA_c = CRC_RNA[row.names(CRC_RNA) %in% SMGs,]
CRC_RNA_diff_exp = find_diff_exp(CRC_mut_g,CRC_RNA_c,name="CRC RNA")

### Proteome ###
CRC_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_PRO_formatted_normalized.txt",sep=""))
CRC_Pro_c = CRC_Pro[row.names(CRC_Pro) %in% SMGs,]
CRC_Pro_diff_exp = find_diff_exp(CRC_mut_g,CRC_Pro_c,name="CRC Proteome")

##### merge everything and plot heatmap #####

### fold change
BRCA_CNV_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_CNV_diff_exp$fold_change)))
BRCA_CNV_diff_exp_fc_m$Cancer = "BRCA"
BRCA_CNV_diff_exp_fc_m$Data = "CNV"
BRCA_RNA_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_RNA_diff_exp$fold_change)))
BRCA_RNA_diff_exp_fc_m$Cancer = "BRCA"
BRCA_RNA_diff_exp_fc_m$Data = "RNA"
BRCA_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_Pro_diff_exp$fold_change)))
BRCA_Pro_diff_exp_fc_m$Cancer = "BRCA"
BRCA_Pro_diff_exp_fc_m$Data = "Protein"
BRCA_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$fold_change)))
BRCA_Pho_diff_exp_fc_m$Cancer = "BRCA"
BRCA_Pho_diff_exp_fc_m$Data = "Phosphoprotein"

OV_CNV_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_CNV_diff_exp$fold_change)))
OV_CNV_diff_exp_fc_m$Cancer = "OV"
OV_CNV_diff_exp_fc_m$Data = "CNV"
OV_RNA_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_RNA_diff_exp$fold_change)))
OV_RNA_diff_exp_fc_m$Cancer = "OV"
OV_RNA_diff_exp_fc_m$Data = "RNA"
OV_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_Pro_diff_exp$fold_change)))
OV_Pro_diff_exp_fc_m$Cancer = "OV"
OV_Pro_diff_exp_fc_m$Data = "Protein"
OV_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_Pho_diff_exp$fold_change)))
OV_Pho_diff_exp_fc_m$Cancer = "OV"
OV_Pho_diff_exp_fc_m$Data = "Phosphoprotein"

CRC_CNV_diff_exp_fc_m = melt(as.matrix(do.call(cbind,CRC_CNV_diff_exp$fold_change)))
CRC_CNV_diff_exp_fc_m$Cancer = "CRC"
CRC_CNV_diff_exp_fc_m$Data = "CNV"
CRC_RNA_diff_exp_fc_m = melt(as.matrix(do.call(cbind,CRC_RNA_diff_exp$fold_change)))
CRC_RNA_diff_exp_fc_m$Cancer = "CRC"
CRC_RNA_diff_exp_fc_m$Data = "RNA"
CRC_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,CRC_Pro_diff_exp$fold_change)))
CRC_Pro_diff_exp_fc_m$Cancer = "CRC"
CRC_Pro_diff_exp_fc_m$Data = "Protein"

# do this by function call
diff_exp_matrix_fc_list = list(BRCA_CNV_diff_exp_fc_m,BRCA_RNA_diff_exp_fc_m,BRCA_Pro_diff_exp_fc_m,BRCA_Pho_diff_exp_fc_m,
                            OV_CNV_diff_exp_fc_m,OV_RNA_diff_exp_fc_m,OV_Pro_diff_exp_fc_m,OV_Pho_diff_exp_fc_m,
                            CRC_CNV_diff_exp_fc_m,CRC_RNA_diff_exp_fc_m,CRC_Pro_diff_exp_fc_m) 

all_diff_exp_fc_matrix = do.call(rbind,diff_exp_matrix_fc_list)
colnames(all_diff_exp_fc_matrix)[1:3] = c("Gene","Mutated_gene","Fold_change")

### FDR
BRCA_CNV_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_CNV_diff_exp$t_fdr)))
BRCA_CNV_diff_exp_fdr_m$Cancer = "BRCA"
BRCA_CNV_diff_exp_fdr_m$Data = "CNV"
BRCA_RNA_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_RNA_diff_exp$t_fdr)))
BRCA_RNA_diff_exp_fdr_m$Cancer = "BRCA"
BRCA_RNA_diff_exp_fdr_m$Data = "RNA"
BRCA_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_Pro_diff_exp$t_fdr)))
BRCA_Pro_diff_exp_fdr_m$Cancer = "BRCA"
BRCA_Pro_diff_exp_fdr_m$Data = "Protein"
BRCA_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$t_fdr)))
BRCA_Pho_diff_exp_fdr_m$Cancer = "BRCA"
BRCA_Pho_diff_exp_fdr_m$Data = "Phosphoprotein"

OV_CNV_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_CNV_diff_exp$t_fdr)))
OV_CNV_diff_exp_fdr_m$Cancer = "OV"
OV_CNV_diff_exp_fdr_m$Data = "CNV"
OV_RNA_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_RNA_diff_exp$t_fdr)))
OV_RNA_diff_exp_fdr_m$Cancer = "OV"
OV_RNA_diff_exp_fdr_m$Data = "RNA"
OV_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_Pro_diff_exp$t_fdr)))
OV_Pro_diff_exp_fdr_m$Cancer = "OV"
OV_Pro_diff_exp_fdr_m$Data = "Protein"
OV_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_Pho_diff_exp$t_fdr)))
OV_Pho_diff_exp_fdr_m$Cancer = "OV"
OV_Pho_diff_exp_fdr_m$Data = "Phosphoprotein"

CRC_CNV_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,CRC_CNV_diff_exp$t_fdr)))
CRC_CNV_diff_exp_fdr_m$Cancer = "CRC"
CRC_CNV_diff_exp_fdr_m$Data = "CNV"
CRC_RNA_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,CRC_RNA_diff_exp$t_fdr)))
CRC_RNA_diff_exp_fdr_m$Cancer = "CRC"
CRC_RNA_diff_exp_fdr_m$Data = "RNA"
CRC_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,CRC_Pro_diff_exp$t_fdr)))
CRC_Pro_diff_exp_fdr_m$Cancer = "CRC"
CRC_Pro_diff_exp_fdr_m$Data = "Protein"

diff_exp_matrix_fdr_list = list(BRCA_CNV_diff_exp_fdr_m,BRCA_RNA_diff_exp_fdr_m,BRCA_Pro_diff_exp_fdr_m,BRCA_Pho_diff_exp_fdr_m,
                               OV_CNV_diff_exp_fdr_m,OV_RNA_diff_exp_fdr_m,OV_Pro_diff_exp_fdr_m,OV_Pho_diff_exp_fdr_m,
                               CRC_CNV_diff_exp_fdr_m,CRC_RNA_diff_exp_fdr_m,CRC_Pro_diff_exp_fdr_m) 

all_diff_exp_fdr_matrix = do.call(rbind,diff_exp_matrix_fdr_list)
colnames(all_diff_exp_fdr_matrix)[1:3] = c("Gene","Mutated_gene","FDR")

fdr_fc = merge(all_diff_exp_fc_matrix,all_diff_exp_fdr_matrix, by=c("Gene","Mutated_gene","Cancer","Data"))
fdr_fc_cis = fdr_fc[fdr_fc$Gene == fdr_fc$Mutated_gene,]

fdr_fc_cis$Fold_change = as.numeric(as.character(fdr_fc_cis$Fold_change))

fdr_fc_cis$Data = factor(fdr_fc_cis$Data, levels = c("Phosphoprotein","Protein","RNA","CNV"))
fdr_fc_cis$Cancer = factor(fdr_fc_cis$Cancer, levels = c("BRCA","OV","CRC"))
### plot ###
fn = paste(pd, 'all_cis_effects.pdf',sep ="_")
fdr.colors=c("NA", "#000000")
min_d = min(fdr_fc_cis$Fold_change), na.rm=T)
max_d = max(fdr_fc_cis$Fold_change), na.rm=T)
bound = max(c(max_d, -min_d))
fdr_fc_cis$sig = as.numeric(as.character(fdr_fc_cis$FDR)) <= 0.1

p = ggplot(data=fdr_fc_cis)
p = p + facet_grid(Cancer ~ ., drop=T,scales = "free_y")
#p = p + facet_grid(Cancer ~ ., drop=T,scales = "free", space = "free")
p = p + geom_tile(aes(x=Gene, y=Data, fill=Fold_change), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-bound,bound))
p = p + geom_tile(aes(x=Gene, y=Data, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
p = p + labs(x="Mutated Gene", y="Fold Change") + theme_bw() + theme_nogrid() +
  theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
p = p + coord_fixed()
p
ggsave(file=fn, useDingbats=FALSE)
#ggsave(file=fn, width=20, limitsize=FALSE, height=20, useDingbats=FALSE)
