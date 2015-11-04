##### find_plot_pQTL.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run pQTL analysis for 3 cancer types and plot the result

##### dependencies #####
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

setwd(paste(baseD,"pan3can_analysis/mutation_impact/", sep=""))
source("/Users/khuang/bin/LIB_exp.R")
source("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/mutation_impact/mutation_impact.R")

pd = "figures/2015-10-17/2015-10-17_KH_all"
# system("mkdir logs")
# logFile = paste("logs/", date, "_mutation_impact_analysis.log", sep="")
# sink(file=logFile)
#sink(file=NULL)
kinaseList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/2014-12-05_CPTAC_Kinase.MATRIX.v3b5_sheet1_genes.list', header=FALSE, stringsAsFactors = F)
drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv_list.txt_list', header=FALSE, stringsAsFactors = F)
cancer_genes127 = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/Kandoth_127SMG_list_brca_refseq_ncbi.txt.list', header=FALSE, stringsAsFactors = F)

kinome = as.vector(t(kinaseList))
druggable = as.vector(t(drugList))
cg127 = as.vector(t(cancer_genes127))
SMGs = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C","TP53", "NF1", "KRAS", "BRCA1", "BRCA2", "CDK12",
         "TP53","KRAS","APC","PIK3CA","SMAD4")

##### BRCA #####
### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep=""))
brcaGenes = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C")
BRCA_mut_g = BRCA_mut[row.names(BRCA_mut) %in% SMGs,]

### Proteome ###
BRCA_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
BRCA_Pro_c = BRCA_Pro[row.names(BRCA_Pro) %in% cg127,]
BRCA_Pro_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_Pro,name="BRCA Proteome")

### Phosphoproteome ###
BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
BRCA_Pho_c = BRCA_Pho[row.names(BRCA_Pho) %in% cg127,]
BRCA_Pho_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_Pho,name="BRCA Phosphoproteome")

### all levels ###
# do this later, may be more benefitial to do the all outlier table and summarize that instead

##### OV #####
### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_SOMATIC_formatted.txt",sep=""))
ovGenes = c("TP53", "NF1", "KRAS", "BRCA1", "BRCA2", "CDK12")
OV_mut_g = OV_mut[row.names(OV_mut) %in% SMGs,]

### Proteome ###
OV_JHU_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_PRO_formatted_normalized.txt",sep=""))
OV_JHU_Pro_c = OV_JHU_Pro[row.names(OV_JHU_Pro) %in% cg127,]
OV_JHU_Pro_diff_exp = find_diff_exp(OV_mut_g,OV_JHU_Pro,name="OV JHU Proteome")

OV_PNNL_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))
OV_PNNL_Pro_c = OV_PNNL_Pro[row.names(OV_PNNL_Pro) %in% cg127,]
OV_PNNL_Pro_diff_exp = find_diff_exp(OV_mut_g,OV_PNNL_Pro,name="OV PNNL Proteome")

### Phosphoproteome ###
OV_PNNL_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))
OV_PNNL_Pho_c = OV_PNNL_Pho[row.names(OV_PNNL_Pho) %in% cg127,]
OV_PNNL_Pho_diff_exp = find_diff_exp(OV_mut_g,OV_PNNL_Pho,name="OV PNNL Phosphoproteome")

### Glycoproteome ###
OV_JHU_Gly = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_GLY_formatted_normalized.txt",sep=""))
OV_JHU_Gly_c = OV_JHU_Gly[row.names(OV_JHU_Gly) %in% cg127,]
OV_JHU_Gly_diff_exp = find_diff_exp(OV_mut_g,OV_JHU_Gly,name="OV JHU Glycoproteome")

### merging the two proteome? ###
### all levels ###

##### CRC #####
### Mutation matrix ###
CRC_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_SOMATIC_formatted.txt",sep=""))
crcGenes = c("TP53","KRAS","APC","PIK3CA","SMAD4")
CRC_mut_g = CRC_mut[row.names(CRC_mut) %in% SMGs,]
### Proteome ###
CRC_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_PRO_formatted_normalized.txt",sep=""))
CRC_Pro_c = CRC_Pro[row.names(CRC_Pro) %in% cg127,]
CRC_Pro_diff_exp = find_diff_exp(CRC_mut_g,CRC_Pro,name="CRC Proteome") # will likely need to change to non-parametric test


##### merge everything and plot #####

### fold change
BRCA_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_Pro_diff_exp$fold_change)))
colnames(BRCA_Pro_diff_exp_fc_m)[3] = "BRCA_PRO"
BRCA_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$fold_change)))
colnames(BRCA_Pho_diff_exp_fc_m)[3] = "BRCA_PHO"
OV_JHU_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_JHU_Pro_diff_exp$fold_change)))
colnames(OV_JHU_Pro_diff_exp_fc_m)[3] = "OV_JHU_PRO"
OV_PNNL_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pro_diff_exp$fold_change)))
colnames(OV_PNNL_Pro_diff_exp_fc_m)[3] = "OV_PNNL_PRO"
OV_JHU_Gly_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_JHU_Gly_diff_exp$fold_change)))
colnames(OV_JHU_Gly_diff_exp_fc_m)[3] = "OV_JHU_GLY"
OV_PNNL_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pho_diff_exp$fold_change)))
colnames(OV_PNNL_Pho_diff_exp_fc_m)[3] = "OV_PNNL_PHO"
CRC_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,CRC_Pro_diff_exp$fold_change)))
colnames(CRC_Pro_diff_exp_fc_m)[3] = "CRC_PRO"

fc = merge(BRCA_Pro_diff_exp_fc_m, BRCA_Pho_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
fc = merge(fc, OV_JHU_Pro_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
fc = merge(fc, OV_PNNL_Pro_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
fc = merge(fc, OV_JHU_Gly_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
fc = merge(fc, OV_PNNL_Pho_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
fc = merge(fc, CRC_Pro_diff_exp_fc_m, by = c("Var1","Var2"), all=T)

fc_m = melt(fc, id.var=c("Var1","Var2"))
colnames(fc_m)[4] = "FC"
## fdr
BRCA_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_Pro_diff_exp$t_fdr)))
colnames(BRCA_Pro_diff_exp_fdr_m)[3] = "BRCA_PRO"
BRCA_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$t_fdr)))
colnames(BRCA_Pho_diff_exp_fdr_m)[3] = "BRCA_PHO"
OV_JHU_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_JHU_Pro_diff_exp$t_fdr)))
colnames(OV_JHU_Pro_diff_exp_fdr_m)[3] = "OV_JHU_PRO"
OV_PNNL_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pro_diff_exp$t_fdr)))
colnames(OV_PNNL_Pro_diff_exp_fdr_m)[3] = "OV_PNNL_PRO"
OV_JHU_Gly_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_JHU_Gly_diff_exp$t_fdr)))
colnames(OV_JHU_Gly_diff_exp_fdr_m)[3] = "OV_JHU_GLY"
OV_PNNL_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pho_diff_exp$t_fdr)))
colnames(OV_PNNL_Pho_diff_exp_fdr_m)[3] = "OV_PNNL_PHO"
CRC_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,CRC_Pro_diff_exp$t_fdr)))
colnames(CRC_Pro_diff_exp_fdr_m)[3] = "CRC_PRO"

fdr = merge(BRCA_Pro_diff_exp_fdr_m, BRCA_Pho_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
fdr = merge(fdr, OV_JHU_Pro_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
fdr = merge(fdr, OV_PNNL_Pro_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
fdr = merge(fdr, OV_JHU_Gly_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
fdr = merge(fdr, OV_PNNL_Pho_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
fdr = merge(fdr, CRC_Pro_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)

fdr_m = melt(fdr, id.var=c("Var1","Var2"))
colnames(fdr_m)[4] = "FDR"

fc_fdr = merge(fc_m, fdr_m, by=c("Var1","Var2","variable"))

## plot
# fn = paste(pd, 'merged_diff_exp_127cgenes.pdf',sep ="_")
fdr.colors=c("NA", "#000000")
min_d = min(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
max_d = max(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
bound = max(c(max_d, -min_d))
fc_fdr$sig = as.numeric(as.character(fc_fdr$FDR)) <= 0.1
# 
# p = ggplot(data=fc_fdr)
# p = p + facet_grid(.~variable)
# p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-bound,bound))
# p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
# p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
# p = p + coord_fixed()
# p
# ggsave(file=fn, width=20, limitsize=FALSE, height=20, useDingbats=FALSE)


### plot selected markers ###
fn = paste(pd, 'merged_diff_exp_127_cgenes_fdr0.001in2.pdf',sep ="_")

for (i in 3:9){
  fdr[,i] = as.numeric(as.character(fdr[,i]))
}

markers = unique(fdr[rowSums(fdr[,c(3:9)]<0.001, na.rm=T) >=2,]$Var1)
fc_fdr_s = fc_fdr[fc_fdr$Var1 %in% as.character(markers),]
#fc_fdr_s2 = fc_fdr_s[fc_fdr_s$FDR<=0.2,]
min_d = min(as.numeric(as.character(fc_fdr_s$FC)), na.rm=T)
max_d = max(as.numeric(as.character(fc_fdr_s$FC)), na.rm=T)
bound = max(c(max_d, -min_d))

p = ggplot(data=fc_fdr_s)
p = p + facet_grid(.~variable)
p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-bound,bound))
p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
p = p + coord_fixed()
p
ggsave(file=fn, width=20, limitsize=FALSE, height=20, useDingbats=FALSE)
