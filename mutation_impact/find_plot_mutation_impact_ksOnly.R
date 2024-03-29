##### find_plot_pQTL.R #####
# Kuan-lin Huang @ WashU 2017 Feb
# run pQTL analysis for 3 cancer types and plot the result

##### dependencies #####
baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

setwd(paste(baseD,"pan3can_analysis/mutation_impact/", sep=""))
source("/Users/khuang/bin/LIB_exp.R")
source("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/mutation_impact/mutation_impact.R")

##### BRCA #####
### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted.txt",sep=""))
threshold = 0.02*dim(BRCA_mut)[2]
BRCA_mut_SMG = row.names(BRCA_mut[rowSums(!(BRCA_mut == "wt" | BRCA_mut == "silent" | BRCA_mut == "intronic")) > threshold,])
BRCA_mut_SMG_inKS = BRCA_mut_SMG[(BRCA_mut_SMG %in% c(k_s_table$GENE,k_s_table$SUB_GENE)) & (BRCA_mut_SMG %in% cancer_genes)]
BRCA_mut_SMG_inKS = BRCA_mut_SMG_inKS[BRCA_mut_SMG_inKS %in% c("BRCA1","BRCA2","CDH1","ERBB2","GATA3",
                                        "MAP3K1","PIK3CA","PIK3CB","TP53")]
k_s_table_neighbor = k_s_table[(k_s_table$GENE %in% BRCA_mut_SMG_inKS) | (k_s_table$SUB_GENE %in% BRCA_mut_SMG_inKS),]
k_s_table_neighbor_gene = unique(c(k_s_table_neighbor$GENE,k_s_table_neighbor$SUB_GENE))

BRCA_mut_g = BRCA_mut[row.names(BRCA_mut) %in% BRCA_mut_SMG_inKS,]

### Proteome ###
BRCA_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
BRCA_Pro_c = BRCA_Pro[row.names(BRCA_Pro) %in% k_s_table_neighbor_gene,]
BRCA_Pro_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_Pro_c,name="BRCA_Proteome")

### Phosphoproteome ###
BRCA_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
BRCA_Pho_c = BRCA_Pho[row.names(BRCA_Pho) %in% k_s_table_neighbor_gene,]
BRCA_Pho_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_Pho_c,name="BRCA_Phosphoproteome")

### all levels ###
# do this later, may be more benefitial to do the all outlier table and summarize that instead

##### OV #####
# ### Proteome ###
# OV_JHU_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_PRO_formatted_normalized.txt",sep=""))
# OV_JHU_Pro_c = OV_JHU_Pro[row.names(OV_JHU_Pro) %in% k_s_table_neighbor_gene,]
# OV_JHU_Pro_diff_exp = find_diff_exp(OV_mut_g,OV_JHU_Pro_c,name="OV_JHU_Proteome")
# 
OV_PNNL_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))

### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_SOMATIC_formatted.txt",sep=""))
OV_PNNL_mut = OV_mut[,colnames(OV_mut) %in% colnames(OV_PNNL_Pro)] # only 68 samples but oh well
threshold = 0.02*dim(OV_PNNL_mut)[2]
OV_PNNL_mut_SMG = row.names(OV_PNNL_mut[rowSums(!(OV_PNNL_mut == "wt" | OV_PNNL_mut == "silent" | OV_PNNL_mut == "intronic")) > threshold,])
OV_PNNL_mut_SMG_inKS = OV_PNNL_mut_SMG[OV_PNNL_mut_SMG %in% c(k_s_table$GENE,k_s_table$SUB_GENE) & (OV_PNNL_mut_SMG %in% cancer_genes)]
OV_PNNL_mut_SMG_inKS = OV_PNNL_mut_SMG_inKS[OV_PNNL_mut_SMG_inKS %in% c("BRCA1","BRCA2","KRAS","MTOR","NF1","RB1","TP53")]
k_s_table_neighbor = k_s_table[(k_s_table$GENE %in% OV_PNNL_mut_SMG_inKS) | (k_s_table$SUB_GENE %in% OV_PNNL_mut_SMG_inKS),]
k_s_table_neighbor_gene = unique(c(k_s_table_neighbor$GENE,k_s_table_neighbor$SUB_GENE))

OV_PNNL_mut_g = OV_PNNL_mut[row.names(OV_PNNL_mut) %in% OV_PNNL_mut_SMG_inKS,]

OV_PNNL_Pro_c = OV_PNNL_Pro[row.names(OV_PNNL_Pro) %in% k_s_table_neighbor_gene,]
OV_PNNL_Pro_diff_exp = find_diff_exp(OV_PNNL_mut_g,OV_PNNL_Pro_c,name="OV_PNNL_Proteome")

### Phosphoproteome ###
OV_PNNL_Pho = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))
OV_PNNL_Pho_c = OV_PNNL_Pho[row.names(OV_PNNL_Pho) %in% k_s_table_neighbor_gene,]
OV_PNNL_Pho_diff_exp = find_diff_exp(OV_PNNL_mut_g,OV_PNNL_Pho_c,name="OV_PNNL_Phosphoproteome")

# ### Glycoproteome ###
# OV_JHU_Gly = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/OV/OV_JHU_GLY_formatted_normalized.txt",sep=""))
# OV_JHU_Gly_c = OV_JHU_Gly[row.names(OV_JHU_Gly) %in% druggable,]
# OV_JHU_Gly_diff_exp = find_diff_exp(OV_mut_g,OV_JHU_Gly_c,name="OV JHU Glycoproteome")

### merging the two proteome? ###
### all levels ###

# ##### CRC #####
# ### Mutation matrix ###
# CRC_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_SOMATIC_formatted.txt",sep=""))
# crcGenes = c("TP53","KRAS","APC","PIK3CA","SMAD4")
# k_s_table_neighbor = k_s_table[(k_s_table$GENE %in% crcGenes) | (k_s_table$SUB_GENE %in% crcGenes),]
# k_s_table_neighbor_gene = unique(c(k_s_table_neighbor$GENE,k_s_table_neighbor$SUB_GENE))
# CRC_mut_g = CRC_mut[row.names(CRC_mut) %in% crcGenes,]
# ### Proteome ###
# CRC_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/CRC/CRC_PRO_formatted_normalized.txt",sep=""))
# CRC_Pro_c = CRC_Pro[row.names(CRC_Pro) %in% k_s_table_neighbor_gene,]
# CRC_Pro_diff_exp = find_diff_exp(CRC_mut_g,CRC_Pro_c,name="CRC Proteome") # will likely need to change to non-parametric test
# 

##### merge everything and plot #####

### fold change
BRCA_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_Pro_diff_exp$fold_change)))
colnames(BRCA_Pro_diff_exp_fc_m)[3] = "BRCA_PRO"
BRCA_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$fold_change)))
colnames(BRCA_Pho_diff_exp_fc_m)[3] = "BRCA_PHO"
# OV_JHU_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_JHU_Pro_diff_exp$fold_change)))
# colnames(OV_JHU_Pro_diff_exp_fc_m)[3] = "OV_JHU_PRO"
OV_PNNL_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pro_diff_exp$fold_change)))
colnames(OV_PNNL_Pro_diff_exp_fc_m)[3] = "OV_PRO"
# OV_JHU_Gly_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_JHU_Gly_diff_exp$fold_change)))
# colnames(OV_JHU_Gly_diff_exp_fc_m)[3] = "OV_JHU_GLY"
OV_PNNL_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pho_diff_exp$fold_change)))
colnames(OV_PNNL_Pho_diff_exp_fc_m)[3] = "OV_PHO"
# CRC_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,CRC_Pro_diff_exp$fold_change)))
# colnames(CRC_Pro_diff_exp_fc_m)[3] = "CRC_PRO"

fc = merge(BRCA_Pro_diff_exp_fc_m, BRCA_Pho_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
# fc = merge(fc, OV_JHU_Pro_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
fc = merge(fc, OV_PNNL_Pro_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
#fc = merge(fc, OV_JHU_Gly_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
fc = merge(fc, OV_PNNL_Pho_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
# fc = merge(fc, CRC_Pro_diff_exp_fc_m, by = c("Var1","Var2"), all=T)

fc_m = melt(fc, id.var=c("Var1","Var2"))
colnames(fc_m)[4] = "FC"
## fdr
BRCA_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_Pro_diff_exp$t_fdr)))
colnames(BRCA_Pro_diff_exp_fdr_m)[3] = "BRCA_PRO"
BRCA_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$t_fdr)))
colnames(BRCA_Pho_diff_exp_fdr_m)[3] = "BRCA_PHO"
# OV_JHU_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_JHU_Pro_diff_exp$t_fdr)))
# colnames(OV_JHU_Pro_diff_exp_fdr_m)[3] = "OV_JHU_PRO"
OV_PNNL_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pro_diff_exp$t_fdr)))
colnames(OV_PNNL_Pro_diff_exp_fdr_m)[3] = "OV_PRO"
# OV_JHU_Gly_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_JHU_Gly_diff_exp$t_fdr)))
# colnames(OV_JHU_Gly_diff_exp_fdr_m)[3] = "OV_JHU_GLY"
OV_PNNL_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pho_diff_exp$t_fdr)))
colnames(OV_PNNL_Pho_diff_exp_fdr_m)[3] = "OV_PHO"
# CRC_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,CRC_Pro_diff_exp$t_fdr)))
# colnames(CRC_Pro_diff_exp_fdr_m)[3] = "CRC_PRO"

fdr = merge(BRCA_Pro_diff_exp_fdr_m, BRCA_Pho_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
# fdr = merge(fdr, OV_JHU_Pro_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
fdr = merge(fdr, OV_PNNL_Pro_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
#fdr = merge(fdr, OV_JHU_Gly_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
fdr = merge(fdr, OV_PNNL_Pho_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
# fdr = merge(fdr, CRC_Pro_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)

fdr_m = melt(fdr, id.var=c("Var1","Var2"))
colnames(fdr_m)[4] = "FDR"

fc_fdr = merge(fc_m, fdr_m, by=c("Var1","Var2","variable"))

## fdr
BRCA_Pro_diff_exp_p_m = melt(as.matrix(do.call(cbind,BRCA_Pro_diff_exp$t_p)))
colnames(BRCA_Pro_diff_exp_p_m)[3] = "BRCA_PRO"
BRCA_Pho_diff_exp_p_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$t_p)))
colnames(BRCA_Pho_diff_exp_p_m)[3] = "BRCA_PHO"
OV_PNNL_Pro_diff_exp_p_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pro_diff_exp$t_p)))
colnames(OV_PNNL_Pro_diff_exp_p_m)[3] = "OV_PRO"
OV_PNNL_Pho_diff_exp_p_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pho_diff_exp$t_p)))
colnames(OV_PNNL_Pho_diff_exp_p_m)[3] = "OV_PHO"

p = merge(BRCA_Pro_diff_exp_p_m, BRCA_Pho_diff_exp_p_m, by = c("Var1","Var2"), all=T)
p = merge(p, OV_PNNL_Pro_diff_exp_p_m, by = c("Var1","Var2"), all=T)
p = merge(p, OV_PNNL_Pho_diff_exp_p_m, by = c("Var1","Var2"), all=T)

p_m = melt(p, id.var=c("Var1","Var2"))
colnames(p_m)[4] = "P"

fc_fdr = merge(fc_fdr, p_m, by=c("Var1","Var2","variable"))

fc_fdr$both_FDR = p.adjust(fc_fdr$P, method="BH")

# look for the kinase-substrate relations
fc_fdr$KSrelation = "none"
fc_fdr$KSrelation[paste(fc_fdr$Var2,fc_fdr$Var1) %in% paste(k_s_table$GENE,k_s_table$SUB_GENE)] = "Downstream"
fc_fdr$KSrelation[paste(fc_fdr$Var2,fc_fdr$Var1) %in% paste(k_s_table$SUB_GENE,k_s_table$GENE)] = "Upstream"
fc_fdr$KSrelation[as.character(fc_fdr$Var2) == as.character(fc_fdr$Var1)] = "CIS"

fc_fdr$KSpair = (paste(fc_fdr$Var2,fc_fdr$Var1) %in% paste(k_s_table$GENE,k_s_table$SUB_GENE)) | 
                (paste(fc_fdr$Var2,fc_fdr$Var1) %in% paste(k_s_table$SUB_GENE,k_s_table$GENE)) 
                
fc_fdr$cis = (as.character(fc_fdr$Var2) == as.character(fc_fdr$Var1))
fc_fdr[fc_fdr$Var2=="PIK3CB" & fc_fdr$Var1=="AKT1",]$KSpair = TRUE # add this likely pair
fc_fdr[fc_fdr$Var2=="PIK3CB" & fc_fdr$Var1=="AKT1",]$KSrelation = "Downstream"

fn = paste(pd, "_mutation_differential_exp_BRCA_OV.txt")
write.table(fc_fdr, file = fn, quote=F, sep = '\t', col.names=NA)

## plot
fdr.colors=c("NA", "#000000")
min_d = min(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
max_d = max(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
bound = max(c(max_d, -min_d))
fc_fdr$sig = as.numeric(as.character(fc_fdr$both_FDR)) < 0.05
fc_fdr$FC_2 = as.numeric(fc_fdr$FC)
fc_fdr[!is.na(fc_fdr$FC_2) & fc_fdr$FC_2>=2,]$FC_2=2
fc_fdr[!is.na(fc_fdr$FC_2) & fc_fdr$FC_2<=-2,]$FC_2=-2
fc_fdr.v = fc_fdr[rowSums(is.na(fc_fdr))<4,]
## dropping lvls, doesn't work yet
# fc_fdr.v = droplevels(fc_fdr.v)
# fc_fdr.v$Var1 = factor(fc_fdr.v$Var1)
# fc_fdr.v$Var2 = factor(fc_fdr.v$Var2)
# fc_fdr.v[fc_fdr.v$variable == "BRCA_PRO",]$Var2 = droplevels(fc_fdr.v[fc_fdr.v$variable == "BRCA_PRO",]$Var2) 
# recalculate FDR
fc_fdr_s = fc_fdr.v[fc_fdr.v$KSrelation != "none",]
fc_fdr_s$KSonly_fdr = p.adjust(fc_fdr_s$P, method="BH")

fc_fdr_s$sig = as.numeric(as.character(fc_fdr_s$KSonly_fdr)) < 0.05

fn = paste(pd, "_mutation_differential_exp_BRCA_OV_KS_only.txt")
write.table(fc_fdr_s, file = fn, quote=F, sep = '\t', col.names=NA)

fc_fdr_s$Var2 = factor(fc_fdr_s$Var2, levels=as.vector(unique(fc_fdr_s$Var2))[order(as.vector(unique(fc_fdr_s$Var2)))])

fn = paste(pd, 'merged_diff_exp_cgenes_k-s_cis_p.pdf',sep ="_")
p = ggplot(data=fc_fdr_s)
p = p + facet_grid(variable~., drop=T, scales = "free", space = "free")
#p = p + facet_wrap(~variable, nrow=1, scales = "free_x", space = "free_x", drop=TRUE)
p = p + geom_point(aes(y=Var2, x=as.numeric(FC_2), fill=KSrelation, size=-log10(as.numeric(KSonly_fdr)), color=ifelse(sig, "black",NA)),pch=21) 
p = p + geom_text(aes(y=Var2, x=as.numeric(FC_2),label=ifelse(sig | cis,as.character(Var1),NA)),angle=45,size=2)
#p = p + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
p = p + geom_vline(xintercept = 0, alpha=0.5)
p = p + scale_colour_manual(values=c("black",NA))
p = p + labs(y="Mutated Gene", x="Expression fold change") + theme_bw() 
#p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())
p = p + coord_fixed()
p
ggsave(file=fn, width=8, height=6, useDingbats=FALSE)


### plot selected markers ###
for (i in 3:6){
  fdr[,i] = as.numeric(as.character(fdr[,i]))
}
markers = unique(fdr[rowSums(fdr[,c(3:6)]<0.005, na.rm=T) >=1,]$Var1)
fc_fdr_s = fc_fdr.v[as.character(fc_fdr.v$Var1) %in% as.character(markers),]

fn = paste(pd, 'merged_diff_exp_cgenes_kinase_substrate_fdr0.005in1_p.pdf',sep ="_")
p = ggplot(data=fc_fdr_s)
p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
#p = p + facet_wrap(~variable, nrow=1, scales = "free_x", space = "free_x", drop=TRUE)
p = p + geom_point(aes(x=Var2, y=Var1, fill=as.numeric(FC_2), size=-log10(as.numeric(FDR)), color=ifelse(sig, "black",NA)),pch=21) 
p = p + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
p = p + scale_colour_manual(values=c("black",NA))
p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
p = p + coord_fixed()
p
ggsave(file=fn, width=12, height=8, useDingbats=FALSE)

### for grant ### 
if (FALSE){
  markers = c("AR","BRAF","CDH1","CHEK2","EGFR","ERBB2","ESR1","GATA3","PIK3CA","PDGFRA","TP53")
  fc_fdr_s = fc_fdr.v[fc_fdr.v$Var1 %in% as.character(markers),]
  fc_fdr_s = fc_fdr_s[fc_fdr_s$variable %in% c("BRCA_PRO","BRCA_PHO","OV_PNNL_PRO","OV_PNNL_PHO","CRC_PRO"),]
  fc_fdr_s$variable = fc_fdr_s$variable[drop=T]
  levels(fc_fdr_s$variable) = c("BRCA_PRO","BRCA_PHO","OV_PRO","OV_PHO","CRC_PRO")
  
  fn = paste(pd, 'merged_diff_exp_cgenes_druggable_U24_label.pdf',sep ="_")
  p = ggplot(data=fc_fdr_s)
  #p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
  p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
  p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
  p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
  p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
    theme(axis.title = element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=6), 
          axis.text.y = element_text(colour="black", size=6), axis.ticks = element_blank(),
          strip.text.x = element_text(size = 6))#element_text(colour="black", size=16))
  p = p + coord_equal()
  #p = p + theme(legend.position="none")
  p
  ggsave(file=fn, width=10, height=5, useDingbats=FALSE)
}

### for PTRC grant ###
if (FALSE){
  # exclude TSG
  oncogenes_f = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/Volgestin_Oncogene_list.txt', header=FALSE, stringsAsFactors = F)
  oncogenes = as.vector(t(oncogenes_f))
  fc_fdr_s_onco = fc_fdr_s[fc_fdr_s$Var1 %in% oncogenes & fc_fdr_s$variable != "OV_JHU_PRO",]
  fc_fdr_s_onco$variable = gsub("_.*_","_", fc_fdr_s_onco$variable)
  fc_fdr_s_onco$variable = factor(fc_fdr_s_onco$variable, levels = c("BRCA_PRO","BRCA_PHO","OV_PRO","OV_PHO","CRC_PRO"))
  
  fn = paste(pd, 'merged_diff_exp_cgenes_druggable_fdr0.1in1_oncogenes.pdf',sep ="_")
  p = ggplot(data=fc_fdr_s_onco)
  #p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
  p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
  p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
  p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
  p = p + labs(x="Mutated Gene", y="Expression") + theme_nogrid() + 
    theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
  p = p + coord_fixed()
  p
  ggsave(file=fn, width=8.5, height=4, useDingbats=FALSE)
  
  
  fc_fdr_s = fc_fdr.v[fc_fdr.v$Var1 %in% as.character(oncogenes),]
  fc_fdr_s = fc_fdr_s[fc_fdr_s$variable %in% c("BRCA_PRO","BRCA_PHO","OV_PNNL_PRO","OV_PNNL_PHO","CRC_PRO"),]
  fc_fdr_s$variable = fc_fdr_s$variable[drop=T]
  levels(fc_fdr_s$variable) = c("BRCA_PRO","BRCA_PHO","OV_PRO","OV_PHO","CRC_PRO")
  
  fn = paste(pd, 'merged_diff_exp_cgenes_druggable_oncogenes.pdf',sep ="_")
  p = ggplot(data=fc_fdr_s)
  #p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
  p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
  p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
  p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
  p = p + labs(x="Mutated Gene", y="Expression") + theme_nogrid() + 
    theme(axis.title = element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=6), 
          axis.text.y = element_text(colour="black", size=6), axis.ticks = element_blank(),
          strip.text.x = element_text(size = 6))#element_text(colour="black", size=16))
  p = p + coord_equal()
  #p = p + theme(legend.position="none")
  p
  ggsave(file=fn, width=7, height=5, useDingbats=FALSE)
}
##### plot collective violin plot #####
### make violin plots for selected mutated gene-marker combo###

plot_diff_exp_violin_c = function (Mgene, gene){
  BRCA_mut_gene = data.frame(t(BRCA_mut[Mgene,]))
  BRCA_mut_gene$carrier = as.character(BRCA_mut_gene[,1]) != "wt" & as.character(BRCA_mut_gene[,1]) != "silent"
  BRCA_Pro_gene = t(BRCA_Pro[gene,])
  BRCA_Pho_gene = t(BRCA_Pho[gene,])
  # merge by sample
  BRCA_Pro_merge = merge(BRCA_mut_gene, BRCA_Pro_gene, by = "row.names")
  BRCA_Pro_merge$data = "BRCA_PRO"
  BRCA_Pho_merge = merge(BRCA_mut_gene, BRCA_Pho_gene, by = "row.names")
  BRCA_Pho_merge$data = "BRCA_PHO"
  
  CRC_mut_gene = data.frame(t(CRC_mut[Mgene,]))
  CRC_mut_gene$carrier = as.character(CRC_mut_gene[,1]) != "wt" & as.character(CRC_mut_gene[,1]) != "silent"
  CRC_Pro_gene = t(CRC_Pro[gene,])
  CRC_Pro_merge = merge(CRC_mut_gene, CRC_Pro_gene, by = "row.names")
  CRC_Pro_merge$data = "CRC_PRO"
  
  OV_mut_gene = data.frame(t(OV_mut[Mgene,]))
  OV_mut_gene$carrier = as.character(OV_mut_gene[,1]) != "wt" & as.character(OV_mut_gene[,1]) != "silent"
  OV_JHU_Pro_gene = t(OV_JHU_Pro[gene,])
  OV_JHU_Pro_merge = merge(OV_mut_gene, OV_JHU_Pro_gene, by = "row.names")
  OV_JHU_Pro_merge$data = "OV_JHU_PRO"
  OV_PNNL_Pro_gene = t(OV_PNNL_Pro[gene,])
  OV_PNNL_Pro_merge = merge(OV_mut_gene, OV_PNNL_Pro_gene, by = "row.names")
  OV_PNNL_Pro_merge$data = "OV_PNNL_PRO"
  OV_PNNL_Pho_gene = t(OV_PNNL_Pho[gene,])
  OV_PNNL_Pho_merge = merge(OV_mut_gene, OV_PNNL_Pho_gene, by = "row.names")
  OV_PNNL_Pho_merge$data = "OV_PNNL_PHO"
  
  colnames(BRCA_Pro_merge)[4]=gene
  colnames(BRCA_Pho_merge)[4]=gene
  colnames(CRC_Pro_merge)[4]=gene
  colnames(OV_JHU_Pro_merge)[4]=gene
  colnames(OV_PNNL_Pro_merge)[4]=gene
  colnames(OV_PNNL_Pho_merge)[4]=gene
  
  gene_all_lvl = rbind(BRCA_Pro_merge,BRCA_Pho_merge,CRC_Pro_merge,OV_JHU_Pro_merge,OV_PNNL_Pro_merge,OV_PNNL_Pho_merge)
  colnames(gene_all_lvl) = c("Sample","Mutation_Type","Mutation_Status","Expression","Dataset")
  
  # plot violin plots faceted by marker genes
  fn = paste(pd, Mgene, gene, "mutational_impact_violin.pdf", sep="_")
  p = ggplot(data=gene_all_lvl)
  p = p + facet_grid(.~Dataset)
  p = p + geom_violin(aes(x=Mutation_Status, y=Expression, fill=Mutation_Status),alpha=0.5) + guides(fill=FALSE) 
  p = p + geom_jitter(aes(x=Mutation_Status, y=Expression, color=Mutation_Type)) #+ geom_point(aes(x=Status, y=value)) 
  p = p + labs(x = paste(Mgene,"Mutation Status"), y = paste(gene, "Expression")) + theme_bw()
  p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
  p
  ggsave(file=fn, width=14.5, limitsize=FALSE, useDingbats=FALSE)
}

# plot any pair of interest
if (FALSE){
  plot_diff_exp_violin_c("TP53","TP53")
  plot_diff_exp_violin_c("TP53","ESR1")
  plot_diff_exp_violin_c("TP53","IGF1R")
  plot_diff_exp_violin_c("TP53","GATA3")
  plot_diff_exp_violin_c("TP53","EGFR")
  plot_diff_exp_violin_c("TP53","CDH1")
  plot_diff_exp_violin_c("GATA3","EGFR")
  plot_diff_exp_violin_c("TP53","CHEK2")
  plot_diff_exp_violin_c("NF1","TP53")
  plot_diff_exp_violin_c("CDH1","CDH1")
  plot_diff_exp_violin_c("CDH1","PDGFRB")
  plot_diff_exp_violin_c("KRAS","MAP2K1")
}

gene="MDM2"
