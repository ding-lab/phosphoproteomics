# Yige Wu @ WashU 2017 Jan
# find correlation of phosphorylation level within each protein

# directory and library ---------------------------------------------------
# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"

library(stringr)
library(ggplot2)


# input regardless to cancer type, choose between kinase or phosphotase-------------------------------------------------------------------
### read in the kinase/substrate table/ phosphorylation data ###
K_S_f = paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep="")
k_s_table = read.delim(K_S_f)

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}


unique_kinase <- unique(k_s_table$GENE)
unique_substrate <- as.vector(unique(k_s_table$SUB_GENE))

least_samples <- 5
out_thres <- 1.5
sig <- 0.05

# input phosphorylation data-------------------------------------------------------------------
cohort <- "BRCA"
if (cohort == "BRCA") {
  BRCA_pho_f = paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA_wGpos.txt",sep="")
  pho_data = read.delim(BRCA_pho_f)
  # ordering the columns by sample name
  pho_data <- pho_data[,order(names(pho_data))]
  pho_main <- pho_data[,1:77]
}

if (cohort == "OV") {
  OV_pho_f = paste(baseD,"pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA_wGpos.txt",sep="")
  pho_data = read.delim(OV_pho_f)
  # ordering the columns by sample name
  pho_data <- pho_data[,order(names(pho_data))]
  pho_main <- pho_data[,3:71]
}

#split the SUBSTRATE and SUB_MOD_RSD in the first column
pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))

#covert the SUB_MOD_RSD from lowercase to uppercase
pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")

# input pairwise hotspot3d result -----------------------------------------
pairwise <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/cptac.brca.phosphosites.pairwise.singleprotein.collapsed_for_R.txt",sep = ""))
gene_p_all <- as.vector(unique(pairwise$Gene1))
gene_p <- data.frame(gene_p_all)
colnames(gene_p) <- c("gene")
gene_p$nrsd <- vector(mode = "numeric", length = nrow(gene_p))
gene_p$ntest <- vector(mode = "numeric", length = nrow(gene_p))
rownames(gene_p) <- gene_p$gene

# calculate how many tests are needed if doing all the combination of genes in pairwise file -------------------------------------
for (i in 1:nrow(gene_p)) {
  gene <- gene_p_all[i]
  gene_p$nrsd[i] <- length(which(pho_rsd_split$SUBSTRATE==gene))
  gene_p$ntest[i] <- (gene_p$nrsd[i])^2 - gene_p$nrsd[i]
}
sum(gene_p$ntest)

# initiate ----------------------------------------------------------------
## substrate, rsd1,rsd2, pvalue, coef, genomic.pos1, genomic.pos2, 
template <- data.frame(t(vector(mode = "numeric", length = 9 )))
colnames(template) <- c("GENE1","RSD1","RSD2","fdr","pvalue","coef","size", "pos1", "pos2")
gene1_p <- gene_p_all[gene_p$nrsd >= 2]

# cor.test, only those exact gene pairs in the pairwise file ------------------------------------------------
phosite_corr <- c()
start.time <- Sys.time()
start.time
#for (i in 1:nrow(gene_pair)) {
for (i in 1:10) {
  gene1 <- gene1_p[i]
  n1 <- which(pho_rsd_split$SUBSTRATE==gene1)

  template$GENE1 <- gene1

  for( r1 in 1:(length(n1)-1) ){
    j1 <- n1[r1]
    template$RSD1 <- pho_rsd_split$SUB_MOD_RSD[j1]
    template$pos1 <- pho_data$genomic_pos[j1]
    
    pho_temp1 <- pho_main[j1,]
    pho_temp1_norm <- range01(unlist(pho_temp1),na.rm = T)
    
    for( r2 in r1:length(n1) ){
      j2 <- n1[r2]
      pho_temp2 <- pho_main[j2,]
      pho_temp2_norm <- range01(unlist(pho_temp2),na.rm = T)
      
      corr_stat = try(cor.test(pho_temp1_norm,pho_temp2_norm, method = "pearson"), silent=T)
      
      if (!is(corr_stat,"try-error")){ # if the correlation test is carried out successfully
        template$RSD2 <- pho_rsd_split$SUB_MOD_RSD[j2]
        template$pos2 <- pho_data$genomic_pos[j2]
        
        data1 <- data.frame(pho_temp1_norm,pho_temp2_norm)
        template$size <- nrow(data1[complete.cases(data1),])
        
        template$coef <- corr_stat$estimate
        template$pvalue <-corr_stat$p.value
        
        phosite_corr <- rbind(phosite_corr,template)
      }
    }
  }
  print(i)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# adjust for FDR ----------------------------------------------------------
phosite_corr <- phosite_corr[phosite_corr$pos1 != as.character(phosite_corr$pos2),]
phosite_corr <- phosite_corr[phosite_corr$size >= 5,]
phosite_corr <- unique(phosite_corr)

phosite_corr$fdr <-p.adjust(phosite_corr$pvalue,method = "fdr")
phosite_corr$pair <- paste(phosite_corr$GENE1,phosite_corr$RSD1,phosite_corr$GENE2,phosite_corr$RSD2,sep = ":")
phosite_corr$fdr_log10 <- -log10(phosite_corr$fdr)

## write out tables
tn = paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/",cohort,"_phosphosite_across_protein_correlation.txt", sep="")
write.table(phosite_corr, file=tn, quote=F, sep = '\t', row.names = FALSE)


# fill the cor_stat into pairwise table -----------------------------------
coef_corr <- vector(mode = "numeric", length = nrow(pairwise)) + NA
fdr_corr <- coef_corr
rsd1_list <-  str_split_fixed(pairwise$RSD1,"[p.C]",4)[,3]
rsd2_list <-  str_split_fixed(pairwise$RSD2,"[p.C]",4)[,3]
gene1_list <- as.vector(pairwise$Gene1)
gene2_list <- as.vector(pairwise$Gene2)


for (i in 1:nrow(pairwise)) {
  rows <- which( (phosite_corr$GENE1==gene1_list[i]) & (phosite_corr$GENE2==gene2_list[i]) & (phosite_corr$RSD1==rsd1_list[i]) & (phosite_corr$RSD2==rsd2_list[i]) )
  if (length(rows) > 0) {
    print(rows)
    coef_corr[i] <- phosite_corr$coef[rows]
    fdr_corr[i] <- phosite_corr$fdr[rows]
  }
}
pairwise$coef_corr <- coef_corr
pairwise$fdr_corr <- fdr_corr

## write out tables
tn = paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/",cohort,"_phosphosite_within_protein_distance_and_correlation.txt", sep="")
write.table(pairwise, file=tn, quote=F, sep = '\t', row.names = FALSE)


# volcano plotting module -------------------------------------------------
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

phosite_corr$coef_filtered = remove_outliers(phosite_corr$coef)
phosite_corr_outlier_removed_m = phosite_corr[!is.na(phosite_corr$coef_filtered) & phosite_corr$fdr_log10 < 100 ,]
# phosite_corr_outlier_removed_m = phosite_corr[!is.na(phosite_corr$coef_filtered),]
plot_fdr_scale <- 50


p = ggplot(phosite_corr_outlier_removed_m,aes(x=coef, y=-log10(fdr)))
p = p + geom_point(alpha=0.005)
p = p + geom_text(aes(label= ifelse(-log10(fdr)>plot_fdr_scale, pair, NA)),size=2,alpha=0.2)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient", y="-log10(FDR)")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/hotspot3d/',cohort,'_substrate_phosphosite_correlation_volcano.pdf',sep ="")
ggsave(file=fn, height=6, width=8, useDingbats=FALSE)
dev.off()


# plot correlation between coef_corr and distances ------------------------

p <- ggplot(data = pairwise, aes(x = dis_3d , y = coef_corr)) +
  geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ylim(c(-1,1))
p
# lm_eqn = function(df){
#   m = glm(coef_corr ~ dis_3d, data = df);
#   eq <- substitute(italic(coef_corr) == a + b %.% italic(dis_3d),
#                    list(a = format(coef(m)[1], digits = 2), 
#                         b = format(coef(m)[2], digits = 2)));
#   as.character(as.expression(eq));
# }

# eq <- lm_eqn(pairwise)
# p2 = p + geom_text(data=eq,aes(x = 5, y = 1.9,label=V1, size=5), parse = TRUE, inherit.aes=FALSE)
# p2

fn = paste(baseD,'pan3can_shared_data/analysis_results/hotspot3d/across_protein_',cohort,'_distance_3d_and_coef_corr_correlation.png',sep ="")
ggsave(file=fn, height=6, width=8)

fit1 <- glm(coef_corr ~ dis_3d ,data = pairwise[pairwise$coef_corr>-1,], family=gaussian())
Pvalue=c(coef(summary(fit1))[2,4]); Pvalue
Coef=fit1$coefficients[2]; Coef

pairwise$dis_linear <- as.numeric(pairwise$dis_linear)
p <- ggplot(data = pairwise[!is.na(pairwise$coef_corr),], aes(x = dis_linear , y = coef_corr)) +
  geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ylim(c(-1,2))
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/hotspot3d/distance_linear_and_coef_corr_correlation.png',sep ="")
ggsave(file=fn, height=6, width=8)

fit1 <- glm(coef_corr ~ dis_linear ,data = pairwise, family=gaussian())
Pvalue=c(coef(summary(fit1))[2,4]); Pvalue
Coef=fit1$coefficients[2]; Coef


# overlap with regression result and this ---------------------------------
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/kinase_substrate_regression.txt",sep = ""))
table1 <- table_2can[table_2can$model=="pho_sub~pro_kin",]
rsd1_p <- str_split_fixed(pairwise$RSD1,"[p.C]",4)[,3]
rsd2_p <- str_split_fixed(pairwise$RSD2,"[p.C]",4)[,3]
cancer <- cohort
#for (i in 1:nrow(pairwise)) {
for (i in which(pairwise$fdr_corr<=0.05)) {
  k1 <- which(table1$SUBSTRATE==gene_p[i] & table1$SUB_MOD_RSD==rsd1_p[i] & table1$Cancer==cancer)
  k2 <- which(table1$SUBSTRATE==gene_p[i] & table1$SUB_MOD_RSD==rsd2_p[i] & table1$Cancer==cancer)
  kinase <- unique(table1$KINASE[k1],table1$KINASE[k2])
  if (length(kinase) > 0) {
    for (j in kinase) {
      p1 <- table1$FDR_pro_kin[table1$KINASE==j & table1$SUBSTRATE==gene_p[i] & table1$SUB_MOD_RSD==rsd1_p[i] & table1$Cancer==cancer]
      p2 <- table1$FDR_pro_kin[table1$KINASE==j & table1$SUBSTRATE==gene_p[i] & table1$SUB_MOD_RSD==rsd2_p[i] & table1$Cancer==cancer]
      
      if ( p1 <= sig && p2 <= sig ) {
        print(paste(gene_p[i],rsd1_p[i],rsd2_p[i],pairwise$fdr_corr[i],j,p1,p2,sep = ":"))
      }
    }
  }
}



# extract RTK results -----------------------------------------------------
RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
rtk <- as.vector(t(RTK_file))
rows <- c()
for (gene in rtk) {
  temp <- which(pairwise$Gene1==gene)
  rows <- c(rows,temp)
}
p <- ggplot(data = pairwise[rows,], aes(x = dis_linear , y = coef_corr)) +
  geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ylim(c(-1,2))
p
