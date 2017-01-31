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

# input -------------------------------------------------------------------
cohort <- "OV"
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


# calculate how many tests are needed -------------------------------------
ntest <- 0
for (gene in unique_substrate) {
  temp <- length(which(pho_rsd_split$SUBSTRATE==gene))
  ntest <- ntest + temp*(temp+1)/2
}
ntest

# initiate ----------------------------------------------------------------
## substrate, rsd1,rsd2, pvalue, coef, genomic.pos1, genomic.pos2, 
template <- data.frame(t(vector(mode = "numeric", length = 7 )))
colnames(template) <- c("SUBSTRATE","RSD1","RSD2","fdr","pvalue","coef","size")

# loop around each protein ------------------------------------------------
phosite_corr <- c()
start.time <- Sys.time()
start.time
for (g in 1:length(unique_substrate)) {
#for (g in 1:10) {
  gene <- unique_substrate[g]
  rows <- which(pho_rsd_split$SUBSTRATE==gene)
  n <- length(rows)
  if (n > 1) {
    template1 <- template
    template1$SUBSTRATE <- gene
    
    for (i in 1:(n-1)) {
      ri <- rows[i]
      rsd1 <- pho_rsd_split$SUB_MOD_RSD[ri];#pos1 <- pho_data$genomic_pos[i]
      pho_temp1 <- pho_main[ri,]
      pho_temp1_norm <- range01(unlist(pho_temp1),na.rm = T)
      
      template2 <- template1
      template2$RSD1 <- rsd1;#template2$pos1 <- pos1
    
      for (j in (i+1):n ) {
        rj <- rows[j]
        rsd2 <- pho_rsd_split$SUB_MOD_RSD[rj];#pos2 <- pho_data$genomic_pos[j]
        
        pho_temp2 <- pho_main[rj,]
        pho_temp2_norm <- range01(unlist(pho_temp2),na.rm = T)
        
        #prepare regression data for model1
        data1 <- data.frame(pho_temp1_norm,pho_temp2_norm)
        
        size <- nrow(data1[complete.cases(data1),])
        if( size > least_samples ){#more than 2 complete dataset
          temp <- template2
          temp$RSD2 <- rsd2; #temp$pos2 <- pos2;
          temp$size <- size
          fit1 <- glm(pho_temp1_norm ~ pho_temp2_norm,data = data1, family=gaussian())
          
          temp$pvalue <- c(coef(summary(fit1))[2,4])
          temp$coef <- fit1$coefficients[2]
          phosite_corr <- rbind(phosite_corr,temp)
        }
      }
    }
  }
  print(g)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# adjust for FDR ----------------------------------------------------------
pos1 <- vector(mode = "character", length = nrow(phosite_corr));
pos2 <- pos1

col_substrate <- as.vector(pho_rsd_split$SUBSTRATE)
col_rsd <- as.vector(pho_rsd_split$SUB_MOD_RSD)
col_genomic <- as.vector(pho_data$genomic_pos)
for (i in 1:nrow(pho_data)) {
  rows <- which((phosite_corr$SUBSTRATE==col_substrate[i]) & (phosite_corr$RSD1==col_rsd[i]))
  if ( length(rows) >0 ) {
    pos1[rows] <- col_genomic[i]
  }
}
for (i in 1:nrow(pho_data)) {
  rows <- which((phosite_corr$SUBSTRATE==col_substrate[i]) & (phosite_corr$RSD2==col_rsd[i]))
  if ( length(rows) >0 ) {
    pos2[rows] <- col_genomic[i]
  }
}

phosite_corr <- cbind(phosite_corr,pos1,pos2)
phosite_corr <- phosite_corr[phosite_corr$pos1 != as.character(phosite_corr$pos2),]

phosite_corr$fdr <-p.adjust(phosite_corr$pvalue,method = "fdr")
phosite_corr$pair <- paste(phosite_corr$SUBSTRATE,phosite_corr$RSD1,phosite_corr$RSD2,sep = ":")
phosite_corr$fdr_log10 <- -log10(phosite_corr$fdr)

## write out tables
tn = paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/",cohort,"_phosphosite_correlation.txt", sep="")
write.table(phosite_corr, file=tn, quote=F, sep = '\t', row.names = FALSE)

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


# input pairwise hotspot3d result -----------------------------------------
pairwise <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/cptac.brca.phosphosites.pairwise.singleprotein.collapsed_for_R.txt",sep = ""))
coef_corr <- vector(mode = "numeric", length = nrow(pairwise)) + NA
fdr_corr <- coef_corr
gene_p <- as.vector(pairwise$Gene1)
rsd1_list <- paste("p.",phosite_corr$RSD1,"C",sep = "")
rsd2_list <- paste("p.",phosite_corr$RSD2,"C",sep = "")

for (i in 1:nrow(pairwise)) {
  rows <- which( (phosite_corr$SUBSTRATE==gene_p[i]) & (rsd1_list==pairwise$RSD1[i]) & (rsd2_list==pairwise$RSD2[i]))
  if (length(rows) > 0) {
    coef_corr[i] <- phosite_corr$coef[rows]
    fdr_corr[i] <- phosite_corr$fdr[rows]
  }
}
pairwise <- cbind(pairwise,coef_corr,fdr_corr)


# plot correlation between coef_corr and distances ------------------------

p <- ggplot(data = pairwise, aes(x = dis_3d , y = coef_corr)) +
  geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ylim(c(-1,2))
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

fn = paste(baseD,'pan3can_shared_data/analysis_results/hotspot3d/',cohort,'_distance_3d_and_coef_corr_correlation.png',sep ="")
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
