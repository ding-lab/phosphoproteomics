# Yige Wu @ WashU 2017 Jan
# find correlation of phosphorylation level within each protein

# directory and library ---------------------------------------------------
# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"

library(stringr)
library(ggplot2)
library(readr)

# choose cohort and other parameters --------------------------------------
sig <- 0.05
cohort <- "BRCA"
# cohort <- "OV"

# other input-------------------------------------------------------------------
#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

# input phosphorylation data-------------------------------------------------------------------
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
pairwise <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/PDB_cptac.brca.ov_sites.pairwise.singleprotein.collapsed",sep = ""), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
pairwise <- data.frame(pairwise[,c(1,2,6,7,12,13,14,15)])
colnames(pairwise) <- c("Gene1","RSD1","Gene2","RSD2","dis_linear","dis_3d","PDB","pvalue_3d")

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
  if (i==1) {
    gene_p$nrsd_cum[i] <- 0
  }
  if (i > 1) {
    gene_p$nrsd_cum[i] <- gene_p$nrsd_cum[i-1] + gene_p$nrsd[i-1]
  }
}
nr <- sum(gene_p$ntest)

# initiate ----------------------------------------------------------------
## substrate, rsd1,rsd2, pvalue, coef, genomic.pos1, genomic.pos2, 

gene1_p <- gene_p_all[gene_p$nrsd >= 2]
pos_list <- as.vector(pho_data$genomic_pos)

x <-  vector(mode = "character", length = nr )
y <-  vector(mode = "numeric", length = nr )
GENE1 <- x; RSD1 <- x; RSD2 <- x; POS1 <- x; POS2 <- x;
pvalue <- y; coef <- y; size <- y; fdr <- y;

# cor.test, only those exact gene pairs in the pairwise file ------------------------------------------------
start.time <- Sys.time()
start.time
row <- 0
for (i in 1:length(gene1_p)) {
#for (i in 1:10) {
  gene1 <- gene1_p[i]
  n1 <- which(pho_rsd_split$SUBSTRATE==gene1)

  for( r1 in 1:(length(n1)-1) ){
    j1 <- n1[r1]
    rsd1 <- pho_rsd_split$SUB_MOD_RSD[j1]
    pos1 <- pos_list[j1]
    
    pho_temp1 <- pho_main[j1,]
    pho_temp1_norm <- range01(unlist(pho_temp1),na.rm = T)
    
    for( r2 in (r1+1):length(n1) ){
      row <- row + 1
      j2 <- n1[r2]
      
      pho_temp2 <- pho_main[j2,]
      pho_temp2_norm <- range01(unlist(pho_temp2),na.rm = T)
      
      corr_stat = try(cor.test(pho_temp1_norm,pho_temp2_norm, method = "pearson"), silent=T)
      if (!is(corr_stat,"try-error")){ # if the correlation test is carried out successfully
        GENE1[row] <- gene1
        RSD1[row] <- rsd1
        POS1[row] <- pos1
        RSD2[row] <- pho_rsd_split$SUB_MOD_RSD[j2]
        POS2[row] <- pos_list[j2]
        size[row] <- length(which(!is.na(pho_temp1_norm) & !is.na(pho_temp2_norm)))
        
        # data1 <- data.frame(pho_temp1_norm,pho_temp2_norm)
        # template$size <- nrow(data1[complete.cases(data1),])
        
        coef[row] <- corr_stat$estimate
        pvalue[row] <-corr_stat$p.value
      }
    }
  }
  print(i)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
phosite_corr <- data.frame(GENE1,RSD1,RSD2,fdr,coef,pvalue,size,POS1,POS2)
colnames(phosite_corr) <- c("GENE1","RSD1","RSD2","fdr","coef","pvalue","size", "pos1", "pos2")

# adjust for FDR ----------------------------------------------------------
phosite_corr <- phosite_corr[phosite_corr$pos1 != as.character(phosite_corr$pos2),]
phosite_corr <- phosite_corr[phosite_corr$size >= 5,]
phosite_corr <- unique(phosite_corr)

phosite_corr$fdr <-p.adjust(phosite_corr$pvalue,method = "fdr")
phosite_corr$pair <- paste(phosite_corr$GENE1,phosite_corr$RSD1,phosite_corr$RSD2,sep = ":")
phosite_corr$fdr_log10 <- -log10(phosite_corr$fdr)

## write out tables
tn = paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/",cohort,"_phosphosite_within_protein_correlation.txt", sep="")
write.table(phosite_corr, file=tn, quote=F, sep = '\t', row.names = FALSE)


# fill the cor_stat into pairwise table -----------------------------------
coef_corr <- vector(mode = "numeric", length = nrow(pairwise)) + NA
fdr_corr <- coef_corr;
rsd1 <- vector(mode = "character", length = (nrow(pairwise))); rsd2 <- rsd1
rsd1_list <-  str_split_fixed(pairwise$RSD1,"[p.XSTYCZ]",5)[,4]
rsd2_list <-  str_split_fixed(pairwise$RSD2,"[p.XSTYCZ]",5)[,4]
gene1_list <- as.vector(pairwise$Gene1)
pairwise$dis_lin <- abs(as.numeric(rsd2_list)-as.numeric(rsd1_list))
RSD1_list <- as.vector(phosite_corr$RSD1); RSD2_list <- as.vector(phosite_corr$RSD2)

rsd1_corr <- str_split_fixed(phosite_corr$RSD1,"[STY]",3)[,2]
rsd2_corr <- str_split_fixed(phosite_corr$RSD2,"[STY]",3)[,2]
  
for (i in 1:nrow(pairwise)) {
  rows1 <- which( (phosite_corr$GENE1==gene1_list[i])  & (rsd1_corr==rsd1_list[i]) & (rsd2_corr==rsd2_list[i]) )
  rows2 <- which( (phosite_corr$GENE1==gene1_list[i])  & (rsd1_corr==rsd2_list[i]) & (rsd2_corr==rsd1_list[i]) )
  rows <- c(rows1,rows2)
  if (length(rows) > 0) {
    coef_corr[i] <- phosite_corr$coef[rows]
    fdr_corr[i] <- phosite_corr$fdr[rows]
    rsd1[i] <- RSD1_list[rows]
    rsd2[i] <- RSD2_list[rows]
  }
}
pairwise$coef_corr <- coef_corr
pairwise$fdr_corr <- fdr_corr
pairwise$rsd1 <- rsd1
pairwise$rsd2 <- rsd2
pairwise$pair <- paste(gene1_list,rsd1,rsd2, sep = ":")

## write out tables
tn = paste(baseD,"pan3can_shared_data/analysis_results/hotspot3d/table/",cohort,"_phosphosite_within_protein_distance_and_correlation.txt", sep="")
write.table(pairwise, file=tn, quote=F, sep = '\t', row.names = FALSE)

pairwise_pos <- pairwise[!is.na(pairwise$fdr_corr) & pairwise$coef_corr > 0, ]
pairwise_neg <- pairwise[!is.na(pairwise$fdr_corr) & pairwise$coef_corr < 0, ]




