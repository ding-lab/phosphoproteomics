# Yige Wu @ WashU 2016 Nov
# mark protein phosphorylation score outlier samples with 
# mutation, protein expression outlier, upstream kinase phosphorylation score outlier

# choose kinase or phosphotase, outlier threshold and least sample number-------------------------
out_thres <- 1 #threshold for outlier
protein <- "kinase"
cancer <- "BRCA"

# library -----------------------------------------------------------------
library(reshape)
library(stringr)
library(readr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# input overlap file & k_s_table------------------------------------------------------
#overlap <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/kinase_BRCA_reverse_marking_out_thers_1.5.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
overlap <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/kinase_BRCA_reverse_marking_median_out_thers_1.txt",sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)

# initiate table for putting outlier statistics ---------------------------

overlap$is.phos.outlier <- (overlap$kin_phos_cscore>=out_thres) & (!is.na(overlap$kin_phos_cscore))
overlap$is.up.pho.outlier <- (overlap$up_phos_cscore>=out_thres) & (!is.na(overlap$up_phos_cscore))
overlap$is.sub.pho.outlier <- (overlap$sub_phos_ctransscore>=out_thres) & (!is.na(overlap$sub_phos_ctransscore))

phos_outlier_status <- data.frame(table(overlap[,c("kinase","is.phos.outlier")]))
phos_outlier_table <- phos_outlier_status[phos_outlier_status$is.phos.outlier==TRUE,]
phos_outlier_table <- phos_outlier_table[order(phos_outlier_table$Freq, decreasing = T),]

tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/No.phos.outlier_ranked_median_centered_",out_thres,"_IQR.txt", sep="")
write.table(phos_outlier_table, file=tn, quote=F, sep = '\t', row.names = FALSE)


# input cancer genes ------------------------------------------------------
cancer_genes_table <- read_delim(paste(baseD,"pan3can_shared_data/reference_files/Volgestin2013Science_125genes.txt", sep=""), 
                           "\t", escape_double = FALSE, col_types = cols(`Ocogene score*` = col_character()), 
                           trim_ws = TRUE)

# pull out for cancer genes -----------------------------------------------
rows <- c()
for ( gene in cancer_genes_table$`Gene Symbol`) {
  temp <- which(phos_outlier_table$kinase==gene)
  rows <- c(rows,temp)
}
phos_outlier_table4cancer <- phos_outlier_table[rows,]
phos_outlier_table4cancer <- phos_outlier_table4cancer[order(phos_outlier_table4cancer$Freq, decreasing = T),]
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/No.phos.outlier_ranked_median_centered_",out_thres,"_IQR_cancer_genes.txt", sep="")
write.table(phos_outlier_table4cancer, file=tn, quote=F, sep = '\t', row.names = FALSE)


# calculate outlier status ------------------------------------------------
is.phos.outlier <- (overlap$kin_phos_cscore>=out_thres) & (!is.na(overlap$kin_phos_cscore))
is.pro.outlier <- (overlap$pro_score>=out_thres) & (!is.na(overlap$pro_score))
is.up.pho.outlier <- (overlap$up_phos_cscore>=out_thres) & (!is.na(overlap$up_phos_cscore))
is.up.pro.outlier <- (overlap$up_pro_score>=out_thres) & (!is.na(overlap$up_pro_score))
is.cnv.outlier <- (overlap$cnv_score>=out_thres) & (!is.na(overlap$cnv_score))

#overlap$mut <- grepl("ins",overlap$mutation) | grepl("del",overlap$mutation) | grepl("sense",overlap$mutation) | grepl("splice",overlap$mutation)
overlap$outlier_status <- paste(is.phos.outlier,is.pro.outlier,overlap$mut,is.up.pho.outlier,sep = ":")
outlier_status <- data.frame(table(overlap[,c("kinase","outlier_status")]))


# input the gene list want to check out -----------------------------------
RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
rtk <- as.vector(t(RTK_file))

# check rtks --------------------------------------------------------------
gene_list <- rtk

row <- c()
for( gene in gene_list ) {
  row <- c(row,which(outlier_status$kinase==gene))
}
list_outlier_status <- outlier_status[row,]
outlier_status_split <- str_split_fixed(list_outlier_status$outlier_status,":",4)
colnames(list_outlier_status)[2] <- "pho.pro.mut.up"
p = ggplot(list_outlier_status[outlier_status_split[,1]=="TRUE",] ,aes(kinase,Freq))
p = p + geom_bar(aes(fill = pho.pro.mut.up), stat="identity") # , position='stack'
p = p + theme_bw()
p = p + xlab("protein") + ylab("number of samples")
p = p + coord_flip()
p = p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/BRCA_RTK_phos_outlier_split.pdf',sep ="")
ggsave(file=fn, height=10, width=13)


p = ggplot(list_outlier_status[list_outlier_status$pho.pro.mut.up!="FALSE:FALSE:FALSE:FALSE",] ,aes(kinase,Freq))
p = p + geom_bar(aes(fill = pho.pro.mut.up), stat="identity") # , position='stack'
p = p + theme_bw()
p = p + xlab("protein") + ylab("number of samples")
p = p + coord_flip()
p = p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p


# initiate table for putting outlier statistics adding upstream protein expression---------------------------
gene_list <- rtk
is.phos.outlier <- (overlap$kin_phos_cscore>=out_thres) & (!is.na(overlap$kin_phos_cscore))
is.pro.outlier <- (overlap$pro_score>=out_thres) & (!is.na(overlap$pro_score))
is.up.pho.outlier <- (overlap$up_phos_cscore>=out_thres) & (!is.na(overlap$up_phos_cscore))
is.up.pro.outlier <- (overlap$up_pro_score>=out_thres) & (!is.na(overlap$up_pro_score))

overlap$mut <- grepl("ins",overlap$mutation) | grepl("del",overlap$mutation) | grepl("sense",overlap$mutation) | grepl("splice",overlap$mutation)
overlap$outlier_status <- paste(is.phos.outlier,is.pro.outlier,overlap$mut,is.up.pho.outlier,is.up.pro.outlier,sep = ":")
outlier_status <- data.frame(table(overlap[,c("kinase","outlier_status")]))

row <- c()
for( gene in gene_list ) {
  row <- c(row,which(outlier_status$kinase==gene))
}
list_outlier_status <- outlier_status[row,]
outlier_status_split <- str_split_fixed(list_outlier_status$outlier_status,":",5)
colnames(list_outlier_status)[2] <- "pho.pro.mut.uppho.uppro"
p = ggplot(list_outlier_status[outlier_status_split[,1]=="TRUE",] ,aes(kinase,Freq))
p = p + geom_bar(aes(fill = pho.pro.mut.uppho.uppro), stat="identity") # , position='stack'
p = p + theme_bw()
p = p + xlab("protein") + ylab("number of samples")
p = p + coord_flip()
p = p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/BRCA_RTK_phos_outlier_split_4ways.pdf',sep ="")
ggsave(file=fn, height=10, width=13)


p = ggplot(list_outlier_status[list_outlier_status$pho.pro.mut.uppho.uppro!="FALSE:FALSE:FALSE:FALSE:FALSE",] ,aes(kinase,Freq))
p = p + geom_bar(aes(fill = pho.pro.mut.uppho.uppro), stat="identity") # , position='stack'
p = p + theme_bw()
p = p + xlab("protein") + ylab("number of samples")
p = p + coord_flip()
p = p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p



# make table for upset ----------------------------------------------------
overlap$is.phos.outlier <- (overlap$kin_phos_cscore>=out_thres) & (!is.na(overlap$kin_phos_cscore))
overlap$is.pro.outlier <- (overlap$pro_score>=out_thres) & (!is.na(overlap$pro_score))
overlap$is.up.pho.outlier <- (overlap$up_phos_cscore>=out_thres) & (!is.na(overlap$up_phos_cscore))
overlap$is.up.pro.outlier <- (overlap$up_pro_score>=out_thres) & (!is.na(overlap$up_pro_score))
overlap$is.cnv.outlier <- (overlap$cnv_score>=out_thres) & (!is.na(overlap$cnv_score))

row <- c()
for( gene in gene_list ) {
  row <- c(row,which(overlap$kinase==gene))
}
outlier_table <- data.frame(cbind(overlap[row,c("kinase","sample")],1*(overlap[row,c("is.up.pho.outlier", "is.up.pro.outlier", "is.pro.outlier", "is.cnv.outlier", "mut","is.phos.outlier")])))

# plot upset --------------------------------------------------------------
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/RTK_outliers_Upset.pdf',sep ="")
pdf(file=fn, height=5, width=7)
upset(outlier_table, sets.bar.color = "#56B4E9",
      sets = c("is.up.pho.outlier", "is.up.pro.outlier", "is.pro.outlier", "is.cnv.outlier", "mut", "is.phos.outlier"), 
      keep.order = T,
      intersections = list(list("is.phos.outlier","is.cnv.outlier","is.pro.outlier"),
                           list("is.phos.outlier","is.cnv.outlier"),
                           list("is.phos.outlier","mut","is.cnv.outlier","is.pro.outlier"),
                           list("is.phos.outlier","is.pro.outlier"),
                           list("is.phos.outlier"),
                           list("is.phos.outlier","mut"),
                           list("is.phos.outlier","is.up.pro.outlier"),
                           list("is.phos.outlier","is.up.pho.outlier"),
                           list("is.phos.outlier","is.cnv.outlier","is.up.pro.outlier"),
                           list("is.phos.outlier","is.cnv.outlier","is.pro.outlier","is.up.pro.outlier"),
                           list("is.phos.outlier","is.cnv.outlier","is.pro.outlier","is.up.pho.outlier")),
      order.by = "freq")
dev.off()

# initiate table for phospho outlier overlap with other outlier -----------
x <- vector(mode = "numeric", length = length(gene_list))
temp <- data.frame(matrix(rep(x,6), ncol = 6, byrow = T))
phos_outlier_split <- cbind(gene_list,temp)
colnames(phos_outlier_split) <- c("gene","pro.outlier","mutation","cnv","up.pro.outlier","up.pho.outlier","other")
rownames(phos_outlier_split) <- gene_list

# fill in the table -------------------------------------------------------
for( gene in gene_list ) {
  is.gene <- overlap$kinase==gene
  phos_outlier_split[gene,"pro.outlier"] <- length(which( is.gene & is.phos.outlier & is.pro.outlier))
  phos_outlier_split[gene,"mutation"] <- length(which(is.gene & is.phos.outlier & overlap$mut))
  phos_outlier_split[gene,"cnv"] <- length(which(is.gene & is.phos.outlier & is.cnv.outlier))
  phos_outlier_split[gene,"up.pro.outlier"] <- length(which(is.gene & is.phos.outlier & is.up.pro.outlier))
  phos_outlier_split[gene,"up.pho.outlier"] <- length(which(is.gene & is.phos.outlier & is.up.pho.outlier))
  phos_outlier_split[gene,"other"] <- length(which(is.gene & is.phos.outlier & !is.pro.outlier & !overlap$mut & !is.cnv.outlier & !is.up.pro.outlier & !is.up.pho.outlier))
}

# make table for bubble chart ---------------------------------------------
phos_outlier_split_m <- melt(phos_outlier_split, id="gene")
phos_outlier_split_m$value0 <- phos_outlier_split_m$value > 0


# order the variable ------------------------------------------------------
phos_outlier_split_m$order <- 6
phos_outlier_split_m$order[phos_outlier_split_m$variable=="cnv"] <- 5
phos_outlier_split_m$order[phos_outlier_split_m$variable=="pro.outlier"] <- 4
phos_outlier_split_m$order[phos_outlier_split_m$variable=="up.pro.outlier"] <- 3
phos_outlier_split_m$order[phos_outlier_split_m$variable=="up.pho.outlier"] <- 2
phos_outlier_split_m$order[phos_outlier_split_m$variable=="other"] <- 1

phos_outlier_split_m$outlier <- reorder(phos_outlier_split_m$variable,phos_outlier_split_m$order)

# plot bubble chart ------------------------------------------------------------
p = ggplot(phos_outlier_split_m[phos_outlier_split_m$value0,],aes(x=gene, y=outlier))# make this the original ethni
p = p + geom_point(aes(color=outlier, size = value ),pch=19) 
p = p + theme_bw()
p = p + labs(x = "gene", y="collapsed phosphorylation level outlier")
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/BRCA_RTK_phos_outlier_overlap.pdf',sep ="")
ggsave(file=fn, height=5, width=6)


# calculate correlation between kinase phosphorylation score and substrate phosphorylation score --------
corr_stat = try(cor.test(overlap$kin_phos_cscore,overlap$sub_phos_ctransscore, method = "pearson"), silent=T)
corr_stat$estimate
corr_stat$p.value

corr_stat = try(cor.test(overlap$kin_phos_cscore,overlap$up_phos_cscore, method = "pearson"), silent=T)
corr_stat$estimate
corr_stat$p.value

row <- c()
for( gene in rtk ) {
  row <- c(row,which(overlap$kinase==gene))
}
corr_stat = try(cor.test(overlap$kin_phos_cscore[row],overlap$sub_phos_ctransscore[row], method = "pearson"), silent=T)
corr_stat$estimate
corr_stat$p.value

corr_stat = try(cor.test(overlap$kin_phos_cscore[row],overlap$up_phos_cscore[row], method = "pearson"), silent=T)
corr_stat$estimate
corr_stat$p.value
