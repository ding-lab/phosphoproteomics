# Yige Wu @ WashU 2016 Dec
# scatterplot: number of kinase/substrates ~ correlatipn coefficients for number of kinase/substrate-based phosphorylation score

# choose kinase or phosphotase, outlier threshold and least sample number-------------------------
protein <- "kinase"
cancer <- "BRCA"
sig <- 0.05

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
# function -------------------------------------------------------------
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
# input correlation table and k_s_table -------------------------------------------------
# ks_phos_corr <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/reverse/table/",protein,"_",cancer,"_kin:sub-based_phos_score_correlate_protein_phos_score.txt", sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
ks_phos_corr <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_",cancer,"_validated_kin:sub-based_phos_score_correlate_protein_phos_score.txt", sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)

ks_phos_corr_s <- ks_phos_corr[!is.na(ks_phos_corr$s_estimate),c("gene","s_fdr","s_estimate","num_s")]; colnames(ks_phos_corr_s) <- c("gene","fdr","cor","num_ks"); ks_phos_corr_s$score_type <- "substrate-based"
ks_phos_corr_k <- ks_phos_corr[!is.na(ks_phos_corr$k_estimate),c("gene","k_fdr","k_estimate","num_k")]; colnames(ks_phos_corr_k) <- c("gene","fdr","cor","num_ks"); ks_phos_corr_k$score_type <- "kinase-based"

K_S_f = paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep="")
k_s_table = read.delim(K_S_f)


# plot number of substrate/kinase~correlation efficients ------------------
# ks_phos_corr_s$num_ks_filtered <- remove_outliers(ks_phos_corr_s$num_ks)
# ks_phos_corr_k$num_ks_filtered <- remove_outliers(ks_phos_corr_k$num_ks)
ks_phos_corr_table <- rbind(ks_phos_corr_s,ks_phos_corr_k)
ks_phos_corr_table$Cancer <- cancer

# ks_phos_corr_table_m <- ks_phos_corr_table[!is.na(ks_phos_corr_table$num_ks_filtered),]

p = ggplot(ks_phos_corr_table,aes(x=num_ks, y=cor))
p = p + facet_grid(Cancer~score_type,scales = "free_y")#, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
p = p + geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x)
p = p + geom_point(alpha=0.1 ,  stroke = 0 )
p = p + geom_text(aes(label= ifelse(cor>0.6, as.character(gene), NA ), vjust = 1, hjust = -0.2),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + xlim(0,50)
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "number of kinase/substrates", y="Correlation coefficients")
p
# fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/',protein,"_",cancer,'No.kin:sub~corr_kin:sub_phos_score.pdf',sep ="")
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/',protein,"_",cancer,'No.kin:sub~corr_kin:sub_phos_score_validated.pdf',sep ="")
ggsave(file=fn, height=5, width=6)


# barplot -----------------------------------------------------------------
# given the distribution of the number of kinases and substrates , better divide num_k to c(1:10, >10)
# correlation coefficients better divided into -1,-0.7,-0.5,-0.3,0,0.3,0.5,0.7,1
num_range <- c(0:5,Inf); num_name <- c(as.character(1:5),">5")
cor_range <- c(-1,-0.7,-0.5,-0.3,0,0.3,0.5,0.7,1)
for (i in 1:(length(num_range)-1)) {
  rown <- (ks_phos_corr_table$num_ks > num_range[i] & ks_phos_corr_table$num_ks <= num_range[i+1])
  ks_phos_corr_table$num_ks_range[rown] <- as.character(num_name[i])
}
for (j in 1:(length(cor_range)-1) ){
  rowc <- (ks_phos_corr_table$cor > cor_range[j] & ks_phos_corr_table$cor <= cor_range[j+1])
  ks_phos_corr_table$cor_range[rowc] <- paste(cor_range[j],"~",cor_range[j+1],sep ="")
}


p <- ggplot()
p <- p + geom_histogram(data=ks_phos_corr_table, aes(x = cor , fill = num_ks_range ),bins = 50,
                  position='stack')
p <- p + facet_grid(score_type~., drop=T,scales = "free")#, space = "free", scales = "free")
#p <- p + scale_fill_manual(values = c("pos_sig" = "red","none" = "grey"))
p <- p + theme_bw()
#p <- p + xlab("correlation coeffi")+ylab("number of substrate phosphosites")
p <- p + theme(axis.title=element_text(size=10))
#p <- p + coord_flip() 
p = p + geom_vline(xintercept=0)
p <- p + theme(axis.text.x = element_text(colour="black", size=8,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/',protein,"_",cancer,'No.kin:sub~corr_kin:sub_phos_score_hist.pdf',sep ="")
ggsave(file=fn, height=4, width=5)


# num_cor_table <- data.frame(table(ks_phos_corr_table[,c("score_type","num_ks_range","cor_range")]))
# cor_order <- as.numeric(str_split_fixed(num_cor_table$cor_range,"~",2)[,1])
# num_cor_table$Cor_range <- reorder(num_cor_table$cor_range,cor_order)
# 
# p <- ggplot()
# p <- p + geom_bar(data=num_cor_table, aes(y = Freq, x = Cor_range , fill = num_ks_range ), stat="identity",
#                   position='stack')
# p <- p + facet_grid(.~score_type, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
# #p <- p + scale_fill_manual(values = c("pos_sig" = "red","none" = "grey"))
# p <- p + theme_bw()
# #p <- p + xlab("correlation coeffi")+ylab("number of substrate phosphosites")
# p <- p + theme(axis.title=element_text(size=10))
# #p <- p + coord_flip() 
# p <- p + theme(axis.text.x = element_text(colour="black", size=8,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))
# p

# p <- ggplot()
# p <- p + geom_density(data=ks_phos_corr_table, aes(x = cor , fill = num_ks_range ),alpha=0.3)
# p <- p + facet_grid(score_type~., drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
# #p <- p + scale_fill_manual(values = c("pos_sig" = "red","none" = "grey"))
# p <- p + theme_bw()
# #p <- p + xlab("correlation coeffi")+ylab("number of substrate phosphosites")
# p <- p + theme(axis.title=element_text(size=10))
# #p <- p + coord_flip() 
# p = p + geom_vline(xintercept=0)
# p <- p + theme(axis.text.x = element_text(colour="black", size=8,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))
# p
