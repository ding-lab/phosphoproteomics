# adopted from Yige Wu @ WashU 2016 Dec
# Kuan Huang @ Washu 2017 Mar
# draw volcano plots for regression result

# choose kinase/phosphotase, cancer , significance level, model -----------------------------------------------
protein <- "kinase"
# protein <- "phosphotase"
sig <- 0.05
#val_sig = 0.1

# library -----------------------------------------------------------------
library(stringr)
library(ggplot2)
library(reshape)
library(grid)
require(plyr)

# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"

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


# # input regression processed data -----------------------------------------
# # manning kinome allow family-wise analysis
# manning_kinome = read.table(header=TRUE, quote = "", sep="\t",file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Gene_family/knome_render_peerj-01-126-s001_wgene.txt")
# manning_kinome_map = manning_kinome[,c(15,9,10)]
# colnames(manning_kinome_map)[1] = "KINASE"
# manning_kinome_map = manning_kinome_map[!duplicated(manning_kinome_map$KINASE),]

# table_HUMAN <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/",protein,"_substrate_regression_trans_edited.txt",sep = ""))
table_HUMAN_cis = read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_cis.txt",sep = ""))
table_HUMAN_trans = read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_trans.txt",sep = ""))

table_PDX_cis = read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_PDX_", protein,"_substrate_regression_cis.txt",sep = ""))
table_PDX_trans = read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_PDX_", protein,"_substrate_regression_trans.txt",sep = ""))
table_PDX_cis_sig = table_PDX_cis[table_PDX_cis$P_pro_kin < sig & table_PDX_cis$coef_pro_kin > 0,]
table_PDX_trans_sig = table_PDX_trans[table_PDX_trans$P_pho_kin < sig & table_PDX_trans$coef_pho_kin > 0,]

table_HUMAN_cis$validatedInPDX = (table_HUMAN_cis$pair %in% table_PDX_cis_sig$pair) & (table_HUMAN_cis$FDR_pro_kin < sig) 
table_HUMAN_trans$validatedInPDX = (table_HUMAN_trans$pair %in% table_PDX_trans_sig$pair) & (table_HUMAN_trans$FDR_pho_kin < sig) 
table_HUMAN_cis$sig = table_HUMAN_cis$FDR_pro_kin < sig & table_HUMAN_cis$coef_pro_kin > 0
table_HUMAN_trans$sig = table_HUMAN_trans$FDR_pho_kin < sig & table_HUMAN_trans$coef_pho_kin > 0

table_PDX_cis_sig_inHuman = table_PDX_cis_sig[table_PDX_cis_sig$pair %in% table_HUMAN_cis$pair[table_HUMAN_cis$sig],]
table_PDX_trans_sig_inHuman = table_PDX_trans_sig[table_PDX_trans_sig$pair %in% table_HUMAN_trans$pair[table_HUMAN_trans$sig],]
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_PDX_", protein,"_substrate_regression_cis_sig_inHuman_fam.txt",sep = "")
write.table(table_PDX_cis_sig_inHuman, file=tn, quote=F, sep = '\t', row.names = FALSE)
tn = paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_PDX_", protein,"_substrate_regression_trans_sig_inHuman_fam.txt",sep = "")
write.table(table_PDX_trans_sig_inHuman, file=tn, quote=F, sep = '\t', row.names = FALSE)

table(table_HUMAN_cis$validatedInPDX,table_HUMAN_cis$sig)
table(table_HUMAN_trans$validatedInPDX,table_HUMAN_trans$sig)

cat("\nPDX validated count\n")
table(table_PDX_cis_sig_inHuman$KINASE)[order(table(table_PDX_cis_sig_inHuman$KINASE),decreasing = T)][1:20]
table(table_PDX_trans_sig_inHuman$KINASE)[order(table(table_PDX_trans_sig_inHuman$KINASE),decreasing = T)][1:20]

##### volcano #####
# cis volcano plotting module -------------------------------------------------
p = ggplot(table_HUMAN_cis,aes(x=coef_pro_kin, y=-log10(FDR_pro_kin), color = validatedInPDX))
p = p + geom_point(alpha=0.3 , size=1, stroke = 0)
#p = p + geom_point(alpha=0.1 , stroke = 0 , size=0.5)
#p = p + geom_text(aes(label= ifelse((-log10(FDR_pro_kin)>9 | coef_pro_kin > 1.5) & sigInPDX, as.character(pair), NA)),size=1.2,alpha=0.8)
p = p + theme_bw() #+ theme_nogrid()
p = p + scale_color_manual(values = c("TRUE" = "red","FALSE" = "grey"))
p = p + geom_vline(xintercept = 0, alpha=0.3) + xlim(-2,2) + ylim(0,22)
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase protein expression", y="-log10(FDR)")
p = p + theme(legend.position="bottom")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_cis_inPDX_volcano.pdf',sep ="")
ggsave(file=fn, height=3, width=3, useDingbats=FALSE)

# trans volcano plotting module -------------------------------------------------
p = ggplot(table_HUMAN_trans,aes(x=coef_pho_kin, y=-log10(FDR_pho_kin), color =validatedInPDX))
p = p + geom_point(alpha=0.1 , size=0.5)
#p = p + geom_point(alpha=0.05 ,  stroke = 0 , size=0.5)
#p = p + geom_text(aes(label= ifelse((-log10(FDR_pho_kin)>4) & sigInPDX, as.character(pair), NA)),size=1.2,alpha=0.8)
p = p + theme_bw() #+ theme_nogrid()
p = p + scale_color_manual(values = c("TRUE" = "red","FALSE" = "grey"))
p = p + geom_vline(xintercept = 0, alpha=0.3) + xlim(-2,2) + ylim(0,9)
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for kinase phosphorylation level", y="-log10(FDR)")
p = p + theme(legend.position="bottom")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_trans_inPDX_volcano.pdf',sep ="")
ggsave(file=fn, height=3, width=3, useDingbats=FALSE)

# Size~signficance --------------------------------------------------------
p <- ggplot(table_HUMAN_trans, aes(x = -log10(FDR_pho_kin) , y = Size))
p <- p + geom_point(alpha=0.1, size = 0.2)
p = p + theme_bw()
p = p + labs(y = "sample size with data")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_trans_size_vs_fdr.pdf',sep ="")
ggsave(file=fn, height=3, width=3, useDingbats=FALSE)

p <- ggplot(table_HUMAN_cis, aes(x = -log10(FDR_pro_kin) , y = Size))
p <- p + geom_point(alpha=0.1, size = 0.2)
p = p + theme_bw()
p = p + labs(y = "sample size with data")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_cis_size_vs_fdr.pdf',sep ="")
ggsave(file=fn, height=3, width=3, useDingbats=FALSE)


##### barplots #####
# cis
table_HUMAN_cis_genes = data.frame(table(table_HUMAN_cis$KINASE))
table_HUMAN_cis_sig_genes = data.frame(table(table_HUMAN_cis$KINASE[table_HUMAN_cis$FDR_pro_kin < sig]))
table_cis_count = merge(table_HUMAN_cis_genes,table_HUMAN_cis_sig_genes,by="Var1")
colnames(table_cis_count) = c("KINASE","phosphosites","regulated")
table_cis_count$other = table_cis_count$phosphosites - table_cis_count$regulated
table_cis_count_f = table_cis_count[table_cis_count$regulated > 7,-2]

table_cis_count_m = melt(table_cis_count_f,id.vars="KINASE")
table_cis_count_m$KINASE = factor(table_cis_count_m$KINASE,levels=table_cis_count_f$KINASE[order(table_cis_count_f$regulated,decreasing=T)])

p <- ggplot()
p <- p + geom_bar(data=table_cis_count_m, aes(y = value, x = KINASE, fill = variable ), stat="identity",
                  position='stack')
p <- p + scale_fill_manual(values = c("regulated" = "red","other" = "grey"))
p <- p + theme_bw() #+ scale_y_log10()
p <- p + xlab(protein)+ylab("# of phosphosites")
p <- p + theme(axis.title=element_text(size=10))
p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p = p + theme(legend.position="bottom")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_cis_kinase_sites_barplot.pdf',sep ="")
ggsave(file=fn, height=4, width=3, useDingbats=FALSE)

table_cis_count_m = melt(table_cis_count[,-2],id.vars="KINASE")
table_cis_count_m$KINASE = factor(table_cis_count_m$KINASE,levels=table_cis_count$KINASE[order(table_cis_count$regulated,decreasing=T)])
p <- ggplot()
p <- p + geom_bar(data=table_cis_count_m, aes(y = value, x = KINASE, fill = variable ), stat="identity",
                  position='stack')
p <- p + scale_fill_manual(values = c("regulated" = "red","other" = "grey"))
p <- p + theme_nogrid() #+ scale_y_log10()
p = p + geom_vline(xintercept = dim(table_cis_count_f)[1] + 0.5, alpha=0.7)
p <- p + xlab(protein)+ylab("# of phosphosites")
p <- p + theme(axis.title=element_text(size=10), axis.text.x = element_blank(), axis.ticks.x= element_blank()) 
#p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p = p + theme(legend.position="bottom")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_cis_kinase_sites_barplot_all.pdf',sep ="")
ggsave(file=fn, height=3, width=6, useDingbats=FALSE)

table_cis_count$rate = table_cis_count$regulated/table_cis_count$phosphosites
table_cis_count = merge(table_cis_count,manning_kinome_map,by="KINASE",all.x=T)
table_cis_count$GroupName[is.na(table_cis_count$GroupName)] = "Other"
table_cis_count$log_value = log10(table_cis_count$regulated)
p <- ggplot(data=table_cis_count, aes(y = regulated, x = rate, color=GroupName ))
p <- p + geom_point(alpha = 0.2, stroke=0)
#p <- p + scale_fill_manual(values = c("regulated" = "red","other" = "grey"))
p <- p + theme_bw() #+ scale_y_log10(breaks=c(5,30,100,200,1000,10000,50000))
p <- p + xlab("% of associated sites")+ylab("# of phosphosites (10*)")
p = p + geom_text(aes(label= ifelse(regulated > 8, paste(as.character(KINASE),regulated, sep=": "), NA)),size=2)
p <- p + theme(axis.title=element_text(size=10))
p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p = p + theme(legend.position="right")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_cis_kinase_sites_count_ratio_plot.pdf',sep ="")
ggsave(file=fn, height=3, width=4, useDingbats=FALSE)

# trans
table_HUMAN_trans_genes = data.frame(table(table_HUMAN_trans$KINASE))
table_HUMAN_trans_sig_genes = data.frame(table(table_HUMAN_trans$KINASE[table_HUMAN_trans$FDR_pho_kin < sig]))
table_trans_count = merge(table_HUMAN_trans_genes,table_HUMAN_trans_sig_genes,by="Var1")
colnames(table_trans_count) = c("KINASE","phosphosites","regulated")
table_trans_count$other = table_trans_count$phosphosites - table_trans_count$regulated
table_trans_count$rate = table_trans_count$regulated/table_trans_count$phosphosites

table_trans_count_f = table_trans_count[table_trans_count$regulated > 30,-2]

table_trans_count_m = melt(table_trans_count_f,id.vars=c("KINASE","rate"))
table_trans_count_m$KINASE = factor(table_trans_count_m$KINASE,levels=table_trans_count_f$KINASE[order(table_trans_count_f$regulated,decreasing=T)])
table_trans_count_m$log_value = log10(table_trans_count_m$value)
table_trans_count_m$variable = factor(table_trans_count_m$variable,levels=c("other","regulated"))

p <- ggplot()
p <- p + geom_bar(data=table_trans_count_m, aes(y = log_value, x = KINASE, fill = variable ), stat="identity",
                  position='stack')
p <- p + scale_fill_manual(values = c("regulated" = "red","other" = "grey"))
p <- p + theme_bw() #+ scale_y_log10(breaks=c(5,30,100,200,1000,10000,50000))
p <- p + xlab(protein)+ylab("# of phosphosites (10*)")
p <- p + theme(axis.title=element_text(size=10))
p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p = p + theme(legend.position="bottom")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_trans_kinase_sites_barplot.pdf',sep ="")
ggsave(file=fn, height=4, width=3, useDingbats=FALSE)

table_trans_count_m = melt(table_trans_count[,-c(2,5)],id.vars="KINASE")
#table_trans_count_m$KINASE = factor(table_trans_count_m$KINASE,levels=table_trans_count$KINASE[order(table_trans_count$other,decreasing=T)])
table_trans_count_m$KINASE = factor(table_trans_count_m$KINASE,levels=table_trans_count$KINASE[order(table_trans_count$regulated,decreasing=T)])
table_trans_count_m$log_value = log10(table_trans_count_m$value)
table_trans_count_m$log_value[!is.finite(table_trans_count_m$log_value)] = 0
p <- ggplot()
p <- p + geom_bar(data=table_trans_count_m, aes(y = log_value, x = KINASE, fill = variable ), stat="identity",
                  position='stack')
p <- p + scale_fill_manual(values = c("regulated" = "red","other" = "grey"))
p <- p + theme_nogrid() #+ scale_y_log10()
p = p + geom_vline(xintercept = dim(table_trans_count_f)[1] + 0.5, alpha=0.7)
p <- p + xlab(protein)+ylab("# of phosphosites (10*)")
p <- p + theme(axis.title=element_text(size=10), axis.text.x = element_blank(), axis.ticks.x= element_blank()) 
#p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p = p + theme(legend.position="bottom")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_trans_kinase_sites_barplot_all.pdf',sep ="")
ggsave(file=fn, height=3, width=6, useDingbats=FALSE)

table_trans_count = merge(table_trans_count,manning_kinome_map,by="KINASE",all.x=T)
table_trans_count$GroupName[is.na(table_trans_count$GroupName)] = "Other"
table_trans_count$log_value = log10(table_trans_count$regulated)
p <- ggplot(data=table_trans_count, aes(y = log_value, x = rate, color=GroupName ))
p <- p + geom_point(alpha = 0.2, stroke=0)
#p <- p + scale_fill_manual(values = c("regulated" = "red","other" = "grey"))
p <- p + theme_bw() #+ scale_y_log10(breaks=c(5,30,100,200,1000,10000,50000))
p <- p + xlab("% of associated substrate sites")+ylab("# of phosphosites (10*)")
p = p + geom_text(aes(label= ifelse(log_value > 1.7 | (rate > 0.4 & log_value > 0.2), paste(as.character(KINASE),regulated, sep=": "), NA)),size=2)
p <- p + theme(axis.title=element_text(size=10)) + xlim(0,0.25)
p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
p = p + theme(legend.position="right")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_trans_kinase_sites_count_ratio_plot.pdf',sep ="")
ggsave(file=fn, height=3, width=4, useDingbats=FALSE)


protein = "phosphotase"
table_HUMAN_trans = read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/BRCA_HUMAN_", protein,"_substrate_regression_trans.txt",sep = ""))
# trans volcano plotting module -------------------------------------------------
p = ggplot(table_HUMAN_trans,aes(x=coef_pho_kin, y=-log10(FDR_pho_kin)))
p = p + geom_point(alpha=0.2 , size=0.5)
#p = p + geom_point(alpha=0.05 ,  stroke = 0 , size=0.5)
#p = p + geom_text(aes(label= ifelse((-log10(FDR_pho_kin)>4) & sigInPDX, as.character(pair), NA)),size=1.2,alpha=0.8)
p = p + theme_bw() #+ theme_nogrid()
p = p + geom_vline(xintercept = 0, alpha=0.3) + xlim(-2,2)
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "Coefficient for phosphotase phosphorylation level", y="-log10(FDR)")
p = p + theme(legend.position="bottom")
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/figures/BRCA_HUMAN_trans_phosphotase_volcano.pdf',sep ="")
ggsave(file=fn, height=3, width=3, useDingbats=FALSE)