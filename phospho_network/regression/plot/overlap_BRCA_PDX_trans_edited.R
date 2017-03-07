# Yige Wu @ WashU 2016 Nov
# look at overlap between BRCA data and PDX data, positive coefficient

# library -----------------------------------------------------------------
library(ggplot2)
library(grid)
library(dplyr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory


# integrate BRCA & FDX processed data -------------------------------------
protein <- "kinase"
sig <- 0.05
table_PDX <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_PDX_substrate_regression_trans_edited.txt",sep = ""))
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))

time <- Sys.time()
sink(file = paste(baseD,"pan3can_shared_data/analysis_results/tables/No.overlap_bw_PDX_BRCAorOV_fdr_",sig,"_",time,".txt",sep = ""))
for (cancer in c("BRCA","OV")) {
  t0 <- rbind(table_PDX,table_2can[table_2can$Cancer==cancer,]) 
  for (iscis in c("trans","cis")) { # loop around cis and trans
    
    if (iscis == "cis") {
      var <- "pro_kin"
    }
    if (iscis == "trans") {
      var <- "pho_kin"
    }
    fdr_var <- paste("FDR_",var,sep = "")
    coef_var <- paste("coef_",var,sep = "")
    
    table_PDX_sig <- table_PDX[table_PDX$SELF == iscis & table_PDX[,fdr_var] <= sig & table_PDX[,coef_var] > 0,]
    table_can_sig <- table_2can[table_2can$Cancer == cancer & table_2can$SELF == iscis & table_2can[,fdr_var] <= sig & table_2can[,coef_var] > 0,]
    
    ks_overlap <- unique(merge(as.matrix(table_PDX_sig[,c("KINASE","SUBSTRATE")]),
                               as.matrix(table_can_sig[,c("KINASE","SUBSTRATE")])))
    # k_overlap <- as.vector(ks_overlap$KINASE)
    # s_overlap <- as.vector(ks_overlap$SUBSTRATE)
    
    ksp_overlap <- unique(merge(as.matrix(table_PDX_sig[,c("KINASE","SUBSTRATE","SUB_MOD_RSD")]),
                                as.matrix(table_can_sig[,c("KINASE","SUBSTRATE","SUB_MOD_RSD")])))
    k_overlap <- as.vector(ksp_overlap$KINASE)
    s_overlap <- as.vector(ksp_overlap$SUBSTRATE)
    rsd_overlap <- as.vector(ksp_overlap$SUB_MOD_RSD)
    
    cat(paste("PDX overlap with ",cancer," ",iscis," significant regression result (FDR <= ",sig,", coef > 0):\n",
              nrow(ks_overlap)," kinase:substrate pairs (",
              nrow(unique(table_PDX_sig[,c("KINASE","SUBSTRATE")]))," PDX pairs vs ",
              nrow(unique(table_can_sig[,c("KINASE","SUBSTRATE")]))," BRCA pairs)\n",
              
              nrow(ksp_overlap)," kinase:substrate:phosphosite pairs (",
              nrow(table_PDX_sig)," PDX pairs vs ",
              nrow(table_can_sig)," BRCA pairs)\n\n",
              sep = ""))
    
    
    rows <- c()
    for(i in 1:nrow(ks_overlap)){
      r <- which(t0$KINASE==k_overlap[i] & t0$SUBSTRATE==s_overlap[i] & t0$SUB_MOD_RSD==rsd_overlap[i])
      rows <- c(rows,r)
    }
    table <- t0[rows,c("KINASE","SUBSTRATE","SUB_MOD_RSD",fdr_var,paste("coef_",var,sep = ""),"SELF","Cancer","pair")]
    colnames(table)[4:5] <- c("fdr","coef")
    table$sig <- (table$fdr <= sig)
    
    if (nrow(table) > 0){
      lim = max(abs(max(table$coef)),abs(min(table$coef)))
      p = ggplot(table,aes(x=SELF, y=pair))# make this the original ethni
      #  p = p + facet_grid(KINASE~., drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
      p = p + facet_grid(KINASE~Cancer, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
      p = p + geom_point(aes(fill=coef, size =-log10(fdr), color=ifelse(sig, "black",NA)),pch=21) 
      p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
      p = p + scale_colour_manual(values=c("black",NA))
      p = p + theme_bw() #+ theme_nogrid()
      p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
      # p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
      p = p + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
      p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
      p = p + theme(panel.margin.x=unit(0, "lines"), panel.margin.y=unit(0, "lines"))
      p
      fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/PDX_',cancer,'_overlap_',iscis,'_sig_',sig,'.pdf',sep ="")
      if (cancer == "BRCA") {
        ggsave(file=fn, height=10, width=3, useDingbats=FALSE)
      }
      if (cancer == "OV") {
        ggsave(file=fn, height=4.5, width=2.5, useDingbats=FALSE)
      }
    }
  }
}
sink()



# qqplot: -log10(FDR) of BRCA and PDX -------------------------------------
# plot_fdr_scale <- -log10(0.05)
overlap_pos$top <- FALSE
temp <- overlap_pos[order(overlap_pos$scale_FDR.y, decreasing = TRUE),]
top <- 100
for (self in c("cis","trans")) {
  for (cohort in c("BRCA","OV")) {
    fdr_order <- temp[temp$SELF==self & temp$cohort.x==cohort,]
    rows <- overlap_pos$SELF==self & overlap_pos$cohort.x==cohort
    overlap_pos$top[rows] <- overlap_pos$scale_FDR.y[rows] >= fdr_order$scale_FDR.y[min(top,nrow(fdr_order))]
  }
}

p = ggplot(overlap_pos,aes(x=scale_FDR.x, y=scale_FDR.y))
p = p + facet_grid(SELF~cohort.x,scales = "free_y")
p = p + geom_point(alpha=0.05)
p = p + geom_text(aes(label= ifelse( (SELF=="cis" & scale_FDR.x >= -log10(0.0005) & top) | (SELF=="trans" & scale_FDR.x >= -log10(0.01) & top) , pair, NA)),size=2,alpha=0.5)
#p = p + geom_text(aes(label= ifelse( top & scale_FDR.x >= -log10(0.05) , pair, NA)),size=2,alpha=0.5)
p = p + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p = p + labs(x = "-log(FDR)", y="-log(FDR) in PDX")
p = p + xlim(0,10) #+ ylim(0,10)
p
fn = paste(baseD,'pan3can_shared_data/analysis_results/PDX/overlap/BRCA_OV_regression_overlap_PDX_both_coef_positive_top_',top,'PDX.pdf',sep ="")
ggsave(file=fn, height=8, width=8, useDingbats=FALSE)

