# library -----------------------------------------------------------------
library(ggplot2)
library(grid)
library(dplyr)
library(stringr)
library(gridExtra)
library(gtable)
library(readr)

# # for working on Kuan's mac
# baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# for working on Yige's mac
baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# choose kinase or phosphotase, significance level, outlier threshold and least sample number-------------------------
sig <- 0.05 # significance level
protein <- "kinase"
#protein <- "phosphotase"
# cancer <- "BRCA" #choose in which cancer extract the top n rows
# cancer <- "OV"

# input -------------------------------------------------------------------
pho_data = read.delim(paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA_wGpos.txt",sep=""))
pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))
pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")

# calculate mean phosphorylation level in subtypes ------------------------
clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))
pho_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(pho_data)),4), ncol=4, byrow=T))
colnames(pho_subtype_mean) <- c("Her2","LumA","LumB","Basal")
for (cohort in c("Her2","LumA","LumB","Basal")) {
  subtype_sample <- na.omit(colnames(clinical)[clinical[1,]==cohort])
  pho_subtype_mean[,cohort] <- rowMeans(pho_data[,c(subtype_sample)], na.rm = TRUE)
}

OV_pho_g = paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PHO_formatted_normalized_max10NA.txt",sep="")
pho_gdata = read.delim(OV_pho_g)
samples_ov <- colnames(pho_gdata[,-1])
pho_ov = read.delim(paste(baseD,"pan3can_shared_data/OV/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq_abbrev_normalized_max10NA.tsv",sep=""))
pho_rsd_ov <- data.frame(str_split_fixed(pho_ov$X, ":", 3))
pho_rsd_ov[,3] <- toupper(pho_rsd_ov[,3])
colnames(pho_rsd_ov) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")
for (i in 1:nrow(pho_ov)) {
  row <- which(pho_rsd_split$SUBSTRATE==as.character(pho_rsd_ov$SUBSTRATE[i]) & pho_rsd_split$SUB_MOD_RSD==as.character(pho_rsd_ov$SUB_MOD_RSD[i]))
  if (length(row) > 0){
    pho_subtype_mean[row,"OV"] <- rowMeans(pho_ov[i,samples_ov], na.rm = TRUE)
  }
}

# calculate mean protein expression level in subtypes ------------------------
pro_data <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_PRO_formatted_normalized_max10NA.txt",sep=""))

pro_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(pro_data)),4), ncol=4, byrow=T))
colnames(pro_subtype_mean) <- c("Her2","LumA","LumB","Basal")
for (cohort in c("Her2","LumA","LumB","Basal")) {
  subtype_sample <- na.omit(colnames(clinical)[clinical[1,]==cohort])
  pro_subtype_mean[,cohort] <- rowMeans(pro_data[,c(subtype_sample)], na.rm = TRUE)
}
rownames(pro_subtype_mean) <- pro_data$X

pro_ov = read.delim(paste(baseD,"pan3can_shared_data/OV/OV_PNNL_PRO_formatted_normalized_max10NA.txt",sep=""))
colnames(pro_ov) <- str_split_fixed(colnames(pro_ov),"_",2)[,1]
for (i in 1:nrow(pro_ov)) {
  row <- which(pro_data$X==as.character(pro_ov$X[i]))
  if (length(row) > 0){
    pro_subtype_mean[row,"OV"] <- rowMeans(pro_ov[i,samples_ov], na.rm = TRUE)
  }
}

# input processed basal and OV data ---------------------------------------
table_subtypes <- read_delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_BRCA_subtype_trans_edited.txt", sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))
table_2can$cohort <- table_2can$Cancer
table_OV <- table_2can[table_2can$Cancer=="OV" ,colnames(table_subtypes)]
table_2can <- rbind(table_subtypes,table_OV)
table_sig <- table_2can[(!is.na(table_2can$FDR_pro_kin)) & table_2can$FDR_pro_kin<=sig |((!is.na(table_2can$FDR_pho_kin)) & table_2can$FDR_pho_kin<=sig),]

# bubble chart, top XX significant result ---------------------
# choose the parameters for bubble chart
for (iscis in c("cis","trans")) { # loop around cis and trans
  if (iscis == "cis") {
    var <- "pro_kin"
  }
  if (iscis == "trans") {
    var <- "pho_kin"
  }
  fdr_var <- paste("FDR_",var,sep = "")
  table_subtypes_sig <- table_subtypes[!is.na(table_subtypes[,fdr_var]) & table_subtypes[,fdr_var] <= sig,]
  table_OV_sig <- table_OV[!is.na(table_OV[,fdr_var]) & table_OV[,fdr_var] <= sig,]
  
  ks_overlap <- unique(merge(as.matrix(table_subtypes_sig[,c("KINASE","SUBSTRATE")]),
                             as.matrix(table_OV_sig[,c("KINASE","SUBSTRATE")])))
  k_overlap <- as.vector(ks_overlap$KINASE)
  s_overlap <- as.vector(ks_overlap$SUBSTRATE)
  
  rows <- c()
  for(i in 1:nrow(ks_overlap)){
    r <- which(table_sig$KINASE==k_overlap[i] & table_sig$SUBSTRATE==s_overlap[i])
    rows <- c(rows,r)
  }
  table_bubble <- table_sig[rows,c("KINASE","SUBSTRATE","SUB_MOD_RSD",fdr_var,paste("coef_",var,sep = ""),"SELF","cohort","pair")]
  
  if (nrow(table_bubble) > 0){
    colnames(table_bubble)[4:5] <- c("fdr","coef")
    table_bubble$sig <- ( table_bubble$fdr <= sig)
    
    ## bubble plot
    lim = max(abs(max(table_bubble$coef, na.rm = T)),abs(min(table_bubble$coef, na.rm = T)))
    p = ggplot(table_bubble,aes(x=SELF, y=pair))# make this the original ethni
    p = p + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
    p = p + geom_point(aes(fill=coef, size =-log10(fdr), color=ifelse(sig, "black",NA)),pch=21) 
    p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
    p = p + scale_colour_manual(values=c("black",NA))
    p = p + theme_bw() #+ theme_nogrid()
    p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
    #p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
    p = p + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
    p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
    p = p + theme(panel.margin.x=unit(0, "lines"), panel.margin.y=unit(0, "lines"))
    # p
    # fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/test.pdf',sep ="")
    # ggsave(file=fn, height=20, width=3.5, useDingbats=FALSE)
    
    # extract the subtype phosphorylation level according to bubble ch --------
    phosphosites <- unique(table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD")])
    pho_table <- c()
    for (cohort in c("Her2","LumA","LumB","Basal","OV")) {
      temp <- phosphosites
      temp$pair <- paste(temp$SUBSTRATE,temp$SUB_MOD_RSD,sep = ":")
      temp$cohort <- cohort
      temp$pho_subtype <- NA
      for (i in 1:nrow(temp)) {
        pho_temp <- pho_subtype_mean[pho_rsd_split$SUBSTRATE==as.character(temp$SUBSTRATE[i]) & pho_rsd_split$SUB_MOD_RSD==as.character(temp$SUB_MOD_RSD[i]),cohort]
        if (length(pho_temp) > 0) {
          temp$pho_subtype[i] <- pho_temp
        }
      }
      pho_table <- rbind(temp,pho_table)
    }
    pho_table$SELF <- iscis
    
    # plot heatmap for subtype phosphorylation level --------------------------
    lim = max(abs(max(pho_table$pho_subtype, na.rm = T)),abs(min(pho_table$pho_subtype, na.rm = T)))
    plot2 = ggplot(pho_table)
    plot2 = plot2 + geom_tile(aes(x=SELF, y=pair, fill=pho_subtype), color=NA)#, linetype="blank") 
    plot2 = plot2 + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
    plot2 = plot2 + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
    plot2 = plot2 + theme_bw() 
    # plot2 = plot2 + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=10))
    plot2 = plot2 + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
    plot2 = plot2 + theme(axis.ticks=element_blank(),legend.position="bottom")
    plot2 = plot2 + theme(panel.margin.x=unit(0, "lines"), panel.margin.y=unit(0, "lines"))
    
    # extract the subtype prosprorylation level according to bubble ch --------
    pro_table <- c()
    for (cohort in c("Her2","LumA","LumB","Basal","OV")) {
      temp <- phosphosites
      temp$cohort <- cohort
      temp$pro_subtype <- NA
      for (i in 1:nrow(temp)) {
        temp$pro_subtype[i] <- pro_subtype_mean[as.character(temp$KINASE[i]) ,cohort]
      }
      pro_table <- rbind(pro_table,temp)
    }
    pro_table$SELF <- iscis
    
    # plot heatmap for subtype kinase expression level --------------------------
    lim = max(abs(max(pro_table$pro_subtype)),abs(min(pro_table$pro_subtype)))
    plot3 = ggplot(pro_table)
    plot3 = plot3 + geom_tile(aes(x=SELF, y=SELF, fill=pro_subtype), color=NA)#, linetype="blank") 
    plot3 = plot3 + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
    plot3 = plot3 + scale_fill_gradientn(name= "pro_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
    plot3 = plot3 + theme_bw() 
    plot3 = plot3 + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
    plot3 = plot3 + theme(axis.ticks=element_blank(),legend.position="bottom")
    plot3 = plot3 + theme(panel.margin.x=unit(0, "lines"), panel.margin.y=unit(0, "lines"))
    
    # plot together --------------------------------------------------------------------
    fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,"/OV_overlaps_BRCA_subtypes_",iscis,"_",protein,'_table.pdf',sep ="")
    grid.newpage()
    # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
    pdf(fn, height = 20, width = 5.5)
    
    plot123 <- cbind(ggplotGrob(p), ggplotGrob(plot2), ggplotGrob(plot3), size = "last")
    title <- textGrob(paste(cis,"-regulated kinase-substrate pairs",sep = ""),gp=gpar(fontsize=16))
    padding <- unit(5,"mm")
    plottgt <- gtable_add_rows(plot123, 
                               heights = grobHeight(title) + padding,
                               pos = 0)
    plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
    grid.draw(plottgt)
    dev.off()
  }
}





