# Yige Wu @ WashU 2017 Jan
# plotting bubble chart and subtype phosphorylation/protein expression heatmap side by side for BRCA dataset

# library -----------------------------------------------------------------
library(ggplot2)
library(grid)
library(dplyr)
library(stringr)
library(gridExtra)
library(gtable)

# for working on Kuan's mac
baseD = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

# # for working on Yige's mac
# baseD = "/Users/yigewu/Box\ Sync/"


setwd(paste(baseD,"pan3can_analysis/phospho_network",sep=""))

source("../pan3can_aes.R") # aes for general purposes; it should be one directory out of the working directory

# choose kinase or phosphotase, significance level, outlier threshold and least sample number-------------------------
sig <- 0.05 # significance level
protein <- "kinase"
#protein <- "phosphotase"
top <- 50 # choose top n rows for fdr for model1
cancer <- "BRCA" #choose in which cancer extract the top n rows
# cancer <- "OV"

# input -------------------------------------------------------------------
pho_data = read.delim(paste(baseD,"pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized_max10NA_wGpos.txt",sep=""))
pho_rsd_split <- data.frame(str_split_fixed(pho_data$X, ":", 3))
pho_rsd_split[,3] <- toupper(pho_rsd_split[,3])
colnames(pho_rsd_split) <- c("SUBSTRATE","transcript","SUB_MOD_RSD")

#table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/regression/table/kinase_substrate_regression_trans_edited.txt",sep = ""))
table_2can <- read.delim(paste(baseD,"pan3can_shared_data/analysis_results/tables/",protein,"_substrate_regression_trans_edited.txt",sep = ""))

# calculate mean phosphorylation level in subtypes ------------------------
clinical <- read.delim(paste(baseD,"pan3can_shared_data/BRCA/BRCA_clinical_summary.txt",sep=""))
pho_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(pho_data)),4), ncol=4, byrow=T))
colnames(pho_subtype_mean) <- c("Her2","LumA","LumB","Basal")
for (cohort in c("Her2","LumA","LumB","Basal")) {
  subtype_sample <- na.omit(colnames(clinical)[clinical[1,]==cohort])
  pho_subtype_mean[,cohort] <- rowMeans(pho_data[,c(subtype_sample)], na.rm = TRUE)
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


# bubble chart, top XX significant result ---------------------
# choose the parameters for bubble chart
for (cancer in c("BRCA","OV")) {
  for (cis in c("cis","trans")) { # loop around cis and trans
    if (cis == "cis") {
      t0 <- table_2can[table_2can$SELF==cis & table_2can$FDR_pro_kin <= sig & table_2can$model=="pho_sub~pro_kin" & table_2can$Cancer==cancer,]
      var <- "pro_kin"
    }
    if (cis == "trans") {
      t0 <- table_2can[table_2can$SELF==cis & table_2can$FDR_pho_kin <= sig & table_2can$model=="pho_sub~pro_sub+pho_kin" & table_2can$Cancer==cancer,]
      var <- "pho_kin"
    }
    if (nrow(t0) > 0){
      fdr_var <- paste("FDR_",var,sep = "");  coef_var <- paste("coef_",var,sep = "")
      t1 <- t0[order(t0[,fdr_var]),]
      len <- min(top,nrow(t1))
      table_bubble <- t1[1:len,c("KINASE","SUBSTRATE","SUB_MOD_RSD",fdr_var, coef_var,"Cancer","pair")]
      colnames(table_bubble) <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","fdr","coef","Cancer","pair")
      table_bubble$sig <- (table_bubble$fdr <= sig)
      
      ## bubble plot
      lim = max(abs(max(table_bubble$coef, na.rm = T)),abs(min(table_bubble$coef, na.rm = T)))
      p = ggplot(table_bubble,aes(x=Cancer, y=pair))# make this the original ethni
      p = p + facet_grid(KINASE~Cancer, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
      p = p + geom_point(aes(fill=coef, size =-log10(fdr), color = NA),pch=21) 
      p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
      p = p + scale_colour_manual(values=c("black",NA))
      p = p + theme_bw() #+ theme_nogrid()
      p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
      #p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
      p = p + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
      p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
      p = p + theme(panel.margin.x=unit(0, "lines"), panel.margin.y=unit(0, "lines"))
      # p
      # fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/',c,'_sorted_cis_',cis,'_fdr',sig,'_top_',top,'.pdf',sep ="")
      # ggsave(file=fn, height=10, width=6, useDingbats=FALSE)
      
      # extract the subtype phosphorylation level according to bubble ch --------
      pho_table <- c()
      for (cohort in c("Her2","LumA","LumB","Basal")) {
        temp <- table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","Cancer")]
        temp$pair <- paste(temp$SUBSTRATE,temp$SUB_MOD_RSD,sep = ":")
        temp$cohort <- cohort
        temp$pho_subtype <- NA
        for (i in 1:nrow(table_bubble)) {
          pho_temp <- pho_subtype_mean[pho_rsd_split$SUBSTRATE==as.character(temp$SUBSTRATE[i]) & pho_rsd_split$SUB_MOD_RSD==as.character(temp$SUB_MOD_RSD[i]),cohort]
          if (length(pho_temp) > 0) {
            temp$pho_subtype[i] <- pho_temp
          }
        }
        pho_table <- rbind(pho_table,temp)
      }
      
      # plot heatmap for subtype phosphorylation level --------------------------
      lim = max(abs(max(pho_table$pho_subtype, na.rm = T)),abs(min(pho_table$pho_subtype, na.rm = T)))
      plot2 = ggplot(pho_table)
      plot2 = plot2 + geom_tile(aes(x=Cancer, y=pair, fill=pho_subtype), color=NA)#, linetype="blank") 
      plot2 = plot2 + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
      plot2 = plot2 + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
      plot2 = plot2 + theme_bw() 
      # plot2 = plot2 + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=10))
      plot2 = plot2 + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
      plot2 = plot2 + theme(axis.ticks=element_blank(),legend.position="bottom")
      plot2 = plot2 + theme(panel.margin.x=unit(0, "lines"), panel.margin.y=unit(0, "lines"))
      
      # extract the subtype prosprorylation level according to bubble ch --------
      pro_table <- c()
      for (cohort in c("Her2","LumA","LumB","Basal")) {
        temp <- table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","Cancer","pair")]
        temp$cohort <- cohort
        temp$pro_subtype <- NA
        for (i in 1:nrow(table_bubble)) {
          temp$pro_subtype[i] <- pro_subtype_mean[as.character(temp$KINASE[i]) ,cohort]
        }
        pro_table <- rbind(pro_table,temp)
      }
      
      
      # plot heatmap for subtype kinase expression level --------------------------
      lim = max(abs(max(pro_table$pro_subtype)),abs(min(pro_table$pro_subtype)))
      plot3 = ggplot(pro_table)
      plot3 = plot3 + geom_tile(aes(x=Cancer, y=pair, fill=pro_subtype), color=NA)#, linetype="blank") 
      plot3 = plot3 + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
      plot3 = plot3 + scale_fill_gradientn(name= "pro_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
      plot3 = plot3 + theme_bw() 
      plot3 = plot3 + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
      plot3 = plot3 + theme(axis.ticks=element_blank(),legend.position="bottom")
      plot3 = plot3 + theme(panel.margin.x=unit(0, "lines"), panel.margin.y=unit(0, "lines"))
      
      # plot together --------------------------------------------------------------------
      fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,"/",cancer,"/",cancer,"_",protein,'_subtype_phospho_and_expression_level_top',len,"_",cis,'.pdf',sep ="")
      grid.newpage()
      # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
      pdf(fn, height = 12*(len+10)/(top+10), width = 5.5)
      
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
}





# # bubble chart, appoint gene list ---------------------
# # choose the parameters for bubble chart
# gene_list <- c("PIK3CA", "PTEN", "AKT1", "TP53", "GATA3", "CDH1", "RB1", "MLL3", "MAP3K1", "CDKN1B","TBX3", "RUNX1", "CBFB", "AFF2", "PIK3R1", "PTPN22", "PTPRD", "NF1", "SF3B1", "CCND3","NF1", "BRCA1", "BRCA2", "RB1", "CDK12","RB1", "NF1", "FAT3", "CSMD3", "GABRA6", "CDK12")
# gene_list <- unique(gene_list)
# 
# RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
# gene_list <- as.vector(t(RTK_file))
# 
# # gene_list as substrate
# for (cis in c("trans","cis")) { # loop around cis and trans
#   t0 <- table_brca[table_brca$SELF==cis & table_brca$FDR_pro_kin <= sig,]
#   
#   ## corresponding results for other two models are extracted and ordered
#   rows <- c()
#   for(gene in gene_list){
#     r <- unlist(which(t0$SUBSTRATE==gene))
#     rows <- c(rows,r)
#   }
#   table_bubble <- t0[rows,]
#   table_bubble$sig <- (table_bubble$fdr <= sig)
#   table <- table_bubble[table_bubble$Cancer=="BRCA",]
#   #table_2can_brca_top = table_2can[table_2can$pair %in% table$pair,]
#   
#   ## actual plotting
#   lim = max(abs(max(table_bubble$coef, na.rm = T)),abs(min(table_bubble$coef, na.rm = T)))
#   p = ggplot(table,aes(x=Cancer, y=pair))# make this the original ethni
#   #  p = p + facet_grid(KINASE~., drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
#   p = p + facet_grid(KINASE~Cancer, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
#   p = p + geom_point(aes(fill=coef, size =-log10(fdr), color=ifelse(sig, "black",NA)),pch=21) 
#   p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
#   p = p + scale_colour_manual(values=c("black",NA))
#   p = p + theme_bw() #+ theme_nogrid()
#   p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   # p
#   # fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/',c,'_sorted_cis_',cis,'_fdr',sig,'_top_',top,'.pdf',sep ="")
#   # ggsave(file=fn, height=10, width=6, useDingbats=FALSE)
#   
#   # extract the subtype phosphorylation level according to bubble ch --------
#   pho_table <- c()
#   for (cohort in c("Her2","LumA","LumB","Basal")) {
#     temp <- table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","Cancer")]
#     temp$pair <- paste(temp$SUBSTRATE,temp$SUB_MOD_RSD,sep = ":")
#     temp$cohort <- cohort
#     temp$pho_subtype <- NA
#     for (i in 1:nrow(table)) {
#       temp$pho_subtype[i] <- pho_subtype_mean[pho_rsd_split$SUBSTRATE==as.character(temp$SUBSTRATE[i]) & pho_rsd_split$SUB_MOD_RSD==as.character(temp$SUB_MOD_RSD[i]),cohort]
#     }
#     pho_table <- rbind(pho_table,temp)
#   }
#   
#   
#   # plot heatmap for subtype phosphorylation level --------------------------
#   lim = max(abs(max(pho_table$pho_subtype, na.rm = T)),abs(min(pho_table$pho_subtype, na.rm = T)))
#   plot2 = ggplot(pho_table)
#   plot2 = plot2 + geom_tile(aes(x=Cancer, y=pair, fill=pho_subtype), color=NA)#, linetype="blank") 
#   plot2 = plot2 + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
#   plot2 = plot2 + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
#   plot2 = plot2 + theme_bw() 
#   plot2 = plot2 + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   
#   # extract the subtype prosprorylation level according to bubble ch --------
#   pro_table <- c()
#   for (cohort in c("Her2","LumA","LumB","Basal")) {
#     temp <- table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","Cancer","pair")]
#     temp$cohort <- cohort
#     temp$pro_subtype <- NA
#     for (i in 1:nrow(table)) {
#       temp$pro_subtype[i] <- pro_subtype_mean[as.character(temp$KINASE[i]) ,cohort]
#     }
#     pro_table <- rbind(pro_table,temp)
#   }
#   
#   
#   # plot heatmap for subtype kinase expression level --------------------------
#   lim = max(abs(max(pro_table$pro_subtype)),abs(min(pro_table$pro_subtype)))
#   plot3 = ggplot(pro_table)
#   plot3 = plot3 + geom_tile(aes(x=Cancer, y=pair, fill=pro_subtype), color=NA)#, linetype="blank") 
#   plot3 = plot3 + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
#   plot3 = plot3 + scale_fill_gradientn(name= "pro_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
#   plot3 = plot3 + theme_bw() 
#   plot3 = plot3 + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
#   
#   # plot together --------------------------------------------------------------------
#   fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,"/",protein,'_subtype_phospho_and_expression_level_sub_gene_list_',cis,'.pdf',sep ="")
#   grid.newpage()
#   pdf(fn, height = 12, width = 15)
#   plot123 <- cbind(ggplotGrob(p), ggplotGrob(plot2), ggplotGrob(plot3), size = "last")
#   title <- textGrob(paste(cis,"-regulated kinase-substrate pairs",sep = ""),gp=gpar(fontsize=16))
#   padding <- unit(5,"mm")
#   plottgt <- gtable_add_rows(plot123, 
#                              heights = grobHeight(title) + padding,
#                              pos = 0)
#   plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
#   grid.draw(plottgt)
#   dev.off()
# }
# 
# # gene_list as substrate
# for (cis in c("trans","cis")) { # loop around cis and trans
#   t0 <- table_brca[table_brca$SELF==cis & table_brca$FDR_pro_kin <= sig,]
#   
#   ## corresponding results for other two models are extracted and ordered
#   rows <- c()
#   for(gene in gene_list){
#     r <- unlist(which(t0$KINASE==gene))
#     rows <- c(rows,r)
#   }
#   table_bubble <- t0[rows,]
#   table_bubble$sig <- (table_bubble$fdr <= sig)
#   table <- table_bubble[table_bubble$Cancer=="BRCA",]
#   #table_2can_brca_top = table_2can[table_2can$pair %in% table$pair,]
#   
#   ## actual plotting
#   lim = max(abs(max(table_bubble$coef, na.rm = T)),abs(min(table_bubble$coef, na.rm = T)))
#   p = ggplot(table,aes(x=Cancer, y=pair))# make this the original ethni
#   #  p = p + facet_grid(KINASE~., drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
#   p = p + facet_grid(KINASE~Cancer, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
#   p = p + geom_point(aes(fill=coef, size =-log10(fdr), color=ifelse(sig, "black",NA)),pch=21) 
#   p = p + scale_fill_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
#   p = p + scale_colour_manual(values=c("black",NA))
#   p = p + theme_bw() #+ theme_nogrid()
#   p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   # p
#   # fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/',c,'_sorted_cis_',cis,'_fdr',sig,'_top_',top,'.pdf',sep ="")
#   # ggsave(file=fn, height=10, width=6, useDingbats=FALSE)
#   
#   # extract the subtype phosphorylation level according to bubble ch --------
#   pho_table <- c()
#   for (cohort in c("Her2","LumA","LumB","Basal")) {
#     temp <- table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","Cancer")]
#     temp$pair <- paste(temp$SUBSTRATE,temp$SUB_MOD_RSD,sep = ":")
#     temp$cohort <- cohort
#     temp$pho_subtype <- NA
#     for (i in 1:nrow(table)) {
#       temp$pho_subtype[i] <- pho_subtype_mean[pho_rsd_split$SUBSTRATE==as.character(temp$SUBSTRATE[i]) & pho_rsd_split$SUB_MOD_RSD==as.character(temp$SUB_MOD_RSD[i]),cohort]
#     }
#     pho_table <- rbind(pho_table,temp)
#   }
#   
#   
#   # plot heatmap for subtype phosphorylation level --------------------------
#   lim = max(abs(max(pho_table$pho_subtype, na.rm = T)),abs(min(pho_table$pho_subtype, na.rm = T)))
#   plot2 = ggplot(pho_table)
#   plot2 = plot2 + geom_tile(aes(x=Cancer, y=pair, fill=pho_subtype), color=NA)#, linetype="blank") 
#   plot2 = plot2 + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
#   plot2 = plot2 + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
#   plot2 = plot2 + theme_bw() 
#   plot2 = plot2 + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   
#   # extract the subtype prosprorylation level according to bubble ch --------
#   pro_table <- c()
#   for (cohort in c("Her2","LumA","LumB","Basal")) {
#     temp <- table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","Cancer","pair")]
#     temp$cohort <- cohort
#     temp$pro_subtype <- NA
#     for (i in 1:nrow(table)) {
#       temp$pro_subtype[i] <- pro_subtype_mean[as.character(temp$KINASE[i]) ,cohort]
#     }
#     pro_table <- rbind(pro_table,temp)
#   }
#   
#   
#   # plot heatmap for subtype kinase expression level --------------------------
#   lim = max(abs(max(pro_table$pro_subtype)),abs(min(pro_table$pro_subtype)))
#   plot3 = ggplot(pro_table)
#   plot3 = plot3 + geom_tile(aes(x=Cancer, y=pair, fill=pro_subtype), color=NA)#, linetype="blank") 
#   plot3 = plot3 + facet_grid(KINASE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
#   plot3 = plot3 + scale_fill_gradientn(name= "pro_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
#   plot3 = plot3 + theme_bw() 
#   plot3 = plot3 + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
#   
#   # plot together --------------------------------------------------------------------
#   fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,"/",protein,'_subtype_phospho_and_expression_level_kin_gene_list_',cis,'.pdf',sep ="")
#   grid.newpage()
#   pdf(fn, height = 12, width = 15)
#   plot123 <- cbind(ggplotGrob(p), ggplotGrob(plot2), ggplotGrob(plot3), size = "last")
#   title <- textGrob(paste(cis,"-regulated kinase-substrate pairs",sep = ""),gp=gpar(fontsize=16))
#   padding <- unit(5,"mm")
#   plottgt <- gtable_add_rows(plot123, 
#                              heights = grobHeight(title) + padding,
#                              pos = 0)
#   plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
#   grid.draw(plottgt)
#   dev.off()
# }
# 
