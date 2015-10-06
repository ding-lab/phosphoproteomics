# outlier analysis pipeline
# nice note for outlier analysis: http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.html

### dependencies ###
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/outlier")
source("/Users/khuang/bin/LIB_exp.R")
system("mkdir logs")
logFile = paste("logs/", date, "_outlier_analysis.log", sep="")
sink(file=logFile)

##### original version #####
# rewriting to look at sample-wise z-score distribution to define outliers
# function to find outliers using the modified zscores. Take a matrix as input
# requires function find_kde_inflection
# requires library reshape
# identify the adjacent local minimum as the threshold to define outliers
# original outlier: using second derivative tests
# find_outlier_s = function(m, name, genes=NULL, minNum = 10, filter=TRUE, plot=FALSE, h=6, w=10){ 
#   print("##### OUTLIER ANALYSIS #####")
#   num = nrow(m)
#   if (filter){
#     m2 = m[rowSums(!is.na(m)) >= minNum, ]
#     m2_genes = m2[row.names(m2) %in% genes, ]
#   } else{ m2 = m; m2_genes=m}
#   num_NA = nrow(m2)
#   num_gene = nrow(m2_genes)
#   m2_genes = as.matrix(m2_genes)
#   print(paste("Looking for outliers in", deparse(substitute(genes)), "of", deparse(substitute(m)), sep=" "))
#   print(paste("Original number of genes:", num, "; NA filtered:", num_NA, "; Gene list filtered:", num_gene, sep=" "))
# 
#   kde_dir_cmd = paste("mkdir figures/", date, "/KDE/", sep="")
#   system(kde_dir_cmd)
#   outlier = matrix(,nrow=dim(m2_genes)[1],ncol=dim(m2_genes)[2])
#   row.names(outlier) = row.names(m2_genes)
#   colnames(outlier) = colnames(m2_genes)
#   outlier_zscore = outlier
#   
#   # gene-wise zscore
#   for (i in 1:nrow(m2_genes)){
#     m = median(m2_genes[i,], na.rm=TRUE)
#     rowMAD=mad(m2_genes[i,], na.rm=TRUE)
#     # modified z-score for outlier: Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
#     outlier_zscore[i,]  = 0.6745*(m2_genes[i,]-m)/rowMAD
#   }
#   
#   # sample-wise outlier
#   for (i in 1:ncol(outlier_zscore)) { 
#     sample = colnames(outlier_zscore)[i]
#     outlier_zscore_s = outlier_zscore[,i,drop=F]
#     
#     #     library(mixtools)
#     #     mixMdl = normalmixEM(outlier_zscore_s, k=4)
#     #     plot(mixMdl,which=2)
#     #     library(fpc)
#     #     pamk(outlier_zscore_s, krange=2:15)
#     kde=find_kde_criticals(outlier_zscore_s)
#     maxs=kde$maxs
#     mins=kde$mins
#     adj_max=kde$adj_max
#     adj_min=kde$adj_min
#     
#     
#     # make sure the adj_min threshold is greater than 95% quantile
#     per95 = quantile(outlier_zscore_s, na.rm=T, probs=0.95)
#     if (!is.na(adj_min)){ if (adj_min < per95){adj_min=mins[mins>per95][1]}}
#     
#     # make a high threshold if adjacent min is missing
#     if (is.na(adj_min)){adj_min = kde$dmode + 3*mad(outlier_zscore_s, na.rm=T)
#       mins = c(mins,adj_min)
#     }
#     
#     ## plot ## 
#     if (plot){
#       df = as.data.frame(outlier_zscore_s)
#       fn = paste("figures/", date, "/KDE/",date, "_", sample, "_", name, '_KDE_mzscore_distribution.pdf', sep = "")
#       p = ggplot(data=df, aes_string(x=sample)) + geom_density()
#       p = p + geom_vline(xintercept=maxs, colour="red", alpha=0.2) + geom_vline(xintercept=mins, colour="blue", alpha=0.2) + geom_vline(xintercept=adj_min, colour="blue")
#       p = p + xlab(paste("modified Z-score in", name, sep=" ")) + 
#         scale_x_continuous(limits=c(-5,5),breaks=c(mins,maxs)) + theme_bw() + 
#         theme(text = element_text(colour="black", size=16), axis.text.x = element_text(angle = 90,colour="black", size=14), axis.text.y = element_text(colour="black", size=14))
#       p
#       ggsave(file=fn)
#     }
#     outlier[,i] = outlier_zscore[,i] >= adj_min
#     print(paste(sample, name, "KDE_outlier_threshold:", adj_min, "Number_of_outliers", sum(outlier[,i], na.rm=T), sep = "  "))
#   }
#   #return(list("outlier_zscore"=outlier_zscore, "outlier"=outlier))
#   ### get the returned list from function find_outlier() and find top 5% z-score and genes within each sample
# 
#   # set up return matrixes
#   zscore=outlier_zscore
#   outlier=outlier
#   num_5percent = dim(zscore)[1]
#   top_outlier_zscore = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
#   top_outlier = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
#   top_outlier_boolean = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
#   top_outlier_raw = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
#   row.names(top_outlier_zscore)=colnames(zscore)
#   row.names(top_outlier)=colnames(zscore)
#   row.names(top_outlier_boolean)=colnames(zscore)
#   row.names(top_outlier_raw)=colnames(zscore)
#   colnames(top_outlier_zscore)=c(1:num_5percent)
#   for (i in 1:num_5percent){colnames(top_outlier_zscore)[i] = paste(name, "top", colnames(top_outlier_zscore)[i], sep=" ")}
#   colnames(top_outlier)=colnames(top_outlier_zscore)
#   colnames(top_outlier_boolean)=colnames(top_outlier_zscore)
#   colnames(top_outlier_raw)= colnames(top_outlier_zscore)
#   # find top 5
#   for (i in colnames(zscore)){
#     whim=zscore[,i]
#     a = whim[order(whim, decreasing=TRUE)][1:num_5percent]
#     top_outlier_zscore[i,] = a
#     whim2 = outlier[,i]
#     top_outlier_boolean[i,] = whim2[order(whim, decreasing=TRUE)][1:num_5percent]
#     whim3 = m2_genes[,i]
#     top_outlier_raw[i,] = whim3[order(whim, decreasing=TRUE)][1:num_5percent]
#     top_outlier[i,] = names(a)
#   }
#   
#   a=rbind(top_outlier, top_outlier_zscore)
#   a = a[order(row.names(a)),]
#   fn = paste(pd,name,'outlier_m-zscore.txt', sep="_")
#   write.table(a, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
#   
# 
#   b=rbind(top_outlier, top_outlier_raw)
#   b = b[order(row.names(b)),]
#   fn = paste(pd,name,'outlier_raw_exp.txt', sep="_")
#   write.table(b, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
#   
#   if (plot){
#     # plot
#     top_outlier.m <- melt(top_outlier[,c(1:5)])
#     top_outlier_zscore.m <- melt(top_outlier_zscore[,c(1:5)])
#     top_outlier_boolean.m <- melt(top_outlier_boolean[,c(1:5)])
#     
#     fn = paste(pd, name, 'top5_mzscore.pdf',sep ="_")
#     YlGnBu = brewer.pal(9, "YlGnBu") 
#     getPalette = colorRampPalette(YlGnBu)
#     outlier.colors=c("NA", "#000000")
#     
#     p = ggplot()
#     p = p + geom_tile(data=top_outlier_zscore.m, aes(x=as.factor(X1), y=X2, fill=value), linetype="blank") + scale_fill_gradientn(name= "modified z-score", colours=getPalette(100))
#     p = p + geom_tile(data=top_outlier_boolean.m, aes(x=as.factor(X1), y=X2, color=value), fill=NA, size=0.5) + scale_colour_manual(name="KDE outlier",values = outlier.colors)
#     p = p + geom_text(data=top_outlier.m,aes(x=as.factor(X1), y=X2, label = value), color="red", size=4, angle=90)
#     p = p + xlab("sample") + ylab("top outlier targets") + theme_bw() + 
#       theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=16), axis.text.y = element_text(colour="black", size=16))
#     p
#     ggsave(file=fn, height=h, width=w)
#   }
#   # return the top outliers
#   return(list("outlier_zscore"=outlier_zscore, "outlier"=outlier,
#     "top_outlier_zscore"=top_outlier_zscore, "top_outlier"=top_outlier, "top_outlier_boolean"=top_outlier_boolean))
# }
grubbs.flag <- function(x) {
  outliers <- NULL
  test <- x
  grubbs.result <- grubbs.test(test)
  pv <- grubbs.result$p.value
  while(pv < 0.05) {
    outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]))
    test <- x[!x %in% outliers]
    grubbs.result <- grubbs.test(test)
    pv <- grubbs.result$p.value
  }
  return(data.frame(X=x,Outlier=(x %in% outliers)))
}



##### modified 20150828: use the z-score distribution of all samples to find p value, and use significance level to define outlier#####
# choose only the genes with SD with at least 1 to test
# try both conventional z-score and modified z-score
#find_outlier_p = function(m, name, genes=FALSE, plot=TRUE, h=6, w=10, minNum = 10){ 
#   cat("##### OUTLIER ANALYSIS #####\n")
#   m = as.matrix(m)
#   num = nrow(m)
#   m2 = as.matrix(m[rowSums(!is.na(m)) >= minNum, ])
#   num_NA= nrow(m2)
#   if (genes[1] != FALSE){
#     m2_genes = m2[row.names(m2) %in% genes, ]
#   } else{ m2 = m; m2_genes=m}
#   num_gene = nrow(m2_genes)
#   m2_genes = as.matrix(m2_genes)
#   cat(paste("Looking for outliers in", deparse(substitute(genes)), "of", deparse(substitute(m)), "\n", sep=" "))
#   cat(paste("Original number of genes:", num, "; NA filtered:", num_NA, "; Gene list filtered:", num_gene, "\n", sep=" "))
#   
#   dis_dir_cmd = paste("mkdir figures/", date, "/Dis/", sep="")
#   system(dis_dir_cmd)
#   outlier = matrix(,nrow=dim(m2)[1],ncol=dim(m2)[2])
#   row.names(outlier) = row.names(m2)
#   colnames(outlier) = colnames(m2)
#   outlier_mzscore = outlier
#   outlier_zscore = outlier
#   outlier_pvalue = outlier
#   sdrow=c()
#   maxrow=c()
#   meanrow=c()
#   
#   # gene-wise zscore
#   for (i in 1:nrow(m2)){
#     # modified z-score for outlier: Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
#     outlier_mzscore[i,]  = 0.6745*(m2[i,]-median(m2[i,], na.rm=TRUE))/mad(m2[i,], na.rm=TRUE)
# #     # conventional z-score
# #     outlier_zscore[i,] = (m2[i,] - mean(m2[i,], na.rm=TRUE))/sd(m2[i,], na.rm=TRUE) 
#     sdrow=c(sdrow,sd(m2[i,], na.rm=TRUE))
#     maxrow=c(maxrow,max(m2[i,], na.rm=TRUE))
#     meanrow=c(meanrow,mean(m2[i,], na.rm=TRUE))
#   }
#   #qualified_genes=row.names(outlier)[sdrow > 1.5]
#   qualified_genes_max=row.names(outlier)[maxrow > 2.5]
#   #qualified_genes_mean=row.names(outlier)[meanrow > 0.5]
#   qualified_intersect = qualified_genes[qualified_genes_mean %in% qualified_genes_max]
#   
#   # cohort-wise outlier: define the null Empirical Cumulative Distribution Function based on the whole cohort
#   outlier_percentile = ecdf(outlier_mzscore)
#   for (i in 1:nrow(outlier_mzscore)){
#     outlier_pvalue[i,] = 1-outlier_percentile(outlier_mzscore[i,])
#   }
#   
#   # test only the qualified genes and calculate FDR
#   q_genes = genes[genes %in% qualified_genes_max]
# 
#   outlier_pvalue_genes = outlier_pvalue[q_genes,]
#   
#   # FDR calc
#   outlier_FDR_genes=p.adjust(outlier_pvalue_genes, method="BH")
#   outlier_FDR_genes=matrix(outlier_FDR_genes, nrow=length(q_genes))
#   row.names(outlier_FDR_genes)=row.names(outlier_pvalue_genes)
#   colnames(outlier_FDR_genes)=colnames(outlier_pvalue_genes)
# 
# ## plot ## 
# if (plot){ # make the genes of interest different color, and mark significance line
#   df = as.data.frame(outlier_zscore_s)
#   fn = paste("figures/", date, "/Dis/",date, "_", sample, "_", name, '_mzscore_distribution.pdf', sep = "")
#   p = ggplot(data=df, aes_string(x=sample)) + geom_density()
#   p = p + geom_vline(xintercept=maxs, colour="red", alpha=0.2) + geom_vline(xintercept=mins, colour="blue", alpha=0.2) + geom_vline(xintercept=adj_min, colour="blue")
#   p = p + xlab(paste("modified Z-score in", name, sep=" ")) + 
#     scale_x_continuous(limits=c(-5,5),breaks=c(mins,maxs)) + theme_bw() + 
#     theme(text = element_text(colour="black", size=16), axis.text.x = element_text(angle = 90,colour="black", size=14), axis.text.y = element_text(colour="black", size=14))
#   p
#   ggsave(file=fn)
# }
# outlier[,i] = outlier_pvalue[,i] <= adj_min
# print(paste(sample, name, "FDR0.05_outlier_threshold:", adj_min, "Number_of_outliers", sum(outlier[,i], na.rm=T), sep = "  "))
# 
# 
#   # set up return matrixes
#   zscore=outlier_zscore
#   outlier=outlier
#   num_5percent = dim(zscore)[1]
#   top_outlier_zscore = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
#   top_outlier = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
#   top_outlier_boolean = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
#   top_outlier_raw = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
#   row.names(top_outlier_zscore)=colnames(zscore)
#   row.names(top_outlier)=colnames(zscore)
#   row.names(top_outlier_boolean)=colnames(zscore)
#   row.names(top_outlier_raw)=colnames(zscore)
#   colnames(top_outlier_zscore)=c(1:num_5percent)
#   for (i in 1:num_5percent){colnames(top_outlier_zscore)[i] = paste(name, "top", colnames(top_outlier_zscore)[i], sep=" ")}
#   colnames(top_outlier)=colnames(top_outlier_zscore)
#   colnames(top_outlier_boolean)=colnames(top_outlier_zscore)
#   colnames(top_outlier_raw)= colnames(top_outlier_zscore)
#   # find top 5
#   for (i in colnames(zscore)){
#     whim=zscore[,i]
#     a = whim[order(whim, decreasing=TRUE)][1:num_5percent]
#     top_outlier_zscore[i,] = a
#     whim2 = outlier[,i]
#     top_outlier_boolean[i,] = whim2[order(whim, decreasing=TRUE)][1:num_5percent]
#     whim3 = m2_genes[,i]
#     top_outlier_raw[i,] = whim3[order(whim, decreasing=TRUE)][1:num_5percent]
#     top_outlier[i,] = names(a)
#   }
#   
#   a=rbind(top_outlier, top_outlier_zscore)
#   a = a[order(row.names(a)),]
#   fn = paste(pd,name,'outlier_m-zscore.txt', sep="_")
#   write.table(a, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
#   
#   
#   b=rbind(top_outlier, top_outlier_raw)
#   b = b[order(row.names(b)),]
#   fn = paste(pd,name,'outlier_raw_exp.txt', sep="_")
#   write.table(b, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
#   
#   if (plot){
#     # plot
#     top_outlier.m <- melt(top_outlier[,c(1:5)])
#     top_outlier_zscore.m <- melt(top_outlier_zscore[,c(1:5)])
#     top_outlier_boolean.m <- melt(top_outlier_boolean[,c(1:5)])
#     
#     fn = paste(pd, name, 'top5_mzscore.pdf',sep ="_")
#     YlGnBu = brewer.pal(9, "YlGnBu") 
#     getPalette = colorRampPalette(YlGnBu)
#     outlier.colors=c("NA", "#000000")
#     
#     p = ggplot()
#     p = p + geom_tile(data=top_outlier_zscore.m, aes(x=as.factor(X1), y=X2, fill=value), linetype="blank") + scale_fill_gradientn(name= "modified z-score", colours=getPalette(100))
#     p = p + geom_tile(data=top_outlier_boolean.m, aes(x=as.factor(X1), y=X2, color=value), fill=NA, size=0.5) + scale_colour_manual(name="KDE outlier",values = outlier.colors)
#     p = p + geom_text(data=top_outlier.m,aes(x=as.factor(X1), y=X2, label = value), color="red", size=4, angle=90)
#     p = p + xlab("sample") + ylab("top outlier targets") + theme_bw() + 
#       theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=16), axis.text.y = element_text(colour="black", size=16))
#     p
#     ggsave(file=fn, height=h, width=w)
#   }
#   # return the top outliers
#   return(list("outlier_zscore"=outlier_zscore, "outlier"=outlier,
#               "top_outlier_zscore"=top_outlier_zscore, "top_outlier"=top_outlier, "top_outlier_boolean"=top_outlier_boolean))
# }
plot=TRUE; h=6; w=10; minNum = 10
##### use the box plot definition of outlier, then rank them by the outlier score ##### 
find_outlier = function(m, name="dataset", plot=TRUE, h=6, w=10, minNum = 10, whim_only=F){ 
  #w=40 for human panels with ~80 samples
  cat("##### OUTLIER ANALYSIS #####\n")
  m = as.matrix(m)
  num = nrow(m)
  m2 = as.matrix(m[rowSums(!is.na(m)) >= minNum, ])
  num_NA= nrow(m2)
  cat(paste("Looking for outliers in", deparse(substitute(genes)), "of", name, "\n", sep=" "))
  cat(paste("Original number of genes:", num, "; NA filtered:", num_NA, "\n", sep=" "))
  
  #   dis_dir_cmd = paste("mkdir figures/", date, "/Dis/", sep="")
  #   system(dis_dir_cmd)
  outlier = matrix(,nrow=dim(m2)[1],ncol=dim(m2)[2])
  row.names(outlier) = row.names(m2)
  colnames(outlier) = colnames(m2)
  outlier_mzscore = outlier
  outlier_box = outlier
  #outlier_box2 = outlier # more stringent outlier definition based on outer fences
  
  # gene-wise outlier and outlier score
  for (i in 1:nrow(m2)){
    # modified z-score for outlier: Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    # outlier_mzscore[i,]  = 0.6745*(m2[i,]-median(m2[i,], na.rm=TRUE))/mad(m2[i,], na.rm=TRUE)
    # box-plot definition of outlier
    IQR = quantile(m2[i,], probs=0.75, na.rm=T) - quantile(m2[i,], probs=0.25, na.rm=T) 
    outlier_box[i,] = (m2[i,] >= quantile(m2[i,], probs=0.75, na.rm=T) + 1.5*IQR)
    outlier_mzscore[i,] = (m2[i,] - quantile(m2[i,], probs=0.75, na.rm=T))/(1.5*IQR) #inner fences
    # outlier_box2[i,] = (m2[i,] >= quantile(m2[i,], probs=0.75, na.rm=T) + 3.5*IQR) #outer fences
  }
  # output the outlier score table
  fn = paste(pd,name,'outlier_score_table.txt', sep="_")
  write.table(outlier_mzscore, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
  
  # extract gene list of interest
  #   outlier_box_d = outlier_box[row.names(outlier_box) %in% genes,]
  #   outlier_box_d[rowSums(outlier_box_d,na.rm=T)>=1,]
  #   #   outlier_box2_d = outlier_box2[row.names(outlier_box2) %in% genes,]
  #   #   outlier_box2_d[rowSums(outlier_box2_d,na.rm=T)>=1,]
  #   outlier_mzscore_d = outlier_mzscore[row.names(outlier_mzscore) %in% genes,]
  
  num_outliers = sum(outlier_box, na.rm=T)
  cat(paste("Number_of_samples:", dim(outlier)[2], "Number_of_outliers:", num_outliers,"; Avg_outlier_per_sample:", num_outliers/dim(outlier)[2], "\n\n", sep = " "))
  
  # set up return matrixes
  zscore=outlier_mzscore
  outlier=outlier_box
  num_genes = dim(zscore)[1] 
  top_outlier_zscore = matrix(,nrow=dim(zscore)[2],ncol=num_genes)
  top_outlier = matrix(,nrow=dim(zscore)[2],ncol=num_genes)
  top_outlier_boolean = matrix(,nrow=dim(zscore)[2],ncol=num_genes)
  top_outlier_raw = matrix(,nrow=dim(zscore)[2],ncol=num_genes)
  row.names(top_outlier_zscore)=colnames(zscore)
  row.names(top_outlier)=colnames(zscore)
  row.names(top_outlier_boolean)=colnames(zscore)
  row.names(top_outlier_raw)=colnames(zscore)
  colnames(top_outlier_zscore)=c(1:num_genes)
  #for (i in 1:num_genes){colnames(top_outlier_zscore)[i] = paste(name, "top", colnames(top_outlier_zscore)[i], sep=" ")}
  colnames(top_outlier)=colnames(top_outlier_zscore)
  colnames(top_outlier_boolean)=colnames(top_outlier_zscore)
  colnames(top_outlier_raw)= colnames(top_outlier_zscore)
  # rank order based on zscore
  for (i in colnames(zscore)){
    whim=zscore[,i]
    a = whim[order(whim, decreasing=TRUE)][1:num_genes]
    top_outlier_zscore[i,] = a
    whim2 = outlier[,i]
    top_outlier_boolean[i,] = whim2[order(whim, decreasing=TRUE)][1:num_genes]
    whim3 = m2[,i]
    top_outlier_raw[i,] = whim3[order(whim, decreasing=TRUE)][1:num_genes]
    top_outlier[i,] = names(a)
  }
  
  a=rbind(top_outlier, top_outlier_zscore)
  a = a[order(row.names(a)),]
  fn = paste(pd,name,'outlier_score.txt', sep="_")
  write.table(a, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
  
  
  b=rbind(top_outlier, top_outlier_raw)
  b = b[order(row.names(b)),]
  fn = paste(pd,name,'outlier_raw_exp.txt', sep="_")
  write.table(b, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
  
  c=rbind(top_outlier, top_outlier_boolean)
  c = c[order(row.names(c)),]
  fn = paste(pd,name,'outlier.txt', sep="_")
  write.table(c, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
  # how many outliers should be shown?
  num_shown=1
  for (i in 1:nrow(top_outlier_boolean)){
    row_outlier = sum(top_outlier_boolean[i,], na.rm=T)
    if (row_outlier > num_shown) {num_shown = row_outlier}
  }
  
  # plot only the outliers: has to do through command line because of not specifying "data" in ggplot2
  if (FALSE){ 
    #num_shown = 3
    if (num_shown > 5) {num_shown=5}
    top_outlier2 = top_outlier[grepl("WHIM",row.names(top_outlier)),]
    top_outlier_zscore2 = top_outlier_zscore[grepl("WHIM",row.names(top_outlier_zscore)),]
    top_outlier_boolean2 = top_outlier_boolean[grepl("WHIM",row.names(top_outlier_boolean)),]
    # plot
    top_outlier.m <- melt(top_outlier2[,c(1:num_shown)])
    top_outlier_zscore.m <- melt(top_outlier_zscore2[,c(1:num_shown)])
    top_outlier_boolean.m <- melt(top_outlier_boolean2[,c(1:num_shown)])
    
    fn = paste(pd, name, 'top5_outlier_only_score.pdf',sep ="_")
    YlOrRd = brewer.pal(9, "YlOrRd") 
    getPalette = colorRampPalette(YlOrRd)
    outlier.colors=c("NA", "#000000")
    top_outlier.m$value = as.character(top_outlier.m$value)
    p = ggplot()
    p = p + geom_tile(aes(x=as.factor(top_outlier_zscore.m$X1), y=top_outlier_zscore.m$X2, fill=ifelse(top_outlier_boolean.m$value,top_outlier_zscore.m$value,NA)), linetype="blank") + scale_fill_gradientn(name= "Outlier score", colours=getPalette(100), na.value=NA, limits=c(0,6))
    p = p + geom_tile(data=top_outlier_boolean.m, aes(x=as.factor(X1), y=X2, color=value), fill=NA, size=0.5) + scale_colour_manual(name="Outlier",values = outlier.colors)
    p = p + geom_text(aes(x=top_outlier.m$X1, y=top_outlier.m$X2, label = ifelse(top_outlier_boolean.m$value,top_outlier.m$value,NA), stringsAsFactors=FALSE), color="black", size=5, angle=90)
    p = p + xlab("Sample") + ylab(paste("Top", name,"outliers", sep=" ")) + theme_bw() + 
      theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=14), axis.text.y = element_blank(),axis.ticks.y = element_blank())#element_text(colour="black", size=14))
    p
    ggsave(file=fn, height=8, width=w)   
  }  
  
  # version that plotted everything
  if (plot){
    # plot
    #num_shown = 4
    top_outlier.m <- melt(top_outlier[,c(1:num_shown)])
    top_outlier_zscore.m <- melt(top_outlier_zscore[,c(1:num_shown)])
    top_outlier_boolean.m <- melt(top_outlier_boolean[,c(1:num_shown)])
    
    fn = paste(pd, name, 'top5_outlier_score_all.pdf',sep ="_")
    #     YlGnBu = brewer.pal(9, "YlGnBu") 
    #     getPalette = colorRampPalette(YlGnBu)
    YlOrRd = brewer.pal(9, "YlOrRd") 
    getPalette = colorRampPalette(YlOrRd)
    outlier.colors=c("NA", "#000000")
    
    p = ggplot()
    p = p + geom_tile(data=top_outlier_zscore.m, aes(x=as.factor(X1), y=X2, fill=value), linetype="blank") + scale_fill_gradientn(name= "Outlier score", colours=getPalette(100))
    p = p + geom_tile(data=top_outlier_boolean.m, aes(x=as.factor(X1), y=X2, color=value), fill=NA, size=0.5) + scale_colour_manual(name="Outlier",values = outlier.colors)
    p = p + geom_text(data=top_outlier.m,aes(x=as.factor(X1), y=X2, label = value), color="black", size=4, angle=90)
    p = p + xlab("Sample") + ylab("Top Druggable Outliers") + theme_bw() + 
      theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=16), axis.text.y = element_blank(),axis.ticks.y = element_blank())#element_text(colour="black", size=16))
    p
    ggsave(file=fn, height=h, width=w)
  }
  
  # plot only the WHIMs
  if (whim_only){
    # get only WHIMs
    top_outlier2 = top_outlier[grepl("WHIM",row.names(top_outlier)),]
    top_outlier_zscore2 = top_outlier_zscore[grepl("WHIM",row.names(top_outlier_zscore)),]
    top_outlier_boolean2 = top_outlier_boolean[grepl("WHIM",row.names(top_outlier_boolean)),]
    # plot
    top_outlier.m <- melt(top_outlier2[,c(1:5)])
    top_outlier_zscore.m <- melt(top_outlier_zscore2[,c(1:5)])
    top_outlier_boolean.m <- melt(top_outlier_boolean2[,c(1:5)])
    
    fn = paste(pd, name, 'top5_outlier_score.pdf',sep ="_")
    YlOrRd = brewer.pal(9, "YlOrRd") 
    getPalette = colorRampPalette(YlOrRd)
    outlier.colors=c("NA", "#000000")
    
    p = ggplot()
    p = p + geom_tile(data=top_outlier_zscore.m, aes(x=as.factor(X1), y=X2, fill=value), linetype="blank") + scale_fill_gradientn(name= "Outlier score", colours=getPalette(100))
    p = p + geom_tile(data=top_outlier_boolean.m, aes(x=as.factor(X1), y=X2, color=value), fill=NA, size=0.5) + scale_colour_manual(name="Outlier",values = outlier.colors)
    p = p + geom_text(data=top_outlier.m,aes(x=as.factor(X1), y=X2, label = value), color="black", size=4, angle=90)
    p = p + xlab("Sample") + ylab("Top Druggable Outliers") + theme_bw() + 
      theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=14), axis.text.y = element_blank())#element_text(colour="black", size=14))
    p
    ggsave(file=fn, height=h, width=w)    
  }
  
  # return the top outliers
  return(list("outlier_zscore"=outlier_mzscore, "outlier"=outlier,
              "top_outlier_zscore"=top_outlier_zscore, "top_outlier"=top_outlier, "top_outlier_boolean"=top_outlier_boolean))
}

##### OUTLIER IN ITRAQ PROTEOME #####
ITRAQ = read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')
ITRAQ.d = ITRAQ[row.names(ITRAQ) %in% druggable,]
ITRAQ_druggable = find_outlier(ITRAQ.d, name = "ITRAQ druggable proteome")

##### OUTLIER IN LFQ#####
LFQ=read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned_collapsed.txt')
LFQ.d = LFQ[row.names(LFQ) %in% druggable,]
LFQ_druggable = find_outlier(LFQ.d, "LFQ druggable proteome")

### overlap between ITRAQ and LFQ
LFQ_top_outlier2=LFQ_druggable$top_outlier[,1:5]
LFQ_top_outlier_zscore2=LFQ_druggable$top_outlier_zscore[,1:5]
LFQ_top_outlier_boolean2=LFQ_druggable$top_outlier_boolean[,1:5]
top_outlier2=ITRAQ_druggable$top_outlier[,1:5]
top_outlier_zscore2=ITRAQ_druggable$top_outlier_zscore[,1:5]
top_outlier_boolean2=ITRAQ_druggable$top_outlier_boolean[,1:5]
top_overlap = merge(LFQ_top_outlier2, top_outlier2, by = "row.names", all.x=T)
colnames(top_overlap)[1] = "WHIM"
top_zscore_overlap = merge(LFQ_top_outlier_zscore2, top_outlier_zscore2, by = "row.names", all.x=T)
colnames(top_zscore_overlap)[1] = "WHIM"
top_boolean_overlap = merge(LFQ_top_outlier_boolean2, top_outlier_boolean2, by = "row.names", all.x=T)
colnames(top_boolean_overlap)[1] = "WHIM"

# plot the overlap 
top_overlap.m <- melt(top_overlap, id.var = "WHIM")
top_zscore_overlap.m <- melt(top_zscore_overlap, id.var = "WHIM")
top_boolean_overlap.m <- melt(top_boolean_overlap, id.var = "WHIM")
colnames(top_boolean_overlap.m)[3]="outlier"
colnames(top_zscore_overlap.m)[3]="outlier_score"
top_overlap.m$variable = sub("druggable proteome ","",top_overlap.m$variable)
top_zscore_overlap.m$variable = sub("druggable proteome ","",top_zscore_overlap.m$variable)
top_boolean_overlap.m$variable = sub("druggable proteome ","",top_boolean_overlap.m$variable)

pdf(paste(pd,'ITRAQ_LFQ_top5_outlier_score_proteome.pdf',sep ="_"), height=6, width=9)
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette = colorRampPalette(YlGnBu)
# YlOrRd = brewer.pal(9, "YlOrRd") 
# getPalette = colorRampPalette(YlOrRd)
outlier.colors=c("NA", "#000000")

p = ggplot()
p = p + geom_tile(data=top_zscore_overlap.m, aes(y=as.factor(WHIM), x=variable, fill=outlier_score), linetype="blank") + scale_fill_gradientn(colours=getPalette(100))
p = p + geom_tile(data=top_boolean_overlap.m, aes(y=as.factor(WHIM), x=variable, color=outlier), fill=NA, size=0.5) + scale_colour_manual(values = outlier.colors)
p = p + geom_text(data=top_overlap.m,aes(y=as.factor(WHIM), x=variable, label = value), color="red", size=3)
p = p + ylab("Sample") + xlab("Top druggable protein") + theme_bw() + 
  theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=12), axis.text.y = element_text(,colour="black", size=12))
p
dev.off()

##### OUTLIER IN ITRAQ PHOSPHO ##### 
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
# 56651 phosphosites
ITRAQpho=ITRAQpho[,-c(1,2)]
cat("Original number of ITRAQ phosphosites: 56651\n")
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQpho.na = ITRAQpho[,-c(17,18,20)]
genes = sub("-NP.*", "", row.names(ITRAQpho.na))
row.names(ITRAQpho.na) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("\\. _.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("\\._.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub(" _.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("_.*","",row.names(ITRAQpho.na)), unique=T)
ITRAQpho.na.d = ITRAQpho.na[genes %in% druggable, ] # 1167 phosphosites
cat("Druggable list filtered ITRAQ phosphosites: 1167\n")

# outlier analysis
ITRAQpho_druggable = find_outlier(ITRAQpho.na.d, "ITRAQ druggable phosphoproteome", h=12)


##### OUTLIER IN LFQ phospho #####
LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)
LFQpho$sites=sub(".*\\(","",LFQpho$phospho_site)
LFQpho$sites=sub("\\)","",LFQpho$site)
LFQpho$sites=paste(LFQpho$gene_name,LFQpho$sites, sep=".")
row.names(LFQpho) = make.names(LFQpho$sites, unique=T)
colnames(LFQpho) = sub("_P", "", colnames(LFQpho))
colnames(LFQpho) = sub("Intensity.W", "WHIM", colnames(LFQpho))
cat("Original number of LFQ phosphosites: 18229\n")
# 18229 phosphosites
LFQpho.na = LFQpho[,-c(19:25)]
genes = sub("\\..*", "", row.names(LFQpho.na))
LFQpho.na.d = LFQpho.na[genes %in% druggable, ] #331 phosphosites
cat("Druggable list filtered LFQ phosphosites: 331\n")

LFQpho_druggable = find_outlier(LFQpho.na.d, "LFQ druggable phosphoproteome", h=12)

### find extreme phosphos that are not indicated by pro 
# phospho (ITRAQpho_outlier_zscore) and pro (outlier_zscore) in iTRAQ of the same sample, see how they correlate in z-score
ITRAQpho_outlier_zscore.dt = as.data.frame(ITRAQpho_druggable$outlier_zscore)
ITRAQpho_outlier_zscore.dt$Gene = sub("\\..*", "", row.names(ITRAQpho_outlier_zscore.dt))
ITRAQpho_outlier_zscore.dt$Site = row.names(ITRAQpho_outlier_zscore.dt)
outlier_zscore.dt = as.data.frame(ITRAQ_druggable$outlier_zscore)
outlier_zscore.dt$Gene = sub("\\..*", "", row.names(outlier_zscore.dt))

ITRAQpho_outlier_zscore.m = melt(ITRAQpho_outlier_zscore.dt, id.var = c("Gene","Site"))
outlier_zscore.dt.m = melt(outlier_zscore.dt, id.var = "Gene")
pro_pho_overlap = merge(ITRAQpho_outlier_zscore.m, outlier_zscore.dt.m, by = c("Gene","variable"), all.x=T)

corPP = cor(pro_pho_overlap$value.x, pro_pho_overlap$value.y, use = 'pairwise.complete.obs', method = "pearson")
cat(paste("\n","##### Pearson correlation for protein and mRNA outlier score:",corPP, "#####\n", sep=" "))

#pro_pho_overlap$gene = rep("Other", nrow(pro_pho_overlap))
#pro_pho_overlap[pro_pho_overlap$Gene == "ERBB2",]$gene="ERBB2"
#pro_pho_overlap[pro_pho_overlap$Gene == "PAK1",]$gene="PAK1"

pdf(paste(pd,'ITRAQ_phospho_vs_pro_score.pdf', sep="_"),height=7, width =7)

p = ggplot(data = pro_pho_overlap, aes(x=value.y, y=value.x, colour=variable))#, shape=gene)) 
p = p + geom_point(alpha=0.3, size=1.5) + xlab("Protein expression outlier score") + ylab("Phosphosite expression outlier score") 
p = p + scale_shape_manual(values=c(0,16,2)) + guides(colour=FALSE)
p = p + theme_bw() + coord_fixed() #+ theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))
p = p + xlim(-5,5) + ylim(-5,5) + geom_vline(xintercept=1.5,alpha=0.2) + geom_hline(yintercept=1.5,alpha=0.2)
#p = p + geom_text(aes(label=ifelse(value.x > 3 & value.x > value.y + 2,Site,'')),hjust=-0.05,just=0,size=3,alpha=0.7)
#p = p + scale_x_continuous(expand=c(0.02,0)) + scale_y_continuous(expand=c(0.02,0))
p = p + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
#p
#dev.off()

p2 = ggplot(data = pro_pho_overlap,aes(x=value.y, colour = variable, fill=variable)) + 
  geom_density(alpha=0.1) + guides(colour=FALSE) + #geom_vline(xintercept=1.5, alpha=0.2)  + 
  xlim(-5,5) + 
  theme_bw() +
  theme0(plot.margin = unit(c(1,0,-1.52,3.5),"lines"))
#theme0(plot.margin = unit(c(1,0,-0.48,2.2),"lines"))

p3 = ggplot(data = pro_pho_overlap,aes(x=value.x, colour = variable, fill=variable)) + 
  geom_density(alpha=0.1) + guides(colour=FALSE) + #geom_vline(xintercept=1.5, alpha=0.2) + 
  xlim(-5,5) + 
  theme_bw() + 
  coord_flip()  +
  theme0(plot.margin = unit(c(0.7,0,2.8,-0.85),"lines"))
#theme0(plot.margin = unit(c(0,1,1.2,-0.48),"lines"))

grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3))

dev.off()
#http://stackoverflow.com/questions/17370460/scatterplot-with-alpha-transparent-histograms-in-r

##### RNA-Seq #####
RNA = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")
row.names(RNA)=make.names(RNA$gene, unique=T)
RNA=RNA[,-1] #16209 genes
RNA = RNA[,colnames(RNA) %in% colnames(ITRAQ)]
RNA.d = RNA[row.names(RNA) %in% druggable,]

RNA_druggable = find_outlier(RNA.d, "RNA druggable genes")#, w=12)

##### find outlier proteome expression not identified by transcriptome #####
outlier_zscore.dt = as.data.frame(unlist(ITRAQ_druggable$outlier_zscore))
outlier_zscore.dt$Gene = row.names(outlier_zscore.dt)
RNA_outlier_zscore.dt = as.data.frame(unlist(RNA_druggable$outlier_zscore))
RNA_outlier_zscore.dt$Gene = row.names(RNA_outlier_zscore.dt)

outlier_zscore.dt.m = melt(outlier_zscore.dt, id.var = "Gene")
RNA_outlier_zscore.dt.m = melt(RNA_outlier_zscore.dt, id.var = "Gene")
pro_rna_overlap = merge(RNA_outlier_zscore.dt.m, outlier_zscore.dt.m, by = c("Gene","variable"))

corPR=cor(pro_rna_overlap$value.x, pro_rna_overlap$value.y, use = 'pairwise.complete.obs', method="pearson")
cat(paste("\n","##### Pearson correlation for protein and mRNA outlier score:",corPR, "#####\n", sep=" "))

#pro_rna_overlap$gene = rep("Other", nrow(pro_rna_overlap))
#pro_rna_overlap[pro_rna_overlap$Gene == "ERBB2",]$gene="ERBB2"
#pro_rna_overlap[pro_rna_overlap$Gene == "AKT2",]$gene="AKT2"
#pro_rna_overlap[pro_rna_overlap$Gene == "AURKA",]$gene="AURKA" #found to be amplified in basal breast tumor

pdf(paste(pd,'rna_vs_pro_outlier_score.pdf', sep="_"),height=7, width =7)
p = ggplot(data = pro_rna_overlap, aes(x=value.x, y=value.y, colour=variable))#, shape=gene)) 
p = p + geom_point(alpha=0.3, size=1.5) + xlab("mRNA expression outlier score") + ylab("Protein expression outlier score")
p = p + scale_shape_manual(values=c(0,2,16)) + guides(colour=FALSE)
p = p + theme_bw() + coord_fixed() #+ theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))
p = p + xlim(-5,5) + ylim(-5,5) + geom_vline(xintercept=1.5,alpha=0.2) + geom_hline(yintercept=1.5,alpha=0.2)
p = p + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))

p2 = ggplot(data = pro_rna_overlap,aes(x=value.x, colour = variable, fill=variable)) + 
  geom_density(alpha=0.1) + guides(colour=FALSE) + #geom_vline(xintercept=1.5, alpha=0.2) + 
  xlim(-5,5) + 
  theme_bw() +
  theme0(plot.margin = unit(c(1,0,-1.52,3.5),"lines"))
  #theme0(plot.margin = unit(c(0,0,0,0),"lines"))

p3 = ggplot(data = pro_rna_overlap,aes(x=value.y, colour = variable, fill=variable)) + 
  geom_density(alpha=0.1) + guides(colour=FALSE) + #geom_vline(xintercept=1.5, alpha=0.2) +  
  xlim(-5,5) + 
  theme_bw() + 
  coord_flip()  +
  theme0(plot.margin = unit(c(0.7,0,2.8,-0.85),"lines"))
  #theme0(plot.margin = unit(c(0,0,0,0),"lines"))

grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3))
dev.off()

# remove_axis =
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank())
# p = p + remove_axis
# p2 = p2 + remove_axis
# p3= p3 + remove_axis
# p4 = p4 + remove_axis
# 
# xyplot.gb = ggplot_build(p)
# densityY.gb = ggplot_build(p2)
# densityX.gb = ggplot_build(p3)
# 
# # densityY.gb$panel$ranges[[1]]$y.range = xyplot.gb$panel$ranges[[1]]$y.range   
# # densityX.gb$panel$ranges[[1]]$x.range = xyplot.gb$panel$ranges[[1]]$x.range
# 
# xyplot.gt = ggplot_gtable(xyplot.gb)
# densityY.gt = ggplot_gtable(densityY.gb)
# densityX.gt = ggplot_gtable(densityX.gb)
# 
# #xyplot.gt_panel = subset(xyplot.gt$layout, name == "panel")
# densityY.gt$heights = xyplot.gt$heights
# densityX.gt$widths = xyplot.gt$widths
# 
# main.grob = arrangeGrob(densityY.gt, p4, xyplot.gt, densityX.gt, ncol=2, nrow=2) #widths=c(0.3,0.7), heights=c(0.7,0.3), 
# grid.draw(main.grob)
# #pdf(file=paste(pd,'rna_vs_pro_outlier_score.pdf', sep="_"), useDingbats=FALSE)
# dev.off()
# 
# pdf(paste(pd,'rna_vs_pro_outlier_score2.pdf', sep="_"))
# grid.arrange(arrangeGrob(densityY.gt,ncol=2,widths=c(3,1)),
#              arrangeGrob(xyplot.gt,densityX.gt,ncol=2,widths=c(3,1)),
#              heights=c(1,3))
##### CNV #####
CNV = read.table(row.names=1,header=TRUE, sep="\t", file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified_normalized.tsv")
CNV.d = CNV[row.names(CNV) %in% druggable,]
# find outliers
CNV_druggable = find_outlier(CNV.d, "CNV druggable genes", w=12)

### all levels: CNV, RNA, proteome, phosphoproteome ### 
# keep the 24 WHIMs in ITRAQ proteome
# all_levels=list("LFQ_druggable"=LFQ_druggable,"RNA_druggable"=RNA_druggable,"CNV_druggable"=CNV_druggable)
# top_overlap_all = merge(ITRAQpho_druggable$top_outlier[,c(1:5)], ITRAQ_druggable$top_outlier[,c(1:5)], by = "row.names", all.x=T)
# for (i in all_levels){
#   rownames(top_overlap_all)=top_overlap_all[,1]
#   top_overlap_all=top_overlap_all[,-1]
#   top_overlap_all = merge(top_overlap_all, i, by = "row.names", all.x=T)
# }
# colnames(top_overlap_all)[1] = "sample"
# 
# all_levels2=list("top_outlier_zscore"=top_outlier_zscore,"LFQ_top_outlier_zscore"=LFQ_top_outlier_zscore,"RNA_top_outlier_zscore"=RNA_top_outlier_zscore,"CNV_top_outlier_zscore"=CNV_top_outlier_zscore)
# ITRAQpho_top_outlier_zscore.n = ITRAQpho_top_outlier_zscore/mean(ITRAQpho_top_outlier_zscore)
# LFQpho_top_outlier_zscore.n = LFQpho_top_outlier_zscore/mean(LFQpho_top_outlier_zscore)
# top_overlap_all_zscore = merge(ITRAQpho_top_outlier_zscore.n, LFQpho_top_outlier_zscore.n, by = "row.names", all.x=T)
# for (i in all_levels2){
#   rownames(top_overlap_all_zscore)=top_overlap_all_zscore[,1]
#   top_overlap_all_zscore=top_overlap_all_zscore[,-1]
#   i.n=i/mean(i)
#   top_overlap_all_zscore = merge(top_overlap_all_zscore, i.n, by = "row.names", all.x=T)
# }
# colnames(top_overlap_all_zscore)[1] = "sample"
# 
# all_levels3=list("top_outlier_boolean"=top_outlier_boolean,"LFQ_top_outlier_boolean"=LFQ_top_outlier_boolean,"RNA_top_outlier_boolean"=RNA_top_outlier_boolean,"CNV_top_outlier_boolean"=CNV_top_outlier_boolean)
# top_overlap_all_boolean = merge(ITRAQpho_top_outlier_boolean, LFQpho_top_outlier_boolean, by = "row.names", all.x=T)
# for (i in all_levels3){
#   rownames(top_overlap_all_boolean)=top_overlap_all_boolean[,1]
#   top_overlap_all_boolean=top_overlap_all_boolean[,-1]
#   top_overlap_all_boolean = merge(top_overlap_all_boolean, i, by = "row.names", all.x=T)
# }
# colnames(top_overlap_all_boolean)[1] = "sample"
# 
# top_overlap_all.m <- melt(top_overlap_all, id.var = "sample")
# top_overlap_all_zscore.m <- melt(top_overlap_all_zscore, id.var = "sample")
# top_overlap_all_boolean.m <- melt(top_overlap_all_boolean, id.var = "sample")
# 
# fn=paste(pd,'all_level_top5_mzscore.pdf',sep ="_")
# YlGnBu = brewer.pal(9, "YlGnBu") 
# getPalette = colorRampPalette(YlGnBu)
# outlier.colors=c("NA", "#000000")
# 
# p = ggplot()
# p = p + geom_tile(data=top_overlap_all_zscore.m, aes(x=as.factor(sample), y=variable, fill=value), linetype="blank") + scale_fill_gradientn(colours=getPalette(100), na.value="white", name="relative modified z-score")# values=rescale(seq(0,6,by=12/99)))
# p = p + geom_tile(data=top_overlap_all_boolean.m, aes(x=as.factor(sample), y=variable, color=value), fill=NA, size=0.5) + scale_colour_manual(values = outlier.colors, name="KDE outlier")
# p = p + geom_text(data=top_overlap_all.m ,aes(x=as.factor(sample), y=variable, label = value), color="red", size=1.5, angle=90)
# p = p + xlab("sample") + ylab("top druggable targets") + theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=12), axis.text.y = element_text(colour="black", size=12))
# p
# ggsave(file=fn, width = 210, height = 400, units = "mm")


##### BRCA77 human data ##### 
###proteome###
BRCA77 = read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt')
BRCA77.d = BRCA77[row.names(BRCA77) %in% druggable,]
BRCA77_druggable = find_outlier(BRCA77.d, "BRCA77 druggable proteome", h=6, w=24)

###phosphoproteome###
BRCA77pho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
# 33239 phosphosites
row.names(BRCA77pho) = make.names(BRCA77pho$Gene.site,unique=T)
BRCA77pho.na = BRCA77pho[,-c(1,2)]
cat("Original number of BRCA77 ITRAQ phosphosites: 33239\n")
genes = sub(".NP.*", "", row.names(BRCA77pho.na))
BRCA77pho.na.d = BRCA77pho.na[genes %in% druggable, ] #651 phosphosites
row.names(BRCA77pho.na.d) = make.names(sub(".NP_\\d+_"," ",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("\\. _.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("\\._.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub(" _.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("_.*","",row.names(BRCA77pho.na.d)), unique=T)
cat("Druggable list filtered BRCA77 ITRAQ phosphosites: 651\n")

# find outliers
BRCA77pho_druggable = find_outlier(BRCA77pho.na.d, "BRCA77 druggable phosphoproteome", h=10, w=24)

# ### find extreme phosphos that are not indicated by pro 
# BRCA77pho_outlier_zscore.dt = as.data.frame(BRCA77pho_outlier_zscore)
# BRCA77pho_outlier_zscore.dt$Gene = sub("\\..*", "", row.names(BRCA77pho_outlier_zscore.dt))
# BRCA77pho_outlier_zscore.dt$Site = row.names(BRCA77pho_outlier_zscore.dt)
# BRCA77_outlier_zscore.dt = as.data.frame(BRCA77_outlier_zscore)
# BRCA77_outlier_zscore.dt$Gene = sub("\\..*", "", row.names(BRCA77_outlier_zscore.dt))
# 
# BRCA77pho_outlier_zscore.m = melt(BRCA77pho_outlier_zscore.dt, id.var = c("Gene","Site"))
# BRCA77_outlier_zscore.dt.m = melt(BRCA77_outlier_zscore.dt, id.var = "Gene")
# BRCA77_pro_pho_overlap = merge(BRCA77pho_outlier_zscore.m, BRCA77_outlier_zscore.dt.m, by = c("Gene","variable"), all.x=T)
# 
# cor(BRCA77_pro_pho_overlap$value.x, BRCA77_pro_pho_overlap$value.y, use = 'pairwise.complete.obs')
# 
# fn = paste(pd,'BRCA77_phospho_vs_pro_mzscore.pdf', sep="_")
# p = ggplot(data = BRCA77_pro_pho_overlap, aes(x=value.y, y=value.x, colour=variable, label=Site)) 
# p = p + geom_point(alpha=0.3) + xlab("druggable proteome z-score") + ylab("druggable phosphosites z-score") + theme_bw()
# p = p + xlim(-5,7.5) + ylim(-5,7.5) + geom_vline(xintercept=2,alpha=0.2) + geom_hline(yintercept=2,alpha=0.2)
# p = p + geom_text(aes(label=ifelse(value.x > 2.5 & value.x > value.y + 2 & Gene=="ERBB2",paste(variable,Site),'')),hjust=-0.05,just=0,size=3,alpha=0.7)
# p
# ggsave(file=fn, height=10, width=10)


##### test for enrichment in overlapping outlier #####
# limit ITRAQ outlier to the set that overlapped with LFQ
ITRAQ_L = ITRAQ_druggable$top_outlier_boolean[row.names(ITRAQ_druggable$top_outlier_boolean) %in% row.names(LFQ_druggable$top_outlier_boolean),]
ITRAQ_o = sum(ITRAQ_L, na.rm=T)
ITRAQ_no = dim(ITRAQ_L)[1]*dim(ITRAQ_L)[2]
LFQ_o = sum(LFQ_druggable$top_outlier_boolean, na.rm=T)
LFQ_no = dim(LFQ_druggable$top_outlier_boolean)[1]*dim(LFQ_druggable$top_outlier_boolean)[2]
overlap = 5

a = matrix(c(5,ITRAQ_o-5,LFQ_o-5,LFQ_no-(ITRAQ_o-5)-(LFQ_o-5)-5),nrow=2,ncol=2)
cat("\n##### Enrichment of overlapping iTRAQ and LFQ outliers#####\n")
fisher.test(a)

##### find outlier using joint TCGA-WHIM cohort #####
merged_proteome=merge(ITRAQ.d, BRCA77.d, by = "row.names") #63 genes
row.names(merged_proteome) = merged_proteome[,1]
merged_proteome = merged_proteome[,-1]
merged_druggable = find_outlier(merged_proteome, name="TCGA WHIM merged druggable proteome", whim_only=T)

sink(file=NULL)