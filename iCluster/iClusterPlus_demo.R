#setwd("/home/nfs/shenr/iClusterPLoSONEsoftware/")
library(iCluster)

data(breast.chr17)
fit=iCluster(datasets=breast.chr17, k=4, lambda=c(0.2,0.2))

plotiCluster(fit)

#a simple data transformation for nice image contrast
cn=breast.chr17[[2]]
cn[cn< -1.5]= -1.5;cn[cn>1.5]=1.5
exp=breast.chr17[[1]]
exp[exp< -2.5]= -2.5; exp[exp>2.5]=2.5

library(gplots)
library(lattice)
plotHeatmap(datasets=list(cn,exp),fit=fit)

#GBM data example including three data matrices: copy number, methylation, and expression.
data(gbm) 
data(coord)
chr=coord[,1]

#A variant of iCluster method with variance weighted shrinkage
fit=iCluster2(datasets=gbm, k=3, lambda=list(0.44,0.33,0.28))

plotiCluster(fit=fit)

compute.pod(fit)

plotHeatmap(fit=fit, datasets=gbm, feature.order=c(FALSE,TRUE,TRUE), sparse=c(FALSE,TRUE,TRUE),plot.chr=c(TRUE,FALSE,FALSE), chr=chr)

#Model tuning (using a simulated data example)
data(simu.datasets)
cv.fit=alist()
for(k in 2:5){
  cat(paste("K=",k,sep=""),'\n')
  cv.fit[[k]]=tune.iCluster2(datasets=simu.datasets, k,nrep=2, n.lambda=8)
}

#Reproducibility index (RI) plot
plotRI(cv.fit)


#Based on the RI plot, k=3 is the best solution
best.fit=cv.fit[[3]]$best.fit
#Try different color schemes
plotHeatmap(fit=best.fit,datasets=simu.datasets,sparse=c(TRUE,TRUE),col.scheme=list(bluered(256), greenred(256)))



# for (k in 1:10){
#   
#   fit.k = iClusterPlus(dt1=tBRCA_mut,dt2=tBRCA_CNV,dt3=tBRCA_RNA,dt4=tBRCA_PRO,
#                        type=c("binomial","gaussian","gaussian","gaussian"),
#                        lambda=c(0.05,0.20,1.00,0.80),K=k,maxiter=30)
#   
#   n = k + 1
#   #plotiCluster = function(fit = ){}
#   bw.col = colorpanel(2,low="white",high="black")
#   col.scheme = alist()
#   col.scheme[[1]] = bw.col
#   col.scheme[[2]] = bluered(256)
#   col.scheme[[3]] = bluered(256)
#   col.scheme[[4]] = bluered(256)
#   fn = paste(pd, n, "clusters_lambda005_020_100_080_HM.pdf", sep="_")
#   pdf(fn, useDingbats=FALSE)
#   plotHeatmap(fit=fit.k,datasets=list(tBRCA_mut,tBRCA_CNV,tBRCA_RNA,tBRCA_PRO),
#               type=c("binomial","gaussian","gaussian","gaussian"), col.scheme = col.scheme,
#               row.order=c(T,T,T,T),chr=chr,plot.chr=c(F,F,F,F),sparse=c(T,T,T,T),cap=c(F,F,F,F))
#   dev.off()
#   
#   fn = paste(pd, n, "-clusters_lambda005_020_100_080_ggHM.pdf", sep="_")
#   plotHeatmap_gg(fit=fit.k,datasets=list(tBRCA_clin,tBRCA_mut,tBRCA_CNV,tBRCA_RNA,tBRCA_PRO), fn = fn,
#                  dataset.names = c("CLINICAL","MUT","CNV","RNA","PRO"),
#                  type=c("clinical","binomial","gaussian","gaussian","gaussian"),
#                  row.order=c(F,T,T,T,T),sparse=c(F,T,T,T,T), row.names=c(T,T,F,F,F))
#   
#   #     clusters=fit.k$clusters
#   #     k=length(unique(clusters))
#   #     sorder=order(clusters)
#   #     row.names(tBRCA_mut)[sorder] # try to see if this ordering is correct
#   #     
#   #     # trouble shoot: use mut matrix as clin to see if key genes match up
#   #     figure = paste(pd,'LFQ_clin_proteome_naMax10_SDi2.pdf', sep="_")
#   #     plot_clin(LFQ_mp, clus_order2 , clin, figure) # change to just input a vector of samples
#   #     
#   BICs[[k]] = fit.k$BIC
#   dev.ratios[[k]] = fit.k$dev.ratio
# }
# 
# 
# BIC_m = do.call(rbind,BICs)
# dev.ratio_m = do.call(rbind,dev.ratios)
# 
# fn = paste(pd, 'icluster_k1-11_dev.ratio_summary.pdf',sep ="_")
# p = ggplot()
# p = p + geom_point(aes(x = c(1:11), y = c(0,dev.ratio_m))) + scale_x_discrete(breaks=c(1:11))
# p = p + labs(x="Number of clusters", y="% explained variation") + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=16))
# p 
# ggsave(file=fn, width=5, height=4, useDingbats=FALSE)
# 
# date()
# feature selection
# features = alist()
# features[[1]] = colnames(tBRCA_mut)
# features[[2]] = colnames(tBRCA_CNV)
# features[[3]] = colnames(tBRCA_RNA)
# features[[4]] = colnames(tBRCA_PRO)
# sigfeatures=alist()
# for(i in 1:4){
#   rowsum=apply(abs(fit.k$beta[[i]]),1, sum)
#   upper=quantile(rowsum,prob=0.75)
#   sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
#   sigfeatures[[i]]
# }
# names(sigfeatures)=c("mutation","copy number","expression","protein")
# #print a few examples of selected features
# sigfeatures[[1]]
# sigfeatures[[2]]
# sigfeatures[[3]]
# sigfeatures[[4]]

##### model tuning #####
#tune_icluster = function(){}
# set.seed(123)
# date()
# for(k in 1:5){
#   cv.fit = tune.iClusterPlus(cpus=12,dt1=tBRCA_mut,dt2=tBRCA_CNV,dt3=tBRCA_RNA,dt4=tBRCA_PRO,
#                              type=c("binomial","gaussian","gaussian","gaussian"),
#                              scale.lambda=c(0.04,0.90,0.90,0.90),n.lambda=307,maxiter=20)
#   save(cv.fit, file=paste("clusterRdata/cv.fit.k",k,".Rdata",sep=""))
# }
# date()


# # model selection
# output=alist()
# files=grep("cv.fit",dir())
# for(i in 1:length(files)){
#   load(dir()[files[i]])
#   output[[i]]=cv.fit
#   }
# nLambda = nrow(output[[1]]$lambda)
# nK = length(output)
# BIC = getBIC(output)
# devR = getDevR(output)
# 
# minBICid = apply(BIC,2,which.min)
# devRatMinBIC = rep(NA,nK)
# for(i in 1:nK){
#   devRatMinBIC[i] = devR[minBICid[i],i]
#   }
# 
# plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
#        ylab="%Explained Variation")
# 
# chr=unlist(strsplit(colnames(gbm.cn),"\\."))
# chr=chr[seq(1,length(chr),by=2)]
# chr=gsub("chr","",chr)
# chr=as.numeric(chr)

# #truncate the values for a better image plot
#   cn.image=gbm.cn
# cn.image[cn.image>1.5]=1.5
# cn.image[cn.image< -1.5]= -1.5
# exp.image=gbm.exp
# exp.image[exp.image>2.5]=2.5
# exp.image[exp.image< -2.5]= -2.5

# ##### example data #####
# data(gbm)
# 
# # select mutation features
# dim(gbm.mut)
# mut.rate=apply(gbm.mut,2,mean)
# gbm.mut2 = gbm.mut[,which(mut.rate>0.02)]
# 
# gbm.cn=gbm.cn[order(rownames(gbm.cn)),]
# # check if all the samples are in the same order for the three data sets
# all(rownames(gbm.mut2)==rownames(gbm.exp))
# 
# fit.single=iClusterPlus(dt1=gbm.mut2,dt2=gbm.exp,
#                           type=c("binomial","gaussian"),
#                           lambda=c(0.04,0.90),K=2,maxiter=10)
# 
# # quick plot to check
# bw.col = colorpanel(2,low="white",high="black")
# col.scheme = alist()
# col.scheme[[1]] = bw.col
# col.scheme[[2]] = bluered(256)
# plotHeatmap(fit=fit.single,datasets=list(gbm.mut2,gbm.exp),
#             type=c("binomial","gaussian"), col.scheme = col.scheme,
#             row.order=c(F,T),chr=chr,plot.chr=c(F,F),sparse=c(T,T),cap=c(F,F))
# 
# # model tuning
# set.seed(123)
# date()
# for(k in 1:5){
#   cv.fit = tune.iClusterPlus(cpus=12,dt1=gbm.mut2,dt2=gbm.cn,dt3=gbm.exp,
#                                type=c("binomial","gaussian","gaussian"),K=k,n.lambda=185,
#                                scale.lambda=c(1,1,1),maxiter=20)
#   save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
#   }
# date()
# 
# model selection
# output=alist()
# files=grep("cv.fit",dir())
# for(i in 1:length(files)){
#   load(dir()[files[i]])
#   output[[i]]=cv.fit
# }
# nLambda = nrow(output[[1]]$lambda)
# nK = length(output)
# BIC = getBIC(output)
# devR = getDevR(output)
# 
# minBICid = apply(BIC,2,which.min)
# devRatMinBIC = rep(NA,nK)
# for(i in 1:nK){
#   devRatMinBIC[i] = devR[minBICid[i],i]
# }
# 
# plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
#        ylab="%Explained Variation")
# 
# chr=unlist(strsplit(colnames(gbm.cn),"\\."))
# chr=chr[seq(1,length(chr),by=2)]
# chr=gsub("chr","",chr)
# chr=as.numeric(chr)
# #truncate the values for a better image plot
#   cn.image=gbm.cn
# cn.image[cn.image>1.5]=1.5
# cn.image[cn.image< -1.5]= -1.5
# exp.image=gbm.exp
# exp.image[exp.image>2.5]=2.5
# exp.image[exp.image< -2.5]= -2.5
# 
# # feature selection
# features = alist()
# features[[1]] = colnames(gbm.mut2)
# features[[2]] = colnames(gbm.cn)
# features[[3]] = colnames(gbm.exp)
# sigfeatures=alist()
# for(i in 1:3){
#   rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
#   upper=quantile(rowsum,prob=0.75)
#   sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
#   }
# names(sigfeatures)=c("mutation","copy number","expression")
# #print a few examples of selected features
#   head(sigfeatures[[1]])
# 
# # plot
# bw.col = colorpanel(2,low="white",high="black")
# col.scheme = alist()
# col.scheme[[1]] = bw.col
# col.scheme[[2]] = bluered(256)
# col.scheme[[3]] = bluered(256)
# plotHeatmap(fit=best.fit,datasets=list(gbm.mut2,cn.image,exp.image),
#               type=c("binomial","gaussian","gaussian"), col.scheme = col.scheme,
#               row.order=c(F,F,T),chr=chr,plot.chr=c(F,T,F),sparse=c(T,F,T),cap=c(F,T,F))