##### iCluster.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# cross level clustering between CNV, RNA, and Proteome data

### resources
setwd("/gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer/pan3can_analysis/iCluster")
source("~/bin/LIB_exp.R")
baseD = "/gscmnt/gc2524/dinglab/Proteomics/projects/CPTAC_pan3Cancer"
cancer_genes = read.table(file='/gscmnt/gc2524/dinglab/Proteomics/projects/reference_files/cgenes_and_druggable.list', header=FALSE, stringsAsFactors = F)
load("clusterRdata/2015-11-13_BRCA.Rdata")

if (FALSE){ # for use on my macpro
  setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/iCluster")
  baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"
  source("/Users/khuang/bin/LIB_exp.R")
  cancer_genes = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/cgenes_and_druggable.list', header=FALSE, stringsAsFactors = F)
}
cgenes = as.vector(t(cancer_genes))

# libraries
library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)
library(gtable)

# # later: input number of k, optimize each k and store the temporary result, separate into an optimization and a plotting script
# args=commandArgs(TRUE)
# k=args[1]


### functions
unfactorize = function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

get.clinical.scale = function() {
  # Set1 colors
  #colors = c("#f0f0f0", "#636363", "#f0f0f0", "#636363", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  # use Perou's intrinsic subtype colors instead
  colors = c(NA, "#101010", NA, "#636363", "red", "pink", "#ffeda0", "#deebf7", "#08306b", "#FFFF33", "#A65628", "#F781BF", "#999999") #positive is dark grey       
  
  color.names = c("wt","mut","negative", "positive", "Basal", "Her2", "CLDN low", "LumB","LumA")
  names(colors) = color.names
  clinical.color.scale = scale_fill_manual(name="mutation", values=colors)
  
  return(clinical.color.scale)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

rbind_gtable <- function(...){
  
  gtl = lapply(list(...), ggplotGrob)
  
  bind2 <- function (x, y) 
  {
    stopifnot(ncol(x) == ncol(y))
    if (nrow(x) == 0) 
      return(y)
    if (nrow(y) == 0) 
      return(x)
    y$layout$t <- y$layout$t + nrow(x)
    y$layout$b <- y$layout$b + nrow(x)
    x$layout <- rbind(x$layout, y$layout)
    x$heights <- gtable:::insert.unit(x$heights, y$heights)
    x$rownames <- c(x$rownames, y$rownames)
    x$widths <- grid::unit.pmax(x$widths, y$widths)
    x$grobs <- append(x$grobs, y$grobs)
    x
  }
  
  Reduce(bind2, gtl)
}

theme1 = function(...) theme( panel.background = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.margin = unit(0,"null"),
                              axis.ticks = element_blank(),
                              axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              axis.ticks.length = unit(0,"null"),
                              axis.ticks.margin = unit(0,"null"),...)

# use to as an alternative to the plotHeatmap function in iCluster package
# use ggplot instead of lattice, allows for flexibility
if (FALSE){ # troubleshooting defaults
  datasets=list(tBRCA_mut,tBRCA_CNV,tBRCA_RNA,tBRCA_PRO);
  dataset.names = c("MUT","CNV","RNA","PRO")
  type=c("binomial","gaussian","gaussian","gaussian")
  row.order=c(F,T,T,T);sparse=c(T,T,T,T)
}
plotHeatmap_gg = function(fit, datasets, dataset.names, type=c("gaussian","binomial","poisson","multinomial"), 
                          sample.order=NULL, row.order=NULL, row.names = NULL, sparse=NULL, threshold=rep(0.25,length(datasets)), 
                          fn="iCluster_heatmap_gg.pdf", height=10, width=14){
  
  m=length(datasets)  
  if(is.null(row.order)){row.order=rep(T,m)}  
#   if(is.null(scale)){scale=rep("none",m)}
#   if(is.null(sparse)){sparse=rep(F,m)}
#   if(is.null(cap)){cap=rep(F,m)}
#   if(is.null(plot.chr)){plot.chr=rep(F,m)}
  
  #get clusters 
  clusters=fit$clusters
  
  k=length(unique(clusters))
  if(is.null(sample.order)){sorder=order(clusters)}else{sorder=sample.order}
  m=length(datasets)
  pp=unlist(lapply(1:m,function(l){dim(datasets[[l]])[2]})) # number of observations per dataset
  n=dim(datasets[[1]])[1] # number of samples
  
  
  #cluster divider
  a=clusters[sorder]
  l=length(a)
  brkpoints=which(a[2:l]!=a[1:(l-1)])
  cluster.start=c(1,brkpoints+1)
  
#   my.panel.levelplot <- function(...) {
#     panel.levelplot(...)
#     panel.abline(v=(cluster.start[-1]-0.5), col="black", lwd=1,lty=1)
#     panel.scales=list(draw=FALSE)
#   }
  
  # list to store plotting data from different levels
  plots = list()
  
  for(i in 1:m){
    
    if (type[i] != "clinical"){
      rowsum=apply(abs(fit$beta[[i-1]]),1, sum)
      if(sum(rowsum)==0)warning(paste("All Lasso coefficients are zero for data type",i), call. = FALSE)
      
      if(mean(rowsum>0)>threshold[i]){upper=quantile(rowsum,prob=(1-threshold[i]))}else{upper=0}
      
      if(sparse[i]==T & sum(rowsum> upper)>1){
        image.data=datasets[[i]][sorder,which(rowsum> upper)]
      }else{image.data=datasets[[i]][sorder,]}
      if(row.order[i]==T){
        diss=1-cor(image.data,use="na.or.complete");
        hclust.fit=hclust(as.dist(diss));
        #hclust.fit=hclust(dist(t(image.data)))
        gorder=hclust.fit$order
        image.data=image.data[,gorder]
      }
    } else{ # for clinical data
      image.data=datasets[[i]][sorder,]
    }

    #     # scaling: may be useful but not using for our purpose
    #     scales=list(draw=F)   
    #     scale.fn=function(x){
    #       x <- sweep(x, 1L, rowMeans(x, na.rm = T), check.margin = T)
    #       sx <- apply(x, 1L, sd, na.rm = T)
    #       x<- sweep(x, 1L, sx, "/", check.margin = T)
    #       return(x)
    #     }
    #     if(scale[i]=='row'){image.data=scale.fn(image.data)}
    #     if(scale[i]=='col'){image.data=scale.fn(t(image.data)); image.data=t(image.data)}
    
    image.data2=as.matrix(rev(as.data.frame(image.data)))   #reverse the rows so the image will be top-down ordering
    image.data_m = melt(image.data2)
    colnames(image.data_m) = c("sample","gene","value")
    #     image.data_m$name = dataset.names[i]
    #     melted_image.data[[i]] = image.data_m
    
    # plot
    if(type[i]=="clinical"){
      image.data_m$value = as.factor(image.data_m$value)
      p = ggplot(image.data_m)
      p = p + geom_tile(aes(x=sample, y=gene, fill=value), color=NA)#, linetype="blank") 
      p = p + get.clinical.scale()
      p = p + theme_bw() 
      p = p + theme1() #+ coord_equal()
      p = p + ylab(dataset.names[i]) + theme(axis.title.y = element_text(size=16, angle=90))
      if (row.names[[i]]){p = p + theme(axis.text.y = element_text(colour="black", size=10))}
      #p = p + geom_vline(xintercept = cluster.start - 0.5)
    } else if(type[i]=="binomial"){
      image.data_m$value = as.factor(image.data_m$value)
      p = ggplot(image.data_m)
      p = p + geom_tile(aes(x=sample, y=gene, fill=value), color=NA)#, linetype="blank") 
      p = p + scale_fill_manual(values=c(NA, "black")) #+ scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024)
      p = p + theme_bw() 
      p = p + theme1() #+ coord_equal()
      p = p + ylab(dataset.names[i]) + theme(axis.title.y = element_text(size=16, angle=90))
      if (row.names[[i]]){p = p + theme(axis.text.y = element_text(colour="black", size=10))}
      p = p + geom_vline(xintercept = cluster.start - 0.5)
    }else{
      max_bound = quantile(image.data_m$value,prob=0.99,na.rm=T)
      min_bound = quantile(image.data_m$value,prob=0.01,na.rm=T)
      image.data_m$truncated_value = image.data_m$value
      image.data_m[image.data_m$truncated_value >= max_bound,]$truncated_value = max_bound
      image.data_m[image.data_m$truncated_value <= min_bound,]$truncated_value = min_bound
      
      p = ggplot(image.data_m)
      p = p + geom_tile(aes(x=sample, y=gene, fill= truncated_value), linetype="blank") + scale_fill_gradientn(name= "Levels", na.value=NA, colours=RdBu1024, limit=c(min_bound,max_bound))
      p = p + theme_bw() 
      p = p + theme1() #+ coord_equal()
      p = p + ylab(dataset.names[i]) + theme(axis.title.y = element_text(size=16, angle=90))
      if (row.names[[i]]){p = p + theme(axis.text.y = element_text(colour="black", size=10))}
      p = p + geom_vline(xintercept = cluster.start - 0.5)
    }
    
    plots[[i]] = p
  }
  
  # bind the plots
  gp = do.call(rbind_gtable, plots)
  # print the integrated plot
  grid.newpage()
  pdf(fn, height=height, width=width,useDingbats = F)
  grid.draw(gp)
  dev.off()
}

##### MAIN CODE #####

##### loop through different ks here ####
BICs = vector("list")
dev.ratios = vector("list")

set.seed(123)
date()

# model selection
output=alist()
files=grep("cv.fit",dir("clusterRdata"))
for(i in 1:length(files)){
  load(paste("clusterRdata/", dir("clusterRdata")[files[i]], sep=""))
  output[[i]]=cv.fit
  
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}

### plot percentage variance explained by different Ks
# plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
#        ylab="%Explained Variation")
fn = paste(pd, 'icluster_Ks_dev.ratio_summary.pdf',sep ="_")
p = ggplot()
p = p + geom_point(aes(x = c(1:(nK+1)), y = c(0,devRatMinBIC))) + scale_x_discrete(breaks=c(1:nK+1))
p = p + labs(x="Number of clusters", y="% explained variation") + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=16))
p 
ggsave(file=fn, width=5, height=4, useDingbats=FALSE)

for (k in 1:nK){
  best.fit=output[[k]]$fit[[which.min(BIC[,k])]]
  n=k+1
  fn = paste(pd, n, "-clusters_lambda005_020_100_080_top_ggHM.pdf", sep="_")
  plotHeatmap_gg(fit=fit.k,datasets=list(tclin,tmut,tCNV,tRNA,tPRO), fn = fn,
                 dataset.names = c("CLINICAL","MUT","CNV","RNA","PRO"),
                 type=c("clinical","binomial","gaussian","gaussian","gaussian"),
                 row.order=c(F,T,T,T,T),sparse=c(F,T,T,T,T), row.names=c(T,T,F,F,F))
}


date()

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
