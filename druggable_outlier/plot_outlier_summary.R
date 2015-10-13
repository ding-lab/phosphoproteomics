##### plot_outlier_summary.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# run outlier analysis for 3 cancer types and plot the result

##### dependencies #####
setwd("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier")
source("/Users/khuang/bin/LIB_exp.R")
#source("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier/outlier.R")
# system("mkdir logs")
# logFile = paste("logs/", date, "_outlier_analysis.log", sep="")
# sink(file=logFile)

summary = read.table(header=T, sep="\t","/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/druggable_outlier/results/2015-10-10_outlier_summary.txt")

fn = paste(pd, 'pan3can_druggable_outliers_summary.pdf',sep ="_")
YlOrRd = brewer.pal(9, "YlOrRd") 
getPalette = colorRampPalette(YlOrRd)
outlier.colors=c("NA", "#000000")

summary$Level = factor(summary$Level, levels = c("CNV", "RNA", "PRO", "PHO", "GLY"))

#levels(summary$Level) = c("CNV", "RNA", "PRO", "PHO", "GLY")
summary$rounded_outlier_percentage = round(as.numeric(summary$Outlier_precentage), digits=1)
summary$Gene = factor(summary$Gene, levels = summary$Gene[order(summary$rounded_outlier_percentage)])
p = ggplot(data=summary)
p = p + facet_grid(.~Cancer)
#p = p + coord_equal()
p = p + geom_tile(aes(x=Level, y=Gene, fill=rounded_outlier_percentage), linetype="blank") + scale_fill_gradientn(name= "Percentage", colours=getPalette(100), na.value=NA, limit =c(0,22))
p = p + geom_text(aes(x=Level, y=Gene, label = rounded_outlier_percentage, stringsAsFactors=FALSE), color="black", size=3)
p = p  + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=5),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
ggsave(file=fn, height=10, width=15, useDingbats=FALSE)

##### short summary ##### 
# note this is a bit problematic at this point as didn't print any values lower than 1.5
fn = paste(pd, 'pan3can_druggable_outliers_summary_short.pdf',sep ="_")
YlOrRd = brewer.pal(9, "YlOrRd") 
getPalette = colorRampPalette(YlOrRd)
outlier.colors=c("NA", "#000000")

# showing less genes
summary2 = summary[summary$rounded_outlier_percentage>=1.5,]

#summary2$Level = factor(summary$Level, levels = c("CNV", "RNA", "PRO", "PHO", "GLY"))

#levels(summary$Level) = c("CNV", "RNA", "PRO", "PHO", "GLY")
#summary$rounded_outlier_percentage = round(as.numeric(summary$Outlier_precentage), digits=1)
#summary2$Gene = factor(summary2$Gene, levels = summary2$Gene[order(summary2$rounded_outlier_percentage)])
p = ggplot(data=na.omit(summary2))
p = p + facet_grid(.~Cancer)
p = p + geom_tile(aes(x=Level, y=Gene, fill=rounded_outlier_percentage), linetype="blank") + scale_fill_gradientn(name= "Percentage", colours=getPalette(100), na.value=NA, limit =c(0,22))
p = p + geom_text(aes(x=Level, y=Gene, label = rounded_outlier_percentage, stringsAsFactors=FALSE), color="black", size=3)
p = p  + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=5),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
ggsave(file=fn, height=10, width=15, useDingbats=FALSE)