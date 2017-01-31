# extract the subtype phosphorylation level according to bubble ch --------
pho_subtype_mean$SUBSTRATE <- pho_rsd_split$SUBSTRATE
pho_subtype_mean$SUB_MOD_RSD <- pho_rsd_split$SUB_MOD_RSD

pho_top <- pho_subtype_mean[order(pho_subtype_mean$Basal,decreasing = TRUE),]


pho_all <- c()
for (cohort in c("Her2","LumA","LumB","Basal")) {
  temp <- data.frame(pho_subtype_mean[,cohort])
  temp <- cbind(temp,pho_rsd_split[,c("SUBSTRATE","SUB_MOD_RSD")])
  temp$cohort <- cohort
  colnames(temp) <- c("pho_subtype","SUBSTRATE","SUB_MOD_RSD","cohort")
  pho_all <- rbind(pho_all,temp)
}
pho_all$X <- "average"
# plot heatmap for subtype phosphorylation level --------------------------
lim = max(abs(max(pho_all$pho_subtype, na.rm = TRUE)),abs(min(pho_all$pho_subtype,na.rm = T)))
plot4 = ggplot(pho_all)
plot4 = plot4 + geom_tile(aes(x=X, y=SUB_MOD_RSD, fill=pho_subtype), color=NA)#, linetype="blank") 
plot4 = plot4 + facet_grid(SUBSTRATE~cohort, drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
plot4 = plot4 + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
plot4 = plot4 + theme_bw() 
plot4 = plot4 + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
plot4