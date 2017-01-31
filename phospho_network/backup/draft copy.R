substrate <- "PPP1CA";
rsd <- "T320"
kinase <- "PRKCA"; 

table <- test[test$kinase==kinase,]
table$pro_level <- t(pro_data[pro_data$X==kinase,table$sample])
table$phosite <- t(pho_data[pho_rsd_split$SUBSTRATE== substrate & pho_rsd_split$SUB_MOD_RSD==rsd,table$sample])

pd <- position_dodge(.65)
ggplot(data = table,aes(x= pro_level, y = phosite , colour = missense.pro_outlier, shape = missense.pro_outlier)) +    
  geom_point(position = pd, size = 4) +
  scale_colour_manual(name = "missense.pro_outlier",
                      labels = c("wt,pro_normal", "wt,pro_outlier", "missense,pro_normal", "missense,pro_outlier"),
                      values = c("blue", "red", "blue", "red")) +   
  scale_shape_manual(name = "missense.pro_outlier",
                     labels = c("wt,pro_normal", "wt,pro_outlier", "missense,pro_normal", "missense,pro_outlier"),
                     values = c(19, 19, 17, 17))
dev.off()

