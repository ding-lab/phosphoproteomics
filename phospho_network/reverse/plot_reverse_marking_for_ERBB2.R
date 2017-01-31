# plot ERBB2 as a example -------------------------------------------------
erbb2 <- k_p_her[k_p_her$kinase=="ERBB2",]
rownames(erbb2) <- erbb2$sample

IQR = quantile(erbb2$ave_sub_phos, probs=0.75, na.rm=T) - quantile(erbb2$ave_sub_phos, probs=0.25, na.rm=T) 
erbb2$phos_outlier = (erbb2$ave_sub_phos >= quantile(erbb2$ave_sub_phos, probs=0.75, na.rm=T) + 1.5*IQR)
erbb2$score <-  (erbb2$ave_sub_phos - quantile(erbb2$ave_sub_phos, probs=0.75, na.rm=T))/IQR

erbb2$pro_level <- vector(mode = "numeric", length = nrow(erbb2))+NaN
erbb2$pro_outlier <- vector(mode = "logical", length = nrow(erbb2))+NaN
for (sam in sample_names) {
  erbb2[sam,"pro_level"] <- pro_outlier["ERBB2",sam]
}
IQR = quantile(erbb2$pro_level, probs=0.75, na.rm=T) - quantile(erbb2$pro_level, probs=0.25, na.rm=T) 
erbb2$pro_outlier = (erbb2$pro_level >= quantile(erbb2$pro_level, probs=0.75, na.rm=T) + 1.5*IQR)

erbb2$Her2_status <- erbb2$Her_status
erbb2$Her2_status[erbb2$pro_outlier==TRUE] <- "pro_outlier"
erbb2$Her2_status[erbb2$pro_outlier==FALSE] <- "pro_normal"
erbb2$Her2_status <- paste(erbb2$Her_status,erbb2$Her2_status,sep = ",")

length(which(erbb2$Her2_status=="Her2-,pro_outlier"))


# y:ave_sub_phos
pd <- position_dodge(.65)
ggplot(data = erbb2,aes(x= kinase, y = ave_sub_phos, colour = Her2_status, shape = Her2_status)) +    
  geom_point(position = pd, size = 4) +
  scale_colour_manual(name = "Her2_status",
                      labels = c("Her2-,pro_normal", "Her2-,pro_outlier", "Her2+,pro_normal", "Her2+,pro_outlier"),
                      values = c("blue", "red", "blue", "red")) +   
  scale_shape_manual(name = "Her2_status",
                     labels = c("Her2-,pro_normal", "Her2-,pro_outlier", "Her2+,pro_normal", "Her2+,pro_outlier"),
                     values = c(19, 19, 17, 17))

## y:score
pd <- position_dodge(.65)
ggplot(data = erbb2,aes(x= kinase, y = score, colour = Her2_status, shape = Her2_status)) +    
  geom_point(position = pd, size = 4) +
  scale_colour_manual(name = "Her2_status",
                      labels = c("Her2-,pro_normal", "Her2-,pro_outlier", "Her2+,pro_normal", "Her2+,pro_outlier"),
                      values = c("blue", "red", "blue", "red")) +   
  scale_shape_manual(name = "Her2_status",
                     labels = c("Her2-,pro_normal", "Her2-,pro_outlier", "Her2+,pro_normal", "Her2+,pro_outlier"),
                     values = c(19, 19, 17, 17))

