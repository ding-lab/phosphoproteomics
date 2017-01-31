# score~pro_kin for RTKs --------------------------------------------------
RTK_file = read.table(paste(baseD,"pan3can_shared_data/reference_files/RTKs_list.txt",sep = ""))
rtk <- as.vector(t(RTK_file))

row <- c()
for (kinase in rtk) {
  temp <- which(overlap$kinase==kinase & !is.na(overlap$score) & !is.na(overlap$pro_level) )
  row <- c(row,temp)
}
table <- overlap[row,]
table$transparency <- vector(mode = "numeric", length = nrow(table))
table$transparency[table$pro_outlier_status] <- table$transparency[table$pro_outlier_status] + 1
table$transparency[table$missense] <- table$transparency[table$missense] + 1

ggplot(data = table,aes(x= pro_level, y = score, colour = kinase , shape = missense.pro_outlier, alpha = (transparency+1)/3)) +    
  geom_point(size = 2) + 
  geom_text(data=subset(table, pro_outlier_status | missense ),
            aes(x= pro_level, y = score,label=kinase),hjust = 0, nudge_x = -0.3 ,nudge_y = 0.1)
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/RTK_BRCA_score~pro_level.png',sep ="")
ggsave(file=fn, height=8, width=10)

ggplot(data = table,aes(x= pro_level, y = score, shape = missense.pro_outlier, alpha = (transparency+3)/5)) +    
  geom_point(size = 2) + 
  geom_text(data=subset(table, aa != "wt" ),
            aes(x= pro_level, y = score,label=aa), hjust = 0, nudge_x = -1.2 ,nudge_y = 0.2, size = 2) +
  facet_wrap(~ kinase, nrow = 4)
fn = paste(baseD,'pan3can_shared_data/analysis_results/reverse/RTK_BRCA_score~pro_level_aa_change.png',sep ="")
ggsave(file=fn, height=10, width=13)


