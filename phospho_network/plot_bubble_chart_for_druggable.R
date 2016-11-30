### plot druggable ones for non-self-regualted
selected <- names #need to input the druggable gene names first
cancer <- "BRCA" #choose the cancer you want to sort
t0 <- table_sig_other

rows_k <- c()
for(kinase in selected){
  r <- unlist(which(t0$KINASE==kinase))
  rows_k <- c(rows_k,r)
}

rows_s <- c()
for(sub in selected){
  r <- unlist(which(t0$SUBSTRATE==sub))
  rows_s <- c(rows_s,r)
}

t1 <- t0[rows_s,]

t1 <- data.table(t1[t1$model=="pho_sub~pro_kinase" & t1$Cancer == cancer,], key="P_pro_kin")

rows <- c()
for(i in 1:top){
  r <- unlist(which(t0$pair==t1$pair[i]))
  rows <- c(rows,r)
}

table <- t0[rows,]
table$order <- nrow(table):1
table$pairs <- reorder(table$pair,table$order)
p = ggplot(table,aes(x=model, y=pairs, colour=-log10(P_pro_kin), size =size))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point() + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
p