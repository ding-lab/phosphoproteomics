### scatterplot_regression_fit.R ### 
# Yige Wu @ WashU 2016 Nov
# scatterplot+regression fit line for kinase/subsrate pair significant regression p-value
# y = substrate phosphorylation at phosphosite
# x = kinase protein expression level
# todo: add another facet for cancer
# todo: add p-value

# names mentioned in kuan's proteomics manuscript
#names <- c("BRCA1","BRCA2","CDH1","CDK12","GATA3","KMT2C","KRAS","MAP3K1","NF1","PIK3CA","TP53","APC","SMAD4","CHEK2","MAPK3","AKT1")
# > intersect(names,table_sig$KINASE)
# [1] "MAPK3" "AKT1" 

# load BRCA_3model.RData

library(ggplot2)
library(grid)
require(plyr)

lm_eqn = function(df){
  m = glm(pho_sub_norm ~ pro_kinase, data = df);
  eq <- substitute(italic(pho_sub_norm) == a + b %.% italic(pro_kinase),
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2)));
  as.character(as.expression(eq));
}
table <- table1

# choose kinase
kinase = "MAPK3"
kinase = "AKT1"
#kinase = "PAK4"
#kinase = "ERBB2"

#kinase="CDK12"; 
#kinase="AKT1"; kinase= "BRCA1"; 
#kinase="TP53"; 
#kinase= "CHEK2"; kinase="APC"; kinase=; kinase="BRCA2"; kinase="MAP3K1" ; kinase="NF1"; kinase="MAPK3" 

# extract the row numbers
rows <- which(table$KINASE == kinase & table$P_pro_kin <= 0.05 & table$model == "pho_sub~pro_kinase")
#rows <- which(table$KINASE == kinase & table$size >= 50 & table$P_pro_kin <= 0.05 & table$model == "pho_sub~pro_kinase")
#rows <- c(7070, 7000, rows[1:11])
rows <- which(table$KINASE == kinase & table$SUBSTRATE == kinase & table$P_pro_kin <= 0.05 & table$model == "pho_sub~pro_kinase")

data_k <- data.frame()
for(i in rows ){
  # extract the data for each pair
  data <- list1[[i]]$DATA
  data$SUBSTRATE <- list1[[i]]$SUBSTRATE
  data$SUB_MOD_RSD <- list1[[i]]$SUB_MOD_RSD
  data_k <- rbind(data_k,data)
}
data_k <- data_k[complete.cases(data_k),]

eq <- ddply(data_k,.(SUBSTRATE,SUB_MOD_RSD),lm_eqn)

p <- ggplot(data = data_k, aes(x = pro_kinase, y = pho_sub_norm)) +
  geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ggtitle(paste("kinase = ",kinase,sep=""))+
  ylim(c(0,1))

p1 = p + facet_grid(SUBSTRATE+SUB_MOD_RSD~.)
p2 = p1 + geom_text(data=eq,aes(x = 0.5*max(data_k$pro_kinase), y = 0.9,label=V1, size=5), parse = TRUE, inherit.aes=FALSE)
p2
#BRCA_cor_pro_kinase_phosphosite_

## one kinase, one substrate

kinase = "CDK2";
sub = "APC";
rows <- which(table$KINASE == kinase & table$SUBSTRATE == sub & table$P_pro_kin <= 0.05 & table$model == "pho_sub~pro_kinase")

data_k <- data.frame()
for(i in rows ){
  # extract the data for each pair
  data <- list1[[i]]$DATA
  data$SUBSTRATE <- list1[[i]]$SUBSTRATE
  data$SUB_MOD_RSD <- list1[[i]]$SUB_MOD_RSD
  data_k <- rbind(data_k,data)
}
data_k <- data_k[complete.cases(data_k),]

eq <- ddply(data_k,.(SUBSTRATE,SUB_MOD_RSD),lm_eqn)

p <- ggplot(data = data_k, aes(x = pro_kinase, y = pho_sub_norm)) +
  geom_smooth(method = "glm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ggtitle(paste("kinase = ",kinase,sep=""))+
  ylim(c(0,1))

p1 = p + facet_grid(SUBSTRATE+SUB_MOD_RSD~.)
p2 = p1 + geom_text(data=eq,aes(x = 0.5*max(data_k$pro_kinase), y = 0.9,label=V1, size=5), parse = TRUE, inherit.aes=FALSE)
p2