### intersect_BRCA_OV.R ### 
# Yige Wu @ WashU 2016 Nov
# interest the kinase/phosphosite pairs with significant regression results in BRCA and OV datasets

## get significant results for BRCA and OV
# BRCA
table_BRCA_sig <- table_BRCA[table_BRCA$P_pro_kin <= 0.05,]
# OV
table_OV_sig <- table_OV[table_OV$P_pro_kin <= 0.05,]

## interset kinase/phosphosite pairs for each model between BRCA and OV
k_ps_intersect <- merge(table_BRCA_sig[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","model")],table_OV_sig[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","model")])
k_ps_int_k <- levels(as.factor(as.character(k_ps_intersect$KINASE)))
as.factor(as.character(k_ps_intersect$SUBSTRATE))
# model1
k_ps_int1 <- k_ps_intersect[k_ps_intersect$model=="pho_sub~pro_kinase",]
as.factor(as.character(k_ps_int1$KINASE))
as.factor(as.character(k_ps_int1$SUBSTRATE))

# model2
k_ps_int2 <- k_ps_intersect[k_ps_intersect$model=="pho_sub~pro_kinase+pro_sub",]
as.factor(as.character(k_ps_int2$KINASE))
as.factor(as.character(k_ps_int2$SUBSTRATE))

# model3
k_ps_int3 <- k_ps_intersect[k_ps_intersect$model=="pho_sub~pro_kinase+pro_sub+pho_kin_grouped",]
as.factor(as.character(k_ps_int3$KINASE))
as.factor(as.character(k_ps_int3$SUBSTRATE))
