# input -------------------------------------------------------------------
library(ggplot2)
setwd("~/Desktop/Figures/volcano_plots")
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Initiate ----------------------------------------------------------------
cancer <- "OV"
mode <- "pho_sub~pro_kinase"
cis <- FALSE

# plot all pairs ----------------------------------------------------------
table <- table_2can[table_2can$Cancer==cancer & table_2can$model==mode & table_2can$self==cis,]
p <- qplot(coef_pro_kin, -log10(FDR_pro_kin), data = table)
p


# plot pairs without outliers ----------------------------------------------------------
x1 <- remove_outliers(table$coef_pro_kin)
p <- qplot(coef_pro_kin, -log10(FDR_pro_kin), data = table[table$coef_pro_kin==x1,])
p


# plot significant pairs --------------------------------------------------
p <- qplot(coef_pro_kin, -log10(FDR_pro_kin), data = table[table$FDR_pro_kin <= sig,])
p
