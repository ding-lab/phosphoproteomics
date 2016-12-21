### common libs and dependencies ### 
# especially for plotting #

# dependencies
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library("gplots")
library(gridExtra)
library(matrixStats)
library("WGCNA")
library(grid)
# library("PerformanceAnalytics")
# library("lattice")
# library(VennDiagram)
# library(grDevices)

# mis
system("mkdir figures")
date=Sys.time()
date = sub(" .*","",date)
figurePath=paste("figures/", date, "/",sep="")
system(paste("mkdir ", figurePath, sep=""))
pd = paste(figurePath, date, "_KH", sep="")

system("mkdir tables")
tablePath=paste("tables/", date, "/",sep="")
system(paste("mkdir ", tablePath, sep=""))
td = paste(tablePath, date, "_KH", sep="")

col_paletteB = colorRampPalette(brewer.pal(9,"Blues"))
col_paletteR = colorRampPalette(brewer.pal(9,"Reds"))
RdBu = brewer.pal(9, "RdBu") 
getPalette = colorRampPalette(rev(RdBu))
RdBu1024 = colorRampPalette(rev(RdBu))(1024)
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette2 = colorRampPalette(YlGnBu)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
set1 = brewer.pal(9,"Set1")
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals) ))

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

theme0 = function(...) theme( legend.position = "none",
                               panel.background = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               #panel.margin = unit(0,"null"),
                               axis.ticks = element_blank(),
                               axis.text.x = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               axis.ticks.length = unit(0,"null"),
                               axis.ticks.margin = unit(0,"null"),
                               panel.border=element_rect(color=NA),...)