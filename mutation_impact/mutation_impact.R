##### mutation_impact.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# find druggable protein downstream of driver mutations
# called by other scripts
# mutation_impact.R 

#setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/mutation_impact/")
#source("/Users/khuang/bin/LIB_exp.R")

diff_exp = function(m, g1, g2, name="data", limma=FALSE){
  g1.n = deparse(substitute(g1))
  g2.n = deparse(substitute(g2))
  m.n = deparse(substitute(m))
  m=as.matrix(m)
  
  # for limma
    g1g2 = c(g1,g2)
    m_g1g2 = m[,colnames(m) %in% g1g2]
    g2.b = as.numeric(colnames(m_g1g2) %in% g2)
    #m_g1g2 = m_g1g2[rowSums(!is.na(m_g1g2[,g1]))>0 & rowSums(!is.na(m_g1g2[,g2]))>0,]

  # for t-test and Wilcoxon test
  # print header
  x = m_g1g2[,colnames(m_g1g2) %in% g1] 
  y = m_g1g2[,colnames(m_g1g2) %in% g2]
  
  stats = matrix(,nrow=nrow(x),ncol=7)
  row.names(stats)=row.names(x)
  colnames(stats) = c(paste(m.n,g1.n,"mean",sep="_"), paste(m.n,g2.n,"mean",sep="_"), "meanChange","t_test_Tstat", "t_test_p", "w_test_Wstat", "w_test_p")
  for (i in 1:nrow(x)){
    x1 = as.numeric(x[i,])
    y1 = as.numeric(y[i,])
    meanx = mean(x1, na.rm=T)
    meany = mean(y1, na.rm=T)
    change = meanx-meany # using iTRAQ data this is analogus to log fold change; as log(exp1/ref) - log(exp2/ref) = log(exp1/exp2)
    if (sum(!is.na(x1))<2 | sum(!is.na(y1))<2){
      stats[i,]=c(meanx,meany,change,rep("NA",4))} 
    else{
      # t-test
      t = try(t.test(x1,y1), silent=T)
      if (is(t,"try-error")){ stats[i,] = c(meanx,meany,change,rep("NA",4)); next}
      t.p = t$p.value
      t.tstat = t$statistic
      t.conf1 = t$conf.int[1]
      t.conf2 = t$conf.int[2]
      # Wilcoxon Rank Sum Test
      w = wilcox.test(x1,y1)
      w.p = w$p.value
      w.Wstat = w$statistic
      # return the results
      stats[i,] = c(meanx, meany, change, t.tstat, t.p, w.Wstat, w.p)
    }
  }
  t_test_fdr=p.adjust(stats[,"t_test_p"], method="BH")
  stats=cbind(stats, t_test_fdr)
  w_test_fdr=p.adjust(stats[,"w_test_p"], method="BH")
  stats=cbind(stats, w_test_fdr)
  if (limma){
    n_genes = dim(m_g1g2)[1]
    fit = lmFit(m_g1g2, design=g2.b)
    fit = eBayes(fit)
    limma_t = topTable(fit, adjust.method="BH", number = n_genes, p.value=1.1)
    stats = merge(stats, limma_t, by = "row.names")
    row.names(stats) = stats$Row.names
    stats=stats[,-1]
  }
  #stats=stats[order(stats$t_test_fdr, stats$adj.P.Val, decreasing=FALSE),]
  #colnames(stats)[1] = "gene"
  #tn = paste(pd,name,m.n,g1.n,"vs",g2.n,"diff_exp.txt", sep="_")
  tn = paste(pd,name,"mutation_impact_diff_exp.txt", sep="_")
  write.table(stats, file=tn, quote=F, sep = '\t', col.names=NA)
  
  return(stats)
} 

#mut = OV_mut_g; exp = OV_PNNL_Pho_c; name="OV PNNL Phosphoproteome" # for testing
### find_diff_exp : find differentially expressed proteins given specific mutations; calls diff_exp ###
# expect a mutation matrix and an expression matrix
find_diff_exp = function(mut, exp, name="data", fdr_cutoff=0.1){
  mut = mut[,colnames(mut) %in% colnames(exp)]
  #   fold_change = matrix(,nrow=nrow(exp),ncol=0)
  #   t_fdr = matrix(,nrow=nrow(exp),ncol=0)
  fold_change = vector("list")
  t_fdr = vector("list")
  geneList2 = c()
  for (gene in row.names(mut)){
    #for (gene in geneList){
    #if (sum(is.na(mut[gene,])) + sum(mut[gene,]=="silent", na.rm=T) + sum(mut[gene,]=="wt", na.rm=T) < ncol(mut)*0.98){
    geneList2 = c(geneList2,gene)
    wt = as.vector(c(colnames(mut[,mut[gene,] == "silent"]),colnames(mut[,mut[gene,] == "wt"]),colnames(mut[,mut[gene,] == "intronic"])))
    mutant= as.vector(colnames(mut[, ! colnames(mut) %in% wt]))
    
    diff_x = diff_exp(m = exp, g1 = mutant, g2 = wt, name=paste(gene,"mut",name,sep="_"))
    #m = exp; g1 = mutant; g2 = wt
    #diff_x = diff_exp(m = exp, g1 = mutant, g2 = wt, g1.n = gene, g2.n = "wt")#,plot=T,pathwayC=T, pathwayT=TRUE)
    
    gene_diff = diff_x[,"meanChange",drop=F]
    colnames(gene_diff) = gene
    fold_change[[gene]] = gene_diff
    #fold_change = cbind(fold_change, gene_diff)
    
    gene_fdr = diff_x[,"t_test_fdr",drop=F]
    colnames(gene_fdr) = gene
    t_fdr[[gene]] = gene_fdr
    #t_fdr = cbind(t_fdr, gene_fdr)
    
    # extract only the significant markers
    #     markers = row.names(gene_fdr)[as.numeric(gene_fdr[,1]) <= fdr_cutoff]
    #     if (length(markers)==0 || is.na(markers[1])){next}
    #     mut_t = t(mut[gene,])
    #     exp_t = t(exp[markers,])
    #     mut_exp = merge(mut_t, exp_t, by="row.names")
    #     mut_exp$Status = mut_exp[,gene]!="wt"
    #     mut_exp = mut_exp[,-c(1,2)]
    #plot_diff_exp_violin(mut_exp, gene, name=name) # don't plot the violin plots for now
    #}
    #}
  }
  fold_change_m = do.call(cbind,fold_change)
  t_fdr_m = do.call(cbind,t_fdr)
  
  plot_diff_exp_heatmap(fold_change_m, t_fdr_m, name=name)
  return(list("fold_change"=fold_change, "t_fdr"=t_fdr))
}

### plot_diff_exp : plot heatmap of differentially expressed genes from find_diff_exp ###
# called from find_diff_exp
plot_diff_exp_heatmap = function(fold_change, t_fdr, fdr_cutoff=0.1, name="data"){
  w = ncol(fold_change) + 1.5
  # plot heatmap 
  fold_change_df = as.data.frame(fold_change)
  t_fdr_df = as.data.frame(t_fdr)
  fold_change_df$gene = rownames(fold_change_df)
  t_fdr_df$gene = rownames(t_fdr_df)
  fold_change.m = melt(fold_change_df, id.var="gene")
  t_fdr.m = melt(t_fdr_df, id.var="gene")
  t_fdr.m$sig = as.numeric(as.character(t_fdr.m$value)) <= fdr_cutoff
  
  fn = paste(pd, name, 'diff_exp.pdf',sep ="_")
  fdr.colors=c("NA", "#000000")
  min_d = min(as.numeric(as.character(fold_change.m$value)), na.rm=T)
  max_d = max(as.numeric(as.character(fold_change.m$value)), na.rm=T)
  bound = max(c(max_d, -min_d))
  
  p = ggplot()
  p = p + geom_tile(data=fold_change.m, aes(x=variable, y=gene, fill=as.numeric(value)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-bound,bound))
  p = p + geom_tile(data=t_fdr.m, aes(x=variable, y=gene, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <=",fdr_cutoff),values = fdr.colors)
  p = p + labs(title = name, x="Mutated Gene", y="Expression") + theme_bw() + 
    theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
  p = p + coord_fixed()
  p
  ggsave(file=fn, width=w, limitsize=FALSE, height=15, useDingbats=FALSE)
}

### plot_diff_exp : plot heatmap of differentially expressed genes from find_diff_exp ###
# called from find_diff_exp
plot_diff_exp_violin = function(mut_exp, gene,name="data"){ ### select only the key mutations and the differentially expressed genes
  # change height based on the number of markers
  h = ncol(mut_exp)*1.5
  mut_exp.m = melt(mut_exp, id.var="Status")
  # plot violin plots faceted by marker genes
  fn = paste(pd, name, gene, "mutational_impact_violin.pdf", sep="_")
  p = ggplot(data=mut_exp.m)
  p = p + facet_grid(variable~.)
  p = p + geom_violin(aes(x=Status, y=value, fill=Status),alpha=0.5) + guides(fill=FALSE) 
  p = p + geom_jitter(aes(x=Status, y=value)) #+ geom_point(aes(x=Status, y=value)) 
  p = p + labs(title = name, x = paste(gene,"Mutation Status"), y = "Expression") + theme_bw()
  p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8, angle = 90))
  p
  ggsave(file=fn, height=h, width=3, limitsize=FALSE, useDingbats=FALSE)
}



