##### summarize_cluster.R #####
# Kuan-lin Huang @ WashU 2017 Jan

setwd("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_analysis/phospho_structure")
source("/Users/khuang/bin/LIB_exp.R")

k_s_table = read.table(header=T, stringsAsFactors = F, quote = "", fill=T, sep = "\t","/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt")
signal_gene = unique(c(k_s_table$GENE,k_s_table$SUB_GENE))
cancer_genes = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/cgenes_and_druggable.list', header=FALSE, stringsAsFactors = F)
c_genes = as.vector(t(cancer_genes))

cluster_f = "cptac.brca.phosphosites.mafish.cptac.brca.phosphosites.pairwise.singleprotein.collapsed.l0.p0.05.r10.clusters"
cluster = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, row.names=NULL, file = cluster_f)
cluster_site_count = data.frame(table(cluster$Cluster))
colnames(cluster_site_count) = c("Cluster_ID","Phosphosite_Count")

cluster_s_f = "cptac.brca.phosphosites.mafish.cptac.brca.phosphosites.pairwise.singleprotein.collapsed.l0.p0.05.r10.clusters.summary"
cluster_s = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, row.names=NULL, file = cluster_s_f)
colnames(cluster_s) <- c(colnames(cluster_s)[-1],"Other")
#colnames(cluster)[1]='clusterID'

cluster_s_c = merge(cluster_s,cluster_site_count,by="Cluster_ID")

cluster_s_c[order(cluster_s_c$Phosphosite_Count,decreasing = T),]
top_genes = cluster[order(cluster$Centrality,decreasing = T),]$gene

p = ggplot(data=cluster_s_c)
p = p + geom_point(aes(x=Centrality, y=Phosphosite_Count),alpha=0.5) + guides(fill=FALSE)
p = p + geom_text(aes(x=Centrality, y=Phosphosite_Count, label = ifelse(Centrality > 2 | Phosphosite_Count >=5,Centroid,NA)),vjust=-0.2,alpha=0.8,size=2)
p = p + theme_bw() + theme_nogrid()
p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
              axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
# p = p + coord_equal() + xlim(0,2) + ylim(0,2)
p
fn = paste(pd,"phospho-cluster_centrality.vs.siteCount.pdf",sep="_")
ggsave(file=fn, useDingbats=FALSE)

xmer_f = "cptac.brca.phospho.xmer"
xmer = read.table(header=F, quote = "", sep="\t", stringsAsFactors = F, row.names=NULL, file =xmer_f)
colnames(xmer) = c("Cluster","PDB_ID","Gene.Chain","nMutations","nResidues","TotalRecurrence","Mutations.Position")
xmer$Gene = gsub("\\|.*","",xmer$Gene.Chain)

p = ggplot(xmer,aes(x = nResidues))
p = p + geom_bar() + theme_bw() + theme_nogrid()
#p = p + labs(x = x_string, y="counts")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p
fn = paste(pd,"phosphosite_PDBstructurePresense_count.pdf",sep="_")
ggsave(file=fn, useDingbats=FALSE)

xmer_cgene = xmer[xmer$Gene %in% c_genes,]
p = ggplot(xmer_cgene,aes(x = nResidues, y = Gene))
p = p + geom_point(alpha=0.3) + theme_bw()
p = p + geom_text(aes(label = ifelse(nResidues > 2,Mutations.Position,NA)),vjust=-0.2,alpha=0.8,size=2)
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p
fn = paste(pd,"phosphosite_PDBstructurePresense_count_cgenes_sites.pdf",sep="_")
ggsave(file=fn, useDingbats=FALSE)

