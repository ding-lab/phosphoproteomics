### convert_uniprot2refseq.R ### 
# Yige Wu @ WashU 2016 Nov
# convert the uniprot id in k_s_table to refseq id in pho_data

# install biomaRt
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library("biomaRt")

listDatasets(useMart("refseq"))
listFilters(ensembl)
listAttributes(ensembl)

ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")

values<- as.vector(head(k_s_table$SUB_ACC_ID))

getBM(attributes=c("refseq_peptide",), filters = "refseq_mrna", values = values, mart= ensembl)

mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
refseq <- c("NM_006945", "NM_152486", "NM_198317")
getBM(filters="refseq_dna", attributes="external_gene_id", values=refseq, mart=mart)

mart = useMart('ensembl')
listDatasets(mart)


