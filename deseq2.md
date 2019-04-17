## Analysis of changes in taxonomic abundance
```
library(plyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library("DESeq2")

otu_table<-read.delim("OTU_table.txt", header=T, row.names=1)
meta<-read.delim("map.txt", header=T)

otu_table[is.na(otu_table)] <- 0

meta$Description<-relevel(meta$Description, ref='control')


bee_deseq<-DESeqDataSetFromMatrix(countData = otu_table, colData=meta, design=~Description)

bee_dds<-DESeq(bee_deseq)

results_de<-results(bee_dds)
results_table<-as.data.frame(results_de)
results_table$OTU<-row.names(results_table)

#add taxonomy to data
annotations<-read.delim("otus_sintax_annotations.txt", header=F)
annotations<-annotations[,c(1,4)]

results_table<-merge(results_table, annotations, by.x='OTU', by.y='V1')


significant_taxa<-which(results_table$pvalue<0.05)

sig_results<-results_table[significant_taxa,]

write.table(sig_results, "significant_deseq_results.txt", sep="\t", row.names=F, quote=F)
```
