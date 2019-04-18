## Analysis of OTUs that respond to diet (relative to control)
```
library(plyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library("DESeq2")

#read in OTU table
otu_table<-read.delim("bee_microbiome/OTU_table.txt", row.names = 1, header=T)
dim(otu_table)
#614 by 78

#read in taxonomic annotations
annotations<-read.delim("bee_microbiome/otus_sintax_annotations.txt", header=F)
row.names(annotations)<-annotations$V1

#split tax annotations
annotations2<-separate(annotations, V4, sep = ",", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#ID chloroplast and mitochondria OTUs
chloroplast_otus<-which(annotations2$Order=='o:Chloroplast' | annotations2$Family=="f:Mitochondria")
length(chloroplast_otus)
#105

chloro_annotations<-annotations2[-chloroplast_otus,]

#remove OTUs from table
otu_table<-otu_table[row.names(otu_table) %in% row.names(chloro_annotations),]
dim(otu_table)
#533 by 78, 81 OTUs lost

#read in metadata
meta<-read.delim("bee_microbiome/map.txt", header=T)

#remove NAs
otu_table[is.na(otu_table)] <- 0

#fix diets
meta$Description<-relevel(meta$Description, ref='control')

#make a deseq object
bee_deseq<-DESeqDataSetFromMatrix(countData = otu_table, colData=meta, design=~Description)

#run DESeq
bee_dds<-DESeq(bee_deseq)

results_de<-results(bee_dds)
results_table<-as.data.frame(results_de)
results_table$OTU<-row.names(results_table)

#add taxonomy annotations to data
annotations<-read.delim("bee_microbiome/otus_sintax_annotations.txt", header=F)
annotations<-annotations[,c(1,4)]
results_table<-merge(results_table, annotations, by.x='OTU', by.y='V1')

#make a table of significant data
significant_taxa<-which(results_table$pvalue<0.05)
sig_results<-results_table[significant_taxa,]

write.table(sig_results, "bee_microbiome/significant_deseq_results.txt", sep="\t", row.names=F, quote=F)
```
