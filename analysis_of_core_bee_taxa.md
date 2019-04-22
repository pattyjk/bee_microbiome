## Core bee taxa analysis

Examining taxa considered to be core honey bee taxa.

Frischella perrara<br/>
  -Orbaceae<br/>
Bartonella apis<br/>
  -Bartonellaceae
Parasaccharibacter apium
 -Acetobacteraceae
Gluconobacter
  -Acetobacteraceae
Lactobacillus Firm-4
Lactobacillus Firm-5
Bifidobacterium asteroides
  -Bifidobacteriaceae
Snodgrassella alvi
  -Neisseriaceae
Gilliamella apicola-
 -Orbaceae

From: https://www.nature.com/articles/nrmicro.2016.43

```
library(vegan)
library(tidyr)
library(reshape2)
library(ggplot2)
library(plyr)

#read in taxonomic annotations
annotations<-read.delim("bee_microbiome/otus_sintax_annotations.txt", header=F)

#split annotations
annotations<-separate(annotations, V4, sep=',', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#extract taxa of interest
core<-which(annotations$Genus == "g:Frischella" | annotations$Genus == "g:Bartonella" | annotations$Genus == "g:Lactobacillus" | annotations$Genus == "g:Gluconobacter" | annotations$Genus == "g:Parasaccharibacter" | annotations$Genus == "g:Bifidobacterium" | annotations$Genus == "g:Snodgrassella" | annotations$Genus == "g:Gilliamella")
core_annotations<-annotations[core,]

##only getting Lactobacillus, so likely not a long enough fragment to get confident calls- relax to family
core_families<-which(annotations$Genus == "g:Lactobacillus" | annotations$Family == "f:Orbaceae" | annotations$Family == "f:Neisseriaceae" | annotations$Family == "f:Bifidobacteriaceae" | annotations$Family == "f:Acetobacteraceae")
core_fam_annotations<-annotations[core_families, ]
row.names(core_fam_annotations)<-core_fam_annotations$V1

#count # OTUs
nrow(core_fam_annotations)
#63 OTUs that match the 'core'

#read in OTU table
otu_table<-read.delim("bee_microbiome/otu_table_no_mito_chloro.txt", row.names=1, header=T)

#rarefy OTU table
otu_table_rare<-rrarefy(t(otu_table), sample=min(rowSums(t(otu_table))))
otu_table_rare<-as.data.frame(t(otu_table_rare))

#make core otu table
core_otu_table<-otu_table_rare[row.names(core_fam_annotations),]
core_otu_table2<-core_otu_table

#transpose table
core_otu_table$OTU<-row.names(core_otu_table)

#read in metadata
meta<-read.delim("bee_microbiome/map.txt", header=T)

#melt mdata and add metadaa
core_m<-melt(core_otu_table)
core_m<-merge(core_m, meta, by.x='variable', by.y='SampleID')

#ggplot(core_m, aes(Week, value, colour=OTU))+
#  theme_bw()+
#  geom_point()+
#  facet_grid(OTU~Description)+
#  geom_line()
```

## Calculate Kruskal-Wallis test
```
#split OTU table by OTU 
otu_split<-split(core_m, core_m$OTU)

#calculate kruskal wallis by diet
core_kruskal<-lapply(otu_split, function(x) kruskal.test(x$value, x$Description))

#convert lists to tables
core_kruskal2<-lapply(core_kruskal, function(x) as.data.frame(t(unlist(x))))

#bind data together
core_kruskal3<-ldply(core_kruskal2, rbind)

#write to file and read back in to fix funkyness
write.table(core_kruskal3, "bee_microbiome/core_kruskal3.txt", sep='\t', quote=F, row.names = F)
core_kruskal3<-read.delim("bee_microbiome/core_kruskal3.txt", header=T)

#correct p-values with FDR
core_kruskal3$fdr_p<-p.adjust(core_kruskal3$p.value, method = 'fdr')
```

## Calculate beta diversity of core taxa
```
#calculate a PCoA of core taxa
core_otu_table2_t<-as.data.frame(t(core_otu_table2))
core.pcoa<-capscale(core_otu_table2_t~1)

#extract coords to a dataframe
core.scores<-scores(core.pcoa)
core.coords<-as.data.frame(core.scores$sites)
core.coords$SampleID<-row.names(core.coords)

#add metadata
core.coords<-merge(core.coords, meta, by.x='SampleID', by.y='SampleID')

#calcualte variation explained
100*round(core.pcoa$CA$eig[1]/sum(core.pcoa$CA$eig), 3)
#39.3
100*round(core.pcoa$CA$eig[2]/sum(core.pcoa$CA$eig), 3)
#18

ggplot(core.coords, aes(MDS1, MDS2, colour=as.factor(Week)))+
  geom_point(aes(size=2))+
  theme_bw()+
  scale_colour_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"))
```

## Asess changes in community compositon of core taxa
```
#read in OTU table
meta<-read.delim('bee_microbiome/map.txt', header=T, row.names=1)
meta$Description<-gsub('06EAA', 'T06EAA', meta$Description)
meta$Description<-gsub('10EAA', 'T10EAA', meta$Description)

#reorder metadata to match OTU table
meta<-meta[row.names(core_otu_table2_t),]

#calculate PERMANOVA
set.seed(441)
adonis(core_otu_table2_t ~ as.factor(Week) + Description, data = meta, permutations =10000, method = 'bray')

#                 Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#as.factor(Week)  8    1.2059 0.150742  1.6884 0.15544 0.0044 ** 
#Description      2    0.5703 0.285153  3.1939 0.07351 0.0005 ***
#Residuals       67    5.9818 0.089281         0.77105           
#Total           77    7.7581                  1.00000           
```
