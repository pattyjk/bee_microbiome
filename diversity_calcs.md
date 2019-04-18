## Beta Diversity analyses 
```
library(tidyr)
library(vegan)

#read in meta data
meta<-read.delim('bee_microbiome/map.txt', header=T)
meta$Description<-gsub('06EAA', 'T06EAA', meta$Description)
meta$Description<-gsub('10EAA', 'T10EAA', meta$Description)

#read in OTU table
otu_table<-read.delim("bee_microbiome/OTU_table.txt", row.names = 1, header=T)
dim(otu_table)
#614 by 78

#read in taxonomic annotations
annotations<-read.delim("bee_microbiome/otus_sintax_annotations.txt", header=F)
row.names(annotations)<-annotations$V1

#split tax annotations
annotations2<-separate(annotations, V4, sep = ",", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#ID chloroplast OTUs
chloroplast_otus<-which(annotations2$Order=='o:Chloroplast' | | annotations2$Family=="f:Mitochondria")
length(chloroplast_otus)
#105

chloro_annotations<-annotations2[-chloroplast_otus,]

#remove OTUs from table
otu_table<-otu_table[row.names(otu_table) %in% row.names(chloro_annotations),]

#transpose and rarefy table
otu_table_t<-t(otu_table)
min(rowSums(otu_table_t))
#10307
otu_rare<-rrarefy(otu_table_t, sample=10307)

#calculate PCoA
bee.pcoa<-capscale(otu_rare~1, distance='bray')
bee.scores<-scores(bee.pcoa)
bee.coords<-as.data.frame(bee.scores$sites)
bee.coords$SampleID<-row.names(bee.coords)

#combine coords w/meta
bee.coords<-merge(bee.coords, meta, by.x='SampleID', by.y='X.SampleID')

#calcualte variation explained
100*round(bee.pcoa$CA$eig[1]/sum(bee.pcoa$CA$eig), 3)
#23
100*round(bee.pcoa$CA$eig[2]/sum(bee.pcoa$CA$eig), 3)
#15.6

#plot it
ggplot(bee.coords, aes(MDS1, MDS2, colour=as.factor(Week)))+
  geom_point(aes(size=2))+
  theme_bw()+
  scale_colour_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"))
```

## Calculate alpha diversity
```
bee.div<-diversity(otu_rare, index='shannon')
bee.rich<-rowSums(otu_rare>0)
s16.div<-as.data.frame(cbind(bee.div, bee.rich))

#add sample IDs
s16.div$SampleID<-row.names(s16.div)

#fix coloumn names
names(s16.div)<-c('Shannon', 'OTUs_Obs', 'SampleID')

#add metadata
s16.div<-merge(s16.div, meta, by.x='SampleID', by.y='X.SampleID')

ggplot(s16.div, aes(as.factor(Week), Shannon, fill=Description))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Week)

ggplot(s16.div, aes(as.factor(Week), OTUs_Obs, fill=Description))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Week)
```

## Test alpha div for significance
```
pairwise.t.test(s16.div$Shannon, s16.div$Description, p.adjust.method = 'fdr')

#        control T06EAA
#T06EAA 0.462   -     
#T10EAA 0.123   0.046

pairwise.t.test(s16.div$OTUs_Obs, s16.div$Description, p.adjust.method = 'fdr')
#       control T06EAA
#T06EAA 0.766   -     
#T10EAA 0.015   0.016 


#ANCOVA
plot(s16.div$Week, s16.div$Shannon)
summary(aov(s16.div$Shannon ~ s16.div$Week + s16.div$Description))
#                     Df Sum Sq Mean Sq F value Pr(>F)  
#s16.div$Week         1  0.013 0.01252   0.279 0.5991  
#s16.div$Description  2  0.281 0.14040   3.125 0.0498 *
#Residuals           74  3.324 0.04492         

plot(s16.div$Week, s16.div$OTUs_Obs)
summary(aov(s16.div$OTUs_Obs ~ as.factor(s16.div$Week) + s16.div$Description))
#                         Df Sum Sq Mean Sq F value  Pr(>F)   
#as.factor(s16.div$Week)  8   6186   773.3   2.934 0.00717 **
#s16.div$Description      2   3008  1504.2   5.708 0.00514 **
#Residuals               67  17656   263.5                  

ggplot(s16.div, aes(Week, OTUs_Obs))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~Description)
```

## PERMANOVA (adonis)
```
#reorder metadata to match OTU table
meta<-read.delim('bee_microbiome/map.txt', header=T)
meta$Description<-gsub('06EAA', 'T06EAA', meta$Description)
meta$Description<-gsub('10EAA', 'T10EAA', meta$Description)
row.names(meta)<-meta$X.SampleID
meta<-meta[row.names(otu_rare),]

set.seed(441)
adonis(otu_rare ~ as.factor(Week) + Description, data = meta, permutations =10000, method = 'bray')

#                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#as.factor(Week)  8    1.9254 0.24067  3.2251 0.25473 9.999e-05 ***
#Description      2    0.6332 0.31661  4.2427 0.08378 9.999e-05 ***
#Residuals       67    4.9998 0.07462         0.66149              
#Total           77    7.5584                 1.00000    

ggplot(bee.coords, aes(Week, MDS1))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~Description)

```





