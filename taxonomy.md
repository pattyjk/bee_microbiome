## Analyze taxonomy
```
library(tidyr)
library(reshape2)
library(plyr)

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
chloroplast_otus<-which(annotations2$Order=='o:Chloroplast' | annotations2$Family=="f:Mitochondria")
length(chloroplast_otus)
#105

chloro_annotations<-annotations2[-chloroplast_otus,]

#remove OTUs from table
otu_table<-otu_table[row.names(otu_table) %in% row.names(chloro_annotations),]
dim(otu_table)
#533 by 78, 81 OTUs lost

#add taxonomic annotations to the OTU table
otu_table$OTU<-row.names(otu_table)
otu_table<-merge(otu_table, annotations[,c(1,4)], by.x='OTU', by.y='V1')

#remove OTU names
otu_table<-otu_table[,-1]

#remvoe taxonomy
tax<-otu_table$V4
otu_table<-otu_table[,-79]

#convert it to relative abundance
otu_table<-sweep(otu_table, 2, colSums(otu_table), '/')

#add tax back in
otu_table<-cbind(otu_table, tax)

#reshape OTU table
otu_m<-melt(otu_table)

#read in meta data
meta<-read.delim("bee_microbiome/map.txt", header=T)

#add meta data to otu table
otu_m<-merge(otu_m, meta, by.x='variable', by.y='X.SampleID')

#sep tax
otu_m<-separate(otu_m, tax, sep=',', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))


#remove NAs
otu_m[is.na(otu_m)]<- 0
otu_m$Class<-gsub("0", "Unclassfied", otu_m$Class)

#plots
pal<-c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")
ggplot(otu_m, aes(variable, value, fill=Class))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_manual(values=pal)+
  facet_wrap(~Description, scales='free')+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Relative Abundance")+
  xlab("")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), text = element_text(size=14))
```
