## Analyze taxonomy
```
library(plyr)
library(reshape2)
library(ggplot2)
library(tidyr)


otu_table<-read.delim("OTU_table.txt", header=T)

annotations<-read.delim("otus_sintax_annotations.txt", header=F)
annotations<-annotations[,c(1,4)]

otu_table<-merge(otu_table, annotations, by.x='X.OTU.ID', by.y='V1')
otu_table<-otu_table[,-1]


tax<-otu_table$V4
otu_table<-otu_table[,-79]

#convert it to relative abundance
otu_table<-sweep(otu_table, 2, colSums(otu_table), '/')

#add tax back in
otu_table<-cbind(otu_table, tax)

#reshape OTU table
otu_m<-melt(otu_table)

#read in meta data
meta<-read.delim("map.txt", header=T)

#add meta data to otu table
otu_m<-merge(otu_m, meta, by.x='variable', by.y='SampleID')

#sep tax
otu_m<-separate(otu_m, tax, sep=',', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#remove NAs
otu_m$Class<-gsub("NA", "C: Unlcassified", otu_m$Class)


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
  
#summarize data
otu_sum<-ddply(otu_m, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Week", "Description"), summarise, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)


otu_sum<-ddply(otu_m, c( "Class", "Description"), summarise, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)

ggplot(otu_sum, aes(Description, Class, fill=mean))+
  geom_tile()+
  theme_bw()

ggplot(otu_sum, aes(Description, mean, fill=Class))+
  geom_bar(stat='identity')+ 
  theme_bw()
```
