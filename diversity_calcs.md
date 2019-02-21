## Diversity values
```
#Read in OTU table
otu_table<-read.delim("OTU_table.txt", header=T, sep='\t', row.names=1)

#read in meta data
meta<-read.delim('map.txt', header=T)
meta$Description<-gsub('06EAA', 'T06EAA', meta$Description)
meta$Description<-gsub('10EAA', 'T10EAA', meta$Description)

#split meta data file
meta_split<-split(meta, meta$Description)

#transpose table
otu_table_t<-t(otu_table)

#sequencing depth
rowSums(otu_table_t)

min(rowSums(otu_table_t))
[1] 13999

#rarefy data
otu_rare<-rrarefy(otu_table_t, sample=13999)
#otu_rare<-otu_rare[row.names(otu_rare) %in% meta_split$T06EAA$SampleID,]

#otu_rare<-rrarefy(otu_table_t, sample=13999)
#otu_rare<-otu_rare[row.names(otu_rare) %in% meta_split$T10EAA$SampleID,]


#calculate NMDS
bee.nmds<-metaMDS(otu_rare, distance='bray')

#calculate PCoA
bee.pcoa<-capscale(otu_rare~1, distance='bray')
bee.scores<-scores(bee.pcoa)
bee.coords<-as.data.frame(bee.scores$sites)
bee.coords$SampleID<-row.names(bee.coords)



#combine coords w/meta
bee.coords<-merge(bee.coords, meta, by='SampleID')

#calcualte variation explained
100*round(bee.pcoa$CA$eig[1]/sum(bee.pcoa$CA$eig), 3)
#29.5
100*round(bee.pcoa$CA$eig[2]/sum(bee.pcoa$CA$eig), 3)
#12

#plot it
ggplot(bee.coords, aes(MDS1, MDS2, colour=as.factor(Week)))+
geom_point(aes(size=2))+
theme_bw()+
scale_colour_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"))

bee.div<-diversity(otu_rare, index='shannon')
bee.rich<-rowSums(otu_rare>0)
s16.div<-as.data.frame(cbind(bee.div, bee.rich))

#add sample IDs
s16.div$SampleID<-row.names(s16.div)

#fix coloumn names
names(s16.div)<-c('Shannon', 'OTUs_Obs', 'SampleID')

#add metadata
s16.div<-merge(s16.div, meta, by='SampleID')

ggplot(s16.div, aes(as.factor(Week), Shannon, fill=Description))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Week)
  ```

