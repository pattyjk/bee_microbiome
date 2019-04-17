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

ggplot(s16.div, aes(as.factor(Week), OTUs_Obs, fill=Description))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Week)
```

## Test alpha div for significance
```
pairwise.t.test(s16.div$Shannon, s16.div$Description, p.adjust.method = 'fdr')
#       control T06EAApairwise.t.test(s16.div$Shannon, s16.div$Description, p.adjust.method = 'fdr')
#T06EAA 0.914   -     
#T10EAA 0.054   0.054 

pairwise.t.test(s16.div$OTUs_Obs, s16.div$Description, p.adjust.method = 'fdr')
#       control T06EAA
#T06EAA 0.814   -     
#T10EAA 0.025   0.025 

#ANCOVA
plot(s16.div$Week, s16.div$Shannon)
summary(aov(s16.div$OTUs_Obs ~ s16.div$Week + s16.div$Description))


summary(aov(s16.div$OTUs_Obs ~ as.factor(s16.div$Week) + s16.div$Description))
#                         Df Sum Sq Mean Sq F value  Pr(>F)   
#as.factor(s16.div$Week)  8   9284  1160.5   3.595 0.00158 **
#s16.div$Description      2   3134  1567.1   4.855 0.01074 * 
#Residuals               67  21629   322.8                   

ggplot(s16.div, aes(Week, OTUs_Obs))+
  geom_point()+
  geom_smooth(method = lm)+
  facet_wrap(~Description)

plot(s16.div$Week, s16.div$OTUs_Obs)

```

## Adonis
```
otu_table<-otu_table[,row.names(meta),drop=F]


meta<-read.delim('map.txt', header=T)
row.names(meta)<-meta$SampleID
meta<-meta[row.names(otu_rare),]

set.seed(441)
adonis(otu_rare ~ as.factor(Week) + Description, data = meta, permutations =10000, method = 'bray')

#adonis(formula = otu_rare ~ as.factor(Week) + Description, data = meta, permutations = 10000, method = "bray") 

#Permutation: free
#Number of permutations: 10000

#Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#as.factor(Week)  8    2.1013 0.26266  3.4598 0.26521 9.999e-05 ***
#Description      2    0.7355 0.36775  4.8440 0.09283 9.999e-05 ***
#Residuals       67    5.0865 0.07592         0.64197              
#Total           77    7.9233                 1.00000   


ggplot(bee.coords, aes(Week, MDS1))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~Description)

ggplot(bee.coords, aes(Description, MDS1))+
  geom_point()+
  geom_smooth()
```

