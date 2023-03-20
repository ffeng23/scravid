#R code to process the clonotype and do diversity
# analysis

#--- started 3/12/2023 Feng

#the data were saved at IgSeq_data.R


library(here)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(SpadeR)
library(iNEXT)


#load data
data.dir<-"Data"

load(here(data.dir,"clonotypeEnclone.RData"))

#now let's clean up

clones<-cloneData %>% select(group_id, group_ncells, clonotype_id,
	exact_subclonotype_id, cid, subs, disease )
clones.group<-clones %>% group_by(group_id,subs,disease) %>%
	summarize(n_cells=unique(group_ncells))
#turn it into wider format

clones.wide <- clones.group %>% select(-disease) %>%
	pivot_wider(names_from=subs, values_from=n_cells, 
		values_fill=0)

#
d.S1<-ChaoSpecies(clones.wide$S1,"abundance", k=3, conf = 0.95)

d.S2<-ChaoSpecies(clones.wide$S6,"abundance", k=3, conf = 0.95)

d<-iNEXT(x=as.matrix(clones.wide[,-1]), q=0, datatype="abundance", 
		size=NULL, endpoint=20000, knots=1000, se=TRUE, 
		conf=0.95, nboot=50) 

d4<-iNEXT(x=as.matrix(clones.wide[,-1]), q=c(0,1,2,3,4), datatype="abundance", 
		size=NULL, endpoint=20000, knots=1000, se=TRUE, 
		conf=0.95, nboot=50) 
#save(d, d4, file="diversity.RData")
load(file=here(data.dir,"diversity.RData"))
ggiNEXT(d4, type=2, facet.var="Order.q")


####estimator
#S1
d4.trans<-d4$AsyEst
index<-4
subj<-1
dataInfo<-d4$DataInfo
rownames(dataInfo)<-dataInfo$Assemblage
d4.trans<-cbind(d4.trans,dataInfo[d4.trans$Assemblage,"n"])
names(d4.trans)[dim(d4.trans)[2]]<-"n"

scale.index<-dim(d4.trans)[2]

for(i in 1:length(unique(d4.trans$Assemblage)))
{
	d4.trans[c(1:3)+3*(i-1),index]<-d4$AsyEst[c(1:3)+3*(i-1),index]/d4.trans[3*(i-1)+1,scale.index]
}


###now let's do stats
#start plotting iNext
#orignial
ggiNEXT(d4, type=1)#, facet.var="Order.q")
group.info<-as.data.frame(unique(clones.group[,c(2:3)]))
rownames(group.info)<-group.info$subs
d4_2<-d4$AsyEst
d4_2$newDiv<-factor(d4_2$Diversity, levels=c("Species richness",
		"Shannon diversity", "Simpson diversity"))
d4_2$newDiv<-as.integer(d4_2$newDiv)
d4_2$disease<-group.info[d4_2$Assemblage,"disease"]
ggplot(d4_2, aes(x=newDiv, y=Estimator,colour=disease))+
	geom_line(aes(linetype=Assemblage),linewidth=1.2)+
	geom_point(aes(shape=Assemblage),size=5)

#transformed
d4.trans_2<-d4.trans
d4.trans_2<-cbind(d4.trans_2, d4_2[,c("newDiv",
		"disease"
		)])

ggplot(d4.trans_2, aes(x=newDiv, y=Estimator,colour=disease))+
	geom_line(aes(linetype=Assemblage),linewidth=1.2)+
	geom_point(aes(shape=Assemblage),size=5)


###diversity using space for estimating high order.
library(SpadeR)
ds1<-Diversity(clones.wide$S1,"abundance", q=seq(1,4,0.5))

ds2<-Diversity(clones.wide$S2,"abundance", q=seq(1,3.5,0.5))

ds3<-Diversity(clones.wide$S3,"abundance", q=seq(1,3,0.5))

ds4<-Diversity(clones.wide$S4,"abundance", q=seq(1,4,0.5))

ds5<-Diversity(clones.wide$S5,"abundance", q=seq(1,3,0.5))

ds6<-Diversity(clones.wide$S6,"abundance", q=seq(1,3,0.5))

