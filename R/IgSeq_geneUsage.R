#R code to do analysis on Gene usage, etc

#Read data from IgSeq_data.R

#
library(dplyr)
library(tidyverse)
library(readr)
library(here)
library(compositions)
library(ggfortify)

#read the saved data
data.dir<-"Data" 
load(here(data.dir,"BCR_vdj.RData"))
	#two variable/data frame read in
	#igData, colData


#let's first do heavy chain
IgHVs<- igData %>%
	filter(chain=="IGH", high_confidence==TRUE) 

#now only do IgHVs for now
IgHVs <- IgHVs %>%
	select( barcode, contig_id, v_gene, subs, cells, disease)

#put it into a wide format
IgHVs.long<- IgHVs %>% 
	group_by(v_gene, subs, cells, disease) %>%
	summarize(freq=n())
#k<-aggregate(IgHVs$barcode, by=list(IgHVs$v_gene, IgHVs$subs, IgHVs$cells), FUN=length) %>% View("f2")
#k[order(k$Group.1),] %>% View("f3")

IgHVs.wide<- IgHVs.long %>%
	pivot_wider(names_from=v_gene, values_from=freq, values_fill=0)


#now let's got them into compositions
#clean up first to get rid of VGenes too few
IgHVs.wide.clean<- IgHVs.wide[,
	-(which(apply(IgHVs.wide[,-c(1,2,3)], 2, sum)<18))]

#now let's get the detection limit 
dl<-apply(IgHVs.wide.clean[,-c(1,2,3)],1, sum)
dl<-1/dl
IgHVs.wide.clean <- IgHVs.wide.clean [,-c(1,2,3)]
IgHVs.wide.clean<-as.data.frame(IgHVs.wide.clean)
rownames(IgHVs.wide.clean)<-paste0(IgHVs.wide$subs,"-", IgHVs.wide$cells)

dl<-rep(dl,dim(IgHVs.wide.clean)[2])
dl<-matrix(dl, ncol=dim(IgHVs.wide.clean)[2], nrow=dim(IgHVs.wide.clean)[1], byrow=F)
IgHVs.comp<-acomp(IgHVs.wide.clean)
IgHVs.comp<-zeroreplace(IgHVs.comp, d=dl, a=1/3)

IgHVs.clr<-clr(IgHVs.comp)
#usage.VGene.ilr<-ilr(usage.VGene.comp)
dt.clr<-data.frame(IgHVs.clr)
IgHVs.dat<-cbind(IgHVs.clr, IgHVs.wide[,c("subs","cells","disease")])
pca.clrd<-prcomp(dt.clr, scale=T)

g<-autoplot(pca.clrd,x=1,y=2,data=IgHVs.dat, 
		colour="subs",shape="cells",
				label=T , label.repel=T, frame=F,  size=4,#frame.type="t"
			loadings=F, loadings.label=F, 
			loadings.label.size=4, loadings.colour="orange", 
			loadings.label.colour="orange"
		)

g34<-autoplot(pca.clrd,x=3,y=4,data=IgHVs.dat, 
		colour="subs",shape="cells",
				label=T, label.repel=T ,frame=F,  size=4,#frame.type="t"
			loadings=F, loadings.label=F, 
			loadings.label.size=4, loadings.colour="orange", 
			loadings.label.colour="orange"
		)

g56<-autoplot(pca.clrd,x=5,y=6,data=IgHVs.dat, 
		colour="subs",shape="cells",
				label=T , label.repel=T,frame=F,  size=4,#frame.type="t"
			loadings=F, loadings.label=F, 
			loadings.label.size=4, loadings.colour="orange", 
			loadings.label.colour="orange"
		)

g78<-autoplot(pca.clrd,x=7,y=8,data=IgHVs.dat, 
		colour="subs",shape="cells",
				label=T , label.repel=T,frame=F,  size=4,#frame.type="t"
			loadings=F, loadings.label=F, 
			loadings.label.size=4, loadings.colour="orange", 
			loadings.label.colour="orange"
		)

####now let's do anova 
	library(car)
	options(contrasts = c("contr.sum", "contr.poly"))
#now we need to do individual ANOVAs on each PCs 
pc.uplim<-8
i<-1
model.ind.anova<-NULL

for(i in 1:pc.uplim)
{
	model.ind<-lm(pca.clrd$x[,i]~IgHVs.dat$disease*IgHVs.dat$cells);
	temp<-cbind(Anova(model.ind, type=2)[c(1:3),],PCs=i, effect=rownames(Anova(model.ind, type=2))[1:3])
    #temp<-cbind(Anova(model.ind, type=3)[c(1:8),],PCs=i, effect=rownames(Anova(model.ind, type=3))[1:8])
	rownames(temp)<-NULL
	temp$effect<-gsub("IgHVs.dat\\$","",temp$effect)
	temp$p.adj.ind<-p.adjust(temp$"Pr(>F)", method="BH")
	model.ind.anova<-rbind(model.ind.anova,temp)
	
}	

#now let's do the individual tests.
IgHVs.dat<-as.data.frame(IgHVs.dat)
model.IGHV.anova<-data.frame()
for(i in 1:(dim(IgHVs.dat)[2]-3))
{
	model.IGHV<-lm(IgHVs.dat[,i]~IgHVs.dat$disease*IgHVs.dat$cells)
	temp<-cbind(Anova(model.IGHV, type=2)[c(1:3),],IgHV=colnames(IgHVs.dat)[i], effect=rownames(Anova(model.ind, type=2))[1:3])
    #temp<-cbind(Anova(model.ind, type=3)[c(1:8),],PCs=i, effect=rownames(Anova(model.ind, type=3))[1:8])
	rownames(temp)<-NULL
	temp$effect<-gsub("IgHVs.dat\\$","",temp$effect)
	temp$p.adj.ind<-p.adjust(temp$"Pr(>F)", method="BH")
	model.IGHV.anova<-rbind(model.IGHV.anova,temp)
}