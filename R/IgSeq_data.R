#R code to read IgSeq/BCR data
# make them ready for other modules

#  -- Feng 03/04/2023
#
#

library(dplyr)
library(tidyverse)
library(readr)
library(here)


#read the data from disk
#assume the data saved in ./Data folder
data.dir<-"Data"
filename<-"filtered_contig_annotations.csv"

#subjects, 6 
subs<-c("S1","S2",
		"S3","S4",
		"S5","S6"	)
subs<-rep(subs, rep(3,6)); # each subject has 3 b cell populations

cells<-c("CD38Low","CD38Mid","CD38High")
cells<-rep(cells,6)

disease<-c("healthy","cvid","cvid","healthy","cvid","healthy")
disease<-rep(disease, rep(3,6))

paths<-c("PJ02/bcr_low_vdj/","PJ02/bcr_mid_vdj/","PJ02/bcr_high_vdj/",
	"PJ03/BCR/PJ03_low/outs/","PJ03/BCR/PJ03_mid/outs/","PJ03/BCR/PJ03_high/outs",
	"PJ04/BCR/PJ04_low/outs","PJ04/BCR/PJ04_mid/outs","PJ04/BCR/PJ04_high/outs",
	"PJ_072629/VDJ_0726/CD38Low/vdj_b/","PJ_072629/VDJ_0726/CD38Mid/vdj_b/","PJ_072629/VDJ_0726/CD38High/vdj_b/",
	"PJ_072629/VDJ_0729/CD38Low/vdj_b/","PJ_072629/VDJ_0729/CD38Mid/vdj_b/","PJ_072629/VDJ_0729/CD38High/vdj_b/",
	"PJ_1115/LowCD38/vdj_b/","PJ_1115/MidCD38/vdj_b/","PJ_1115/HighCD38/vdj_b/")

#loading subject 1. high/low/mid

colData<-data.frame()
igData<-data.frame()
#read all the files
for(i in 1:length(paths))
{
	#reading data set
	cat("reading subs ", subs[i]," B cells ", cells[i], 
			"; dir :", paths[i],"\n")
	x<-read_csv(here(data.dir, paths[i],filename))
	c<-data.frame(barcode=x$barcode, subs=subs[i],
			cells=cells[i],disease=disease[i])
	c$barcode<-paste0(c$barcode,"-",subs[i],"-",cells[i])
	colData<-rbind(colData,c)
	x$barcode<-paste0(x$barcode,"-",subs[i],"-",cells[i])
	x$subs=subs[i]
	x$cells=cells[i]
	x$disease=disease[i]
	igData<-rbind(igData,x)

}

#now ready to save the data
README<-"the Igseq/BCR data. The data were read in at IgSeq_data.R\n"
README<-paste0(README, "two data frames were saved, igData and colData\n")
README<-paste0(README, "igData is everything from filtered_contig_annotations.csv\n")
README<-paste0(README, "The data is the output from cellranger BCR vdj output\n")
README<-paste0(README,"The data are for all six subjects. colData is \n")
README<-paste0(README, "the metadata about sequences.\n")

save(file=here(data.dir,"BCR_vdj.RData"), igData, colData, README)


#now I start to read enclone clonotype data.
# the data were in the Data folder >> sub folder "enclone"
data.dir<-"Data/enclone"
filename<- rep("encloneOut.txt",6)

#subjects, 6 
subs<-c("S1","S2",
		"S3","S4",
		"S5","S6"	)
#subs<-rep(subs, rep(3,6)); # each subject has 3 b cell populations

#cells<-c("CD38Low","CD38Mid","CD38High")
#cells<-rep(cells,6)

disease<-c("healthy","cvid","cvid","healthy","cvid","healthy")
#disease<-rep(disease, rep(3,6))

paths<-c("PJ02",
	"PJ03","PJ04","PJ0726","PJ0729","PJ1115")

#loading subject 1. high/low/mid

rowData<-data.frame()
cloneData<-data.frame()
last.three<-c()
#read all the files
for(i in 1:length(paths))
{
	#reading data set
	cat("reading subs ", subs[i], 
			"; dir :", paths[i],"\n")
	#x<-read_csv(here(data.dir, paths[i],filename[i]))
	x<-read.table(file=here(data.dir, paths[i],filename[i]),
			sep=",", header=T, quote="\"'"
		)
	len<-length(names(x))
	last.three<-c(last.three,names(x)[c(len-2,len-1,len)])
	names(x)[c(len-2,len-1,len)]<-paste0("barcode",c(1:3))
	c<-data.frame(cid=paste0(x$group_id,"-",x$clonotype_id,"-",
			x$exact_subclonotype_id), subs=subs[i],
			cells=cells[i],disease=disease[i])
	c$cid<-paste0(c$cid,"-",subs[i])
	rowData<-rbind(rowData,c)
	x$cid<-c$cid #<-paste0(x$barcode,"-",subs[i])
	x$subs=subs[i]
	#x$cells=cells[i]
	x$disease=disease[i]
	cloneData<-rbind(cloneData,x)

}
#now, let's save them.
README<-"the clonotyping data in a big table. this \n"
README<-paste0(README, "is the output from enclone. \n")
README<-paste0(README, "I so far grouped/pooled all 3\n")
README<-paste0(README, "samples within each subject together\n")
README<-paste0(README, "to increase the chanace of bigger\n")
README<-paste0(README, "clones. I also will use the\n")
README<-paste0(README, "\"group\", which is bigger, to\n")
README<-paste0(README, "to get more clones. this is supposed\n")
README<-paste0(README, "to be different from clonotype,\n")
README<-paste0(README, "but I don't know what it is about\n")
README<-paste0(README, "it is said to be a functional similar\n")
README<-paste0(README, "rather than evolutionary\n")
README<-paste0(README, "check enclone for more (there is \n")
README<-paste0(README, "very little detail about it\n")
README<-paste0(README, "About the data last.three, it is the order\n")
README<-paste0(README, "of the last 3 fields of each enclone output, about the barcodes in 3 samples in our case.\n")
README<-paste0(README, "IN it, it follows the order of s1-s6 in trios\n")
save(file=here(data.dir,"clonotypeEnclone.RData"), cloneData, rowData, last.three, README)
