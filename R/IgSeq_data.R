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
