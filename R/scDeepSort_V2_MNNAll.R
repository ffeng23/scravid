#R code to run python scDeepSort

#

library(reticulate)
library(here)
library(Seurat)
library(dplyr)

###show configuration
py_config()

##load module and run demo
ds<-import("deepsort")
ds$demo()

############
#python code
#
# import deepsort
# deepsort.demo()
#

##########################
#see how can we do real analysis
##process data and save data and make 
## them ready for analysis
data.dir<-"Data"
geneinfo<-readRDS(here(data.dir,"geneinfo.rds"))

geneinfo <- geneinfo %>% filter(species=="Human")



#cluster 2, scRNASeq data.
CVIDagg6.int<-readRDS(here("Output","CVIDagg6_MNN_2K.rds"))
		#loaded the CVIDagg6.Low.Cluster2


cvid.c2<- CVIDagg6.int[['RNA']]@data

# revising gene symbols
genename<- rownames(cvid.c2)
genename1<- genename[genename %in% geneinfo$Symbol]
genename2<- genename[!genename %in% geneinfo$Symbol]
genename3<- genename2[genename2 %in% geneinfo$Synonyms]
genename4<- rep('NA',length(genename3))
for (i in 1:length(genename3)) {
  d1<- geneinfo[geneinfo$Synonyms == genename3[i],]$Symbol
  if(length(d1) == 1){
    genename4[i]<- d1
  }
}
genename3<- c(genename1,genename3)
genename4<- c(genename1,genename4)
genedata<- data.frame(raw_name = genename3,new_name = genename4,stringsAsFactors = F)
genedata<- genedata[!genedata$new_name == 'NA',]
genedata1<- as.data.frame(table(genedata$new_name),stringsAsFactors = F)
genedata1<- genedata1[genedata1$Freq == 1,]
genedata<- genedata[genedata$new_name %in% genedata1$Var1,]
cvid.c2<- cvid.c2[genedata$raw_name,]
all(rownames(cvid.c2) == genedata$raw_name)
rownames(cvid.c2)<- genedata$new_name
all(rownames(cvid.c2) == genedata$new_name)
all(rownames(cvid.c2) %in% geneinfo$Symbol)

cvid.c2<- as.matrix(cvid.c2)
write.csv(cvid.c2,file = here(data.dir,'cvid.c2_data.csv'))

#start calling 

model <- ds$DeepSortPredictor(species='human',
                          tissue='Blood')
test_files = here(data.dir,'cvid.c2_data.csv')
model$predict(test_files, save_path=NULL)->results
table(results$cell_type)
#do this to the other 2 clusters

#also Todo??? train own model??? using celldex data for 
# fine label??? 