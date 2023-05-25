#R code to do pathway analysis using pagoda (scde)
#
#--- started 5/23/2023


library(scde)
library(here)
library(dplyr)
library(org.Hs.eg.db)
library(Seurat)


#start loading data
# for scde/pagoda we need raw count (integer counts)

data.dir<-"Data"

cvid.combined<-readRDS(file=here("Output","CVIDagg6_MNN_2K.rds"))
DefaultAssay(cvid.combined)<-"RNA"
mcd<-GetAssayData(cvid.combined, slot="counts")
mcd<-as.matrix(mcd)

x<-apply(mcd, 2, function(x){
		sum(x>0)
	})
mcd<-apply(mcd,2,
	function(x) {storage.mode(x) <- 'integer'; x})

# get data counts matrix
cd <- clean.counts(mcd)
# check the final dimensions of the read count matrix
dim(cd)

rm(mcd)
gc()

#modify and clean up for cell names
#x <- gsub("^Hi_(.*)_.*", "\\1", colnames(cd))
#l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(factor(x, levels = c("NPC", "GW16", "GW21", "GW21+3")))]

#start fitting for error model!!!
#
knn <- knn.error.models(cd, k = ncol(cd)/30, n.cores = 3, min.count.threshold = 2, min.nonfailed = 5, 
		save.model.plots=T,max.model.plots = 2, verbose=4)
