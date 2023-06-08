#R code to do pathway analysis using pagoda (scde)
#
#--- started 5/23/2023


library(scde)
library(here)
library(dplyr)
library(org.Hs.eg.db)
library(Seurat)
library(SeuratObject)


#start loading data
# for scde/pagoda we need raw count (integer counts)

data.dir<-"Data"

cat("Read data..........\n")

#cvid.combined<-readRDS(file=here("Output","CVIDagg6_MNN_2K.rds"))
cvid.combined<-readRDS(file=here("Output","CVIDagg6.int_2K_singleR.Rds"))
DefaultAssay(cvid.combined)<-"RNA"

#cvid<-subset(cvid.combined, subset= cvid=='CVID')
#cd38_high<-subset(cvid.combined, subset= cells=='CD38_high')

#cat("Processing data.........\n")
#mcd<-GetAssayData(cd38_high, slot="counts")
#mcd<-as.matrix(mcd)

#x<-apply(mcd, 2, function(x){
#		sum(x>0)
#	})
#mcd<-apply(mcd,2,
#	function(x) {storage.mode(x) <- 'integer'; x})
#cat("Clean up the data counts......\n")
# get data counts matrix
#cd <- clean.counts(mcd)
# check the final dimensions of the read count matrix
#dim(cd)

#rm(mcd)
#gc()


#cat("Running error model estimating cd38 high cells........\n")
##modify and clean up for cell names


##start fitting for error model!!!
##newly updated 5/31/2023
##
#system.time(
#knn <- knn.error.models(cd, k = ncol(cd)/15, n.cores = 10, 
#		min.count.threshold = 1,
#		min.size.entries =2000, 
#		min.nonfailed = 10, 
#		save.model.plots=T,max.model.plots = 2, verbose=4)
#)

#saveRDS(file=here("Output","knn.error.models_3rd_cd38high.Rds"),knn)

#cat("Done.........!!\n")

#rm(cd38_high)
#rm(cd)
#gc() 
##############################
#      HC  
##############################
cd38_mid<-subset(cvid.combined, subset= cells=='CD38_mid')

cat("Processing data.........\n")
mcd<-GetAssayData(cd38_mid, slot="counts")
mcd<-as.matrix(mcd)

x<-apply(mcd, 2, function(x){
		sum(x>0)
	})
mcd<-apply(mcd,2,
	function(x) {storage.mode(x) <- 'integer'; x})
cat("Clean up the data counts......\n")
# get data counts matrix
cd <- clean.counts(mcd)
# check the final dimensions of the read count matrix
dim(cd)

rm(mcd)
gc()

#cd<-cd[,1:200]
cat("Running error model estimating cd38 mid cells........\n")
#modify and clean up for cell names
#x <- gsub("^Hi_(.*)_.*", "\\1", colnames(cd))
#l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(factor(x, levels = c("NPC", "GW16", "GW21", "GW21+3")))]

#start fitting for error model!!!
#newly updated 5/31/2023
#
system.time(
knn <- knn.error.models(cd, k = ncol(cd)/15, n.cores = 10, 
		min.count.threshold = 1,
		min.size.entries =2000, 
		min.nonfailed = 10, 
		save.model.plots=T,max.model.plots = 2, verbose=4)
)

saveRDS(file=here("Output","knn.error.models_3rd_cd38mid.Rds"),knn)

cat("Done.........!!\n")
rm(cd38_mid)
rm(cd)
gc() 

#############################
#      HC  
##############################
cd38_low<-subset(cvid.combined, subset= cells=='CD38_low')

cat("Processing data.........\n")
mcd<-GetAssayData(cd38_low, slot="counts")
mcd<-as.matrix(mcd)

x<-apply(mcd, 2, function(x){
		sum(x>0)
	})
mcd<-apply(mcd,2,
	function(x) {storage.mode(x) <- 'integer'; x})
cat("Clean up the data counts......\n")
# get data counts matrix
cd <- clean.counts(mcd)
# check the final dimensions of the read count matrix
dim(cd)

rm(mcd)
gc()

#cd<-cd[,1:200]
cat("Running error model estimating cd38 low cells........\n")
#modify and clean up for cell names
#x <- gsub("^Hi_(.*)_.*", "\\1", colnames(cd))
#l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(factor(x, levels = c("NPC", "GW16", "GW21", "GW21+3")))]

#start fitting for error model!!!
#newly updated 5/31/2023
#
system.time(
knn <- knn.error.models(cd, k = ncol(cd)/15, n.cores = 10, 
		min.count.threshold = 1,
		min.size.entries =2000, 
		min.nonfailed = 10, 
		save.model.plots=T,max.model.plots = 2, verbose=4)
)

saveRDS(file=here("Output","knn.error.models_3rd_cd38low.Rds"),knn)

cat("Done.........!!\n")


