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
cvid.combined<-readRDS(file=here("Output","CVIDagg6.int_2K_singleR.Rds"))
DefaultAssay(cvid.combined)<-"RNA"

cat("Processing data.........\n")
mcd<-GetAssayData(cvid.combined, slot="counts")
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

cat("Running error model estimating........\n")
#modify and clean up for cell names
#x <- gsub("^Hi_(.*)_.*", "\\1", colnames(cd))
#l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(factor(x, levels = c("NPC", "GW16", "GW21", "GW21+3")))]

#start fitting for error model!!!
#
system.time(
knn <- knn.error.models(cd, k = ncol(cd)/80, n.cores = 40, 
		min.count.threshold = 1,
		min.size.entries =2000, 
		min.nonfailed = 50, 
		save.model.plots=T,max.model.plots = 2, verbose=4)
)

saveRDS(file=here("Output","knn.error.models.Rds"),knn)

cat("Done.........!!\n")