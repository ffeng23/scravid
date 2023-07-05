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

#start get the running error models
#we want to do fitting for each different cluster separately.
cluster<-1
for( i in unique(cvid.combined@meta.data$seurat_clusters))
{
cdata<-subset(cvid.combined, subset= seurat_clusters==i)

cat("Processing data.........\n")
mcd<-GetAssayData(cdata, slot="counts")
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

cat("Running error model estimating cluster ", i," cells........\n")

#start fitting for error model!!!
#newly updated 5/31/2023
#

system.time(
knn <- knn.error.models(cd, k = ncol(cd)/3, n.cores = 5, 
min.count.threshold = 1,
min.size.entries =2000, 
min.nonfailed = 5, 
save.model.plots=T,max.model.plots = 2, verbose=4)
)
 filename<-paste0("knn.error.models_4th_cluster_",i,".Rds")
 saveRDS(file=here("Output",filename),knn)
}
cat("Done.........!!\n")
rm(cdata)
rm(cd)
gc() 

cat("Done.........!!\n")


