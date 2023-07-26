#R code to do pathway analysis using pagoda (scde)
#https://davetang.org/muse/2013/12/16/bioconductor-annotation-packages/
#	in this one we do cell analysis 
#------------6/7/2023
# start run pagoda pathway analysis following the
# run of knn.error.model in pathway_pagoda_v1.0.R
#
#============================
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
mcd<-GetAssayData(cvid.combined, slot="counts")
mcd<-as.matrix(mcd)

#x<-apply(mcd, 2, function(x){
#		sum(x>0)
#	})
mcd<-apply(mcd,2,
	function(x) {storage.mode(x) <- 'integer'; x})
cat("Clean up the data counts......\n")
# get data counts matrix
cd <- clean.counts(mcd)
# check the final dimensions of the read count matrix
dim(cd)

rm(mcd)
gc()


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

#skip above error model estimation and loaded the one saved from bumc server theano
knn<-readRDS(file=here("Output","knn.error.models_1st_allCells.Rds"))


#start doing the downstream analysis
# variance normalization
cat("Start doing varnorm.......\n")
#system.time(
#varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, 
#		n.cores = 4, plot = TRUE)
#)
# list top overdispersed genes
#sort(varinfo$arv, decreasing = TRUE)[1:10]

#varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))
#saveRDS(varinfo, file="varinfo_allCells_1st.Rds")
varinfo<-readRDS(file=here("Output","varinfo_allCells_1st.Rds"))

#doing pca 
# translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids 
# convert GO lists from ids to gene names
gos.interest <- unique(ls(org.Hs.egGO2ALLEGS))#[1:100],"GO:0022008","GO:0048699", "GO:0000280", "GO:0007067")) 
go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment

pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 19)

df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = F, 
	z.score = 1.96)

saveRDS(pwpca, file=here("Output","Output_allCells_1st_pwpca.Rds" ))
#pwpca<-readRDS(file=here("Output","Output_allCells_1st_pwpca.Rds" ))

#gene cluster
clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), 
		n.clusters = 100, n.cores = 6, plot = F)

saveRDS(clpca, file=here("Output","Output_allCells_1st_clpca.Rds" ))
