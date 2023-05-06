## Seurat Pipeline Analysis 
#------------------
# updated 4/11/2023
# We are here doing testing
# 1) try to find the best condition to do integration
#	batch correction
# The code were originally from me, but modified by Erik, and I continue
# to use this one. (see the code from CVIDagg6_SeuratPiepline.R)
#==========================
## scRNAseq analysis of 5 samples
## 3 Healthy: ZA03 (Sub1), ZA61 (Sub4), ZA63 (Sub6)
## 3 CVIDa: CC33D (Sub2), CC22D (Sub3), CC27B (Sub5)

#the data were combine and downloaded to local.

#library
library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(patchwork)
library(cowplot)
library(stringr)
library(ggpubr)
library(batchelor)
library(SeuratWrappers)

#CVIDagg6.data <- Read10X(data.dir = "Documents/PhD Thesis Work/Data/Single Cell/CVIDagg6_05Dec2022/agg6_filtered_feature_bc_matrix")
CVIDagg6.data <- Read10X(data.dir = 
	"/home/feng/Windows/windowsD/feng/LAB/MSI/maglione/scRNA_analysis1_20220908/aggr6sub/aggr2212_6sbu/outs/count/filtered_feature_bc_matrix")
CVIDagg6 <- CreateSeuratObject(counts = CVIDagg6.data, project = "CVIDagg6", min.cells = 3, min.features = 200)
ct<-GetAssayData(CVIDagg6, slot="counts")
dat<-GetAssayData(CVIDagg6, slot="data")

#- Add meta.data
mito.genes <- grep(pattern = "^MT-", x = rownames(x = dat), value = TRUE)
percent.mito <- Matrix::colSums(dat[mito.genes, ]) / Matrix::colSums(dat)
CVIDagg6@meta.data$percent.mito <- percent.mito
CVIDagg6@meta.data$sub <- plyr::mapvalues(x = sapply(str_split(rownames(CVIDagg6@meta.data), "[-]"), function(x) x[2]), from = 1:18, to = rep(seq(1,6),c(3,3,3,3,3,3)))
CVIDagg6@meta.data$cvid <- plyr::mapvalues(x = sapply(str_split(rownames(CVIDagg6@meta.data), "[-]"), function(x) x[2]), from = 1:18, to = rep(c("HC","CVID","CVID","HC","CVID","HC"),c(3,3,3,3,3,3)))
CVIDagg6@meta.data$cells <- plyr::mapvalues(x = sapply(str_split(rownames(CVIDagg6@meta.data), "[-]"), function(x) x[2]), from = 1:18, to = rep(c("CD38_high","CD38_mid","CD38_low","CD38_high","CD38_mid","CD38_low","CD38_high","CD38_mid","CD38_low","CD38_high","CD38_low","CD38_mid","CD38_high","CD38_low","CD38_mid","CD38_high","CD38_low","CD38_mid"),c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)))

CVIDagg6@meta.data %>% View("meta.data")

#- QC
CVIDagg6 <- subset(CVIDagg6, subset = nFeature_RNA > 200 & nCount_RNA < 40000 & percent.mito < .15)

CVIDagg6 
##An object of class Seurat 
##21188 features across 46878 samples within 1 assay 
##Active assay: RNA (21188 features, 0 variable features)

#normalize with default setting
CVIDagg6<-NormalizeData(CVIDagg6, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable features
CVIDagg6 <- FindVariableFeatures(object = CVIDagg6, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 20000)

CVIDagg6.orig<-CVIDagg6

#Idents(CVIDagg6.int) <- CVIDagg6.int@meta.data$seurat_clusters
#plot the data before integration/batch correction.
#do scale and PCA
#CVIDagg6.int<-CVIDagg6
CVIDagg6<-ScaleData(CVIDagg6, verbose =T, 
	vars.to.regress = c("nCount_RNA", "percent.mito"), 
	features=rownames(CVIDagg6))
CVIDagg6<-RunPCA(CVIDagg6, npcs = 100, verbose = T)
print(CVIDagg6[["pca"]], dims = 1:10, nfeatures = 5)
#do clustering
ElbowPlot(object = CVIDagg6, ndims = 50)
CVIDagg6<- FindNeighbors(CVIDagg6, reduction = "pca", dims = 1:20)
CVIDagg6<- FindClusters(object = CVIDagg6, #reduction.type = "pca", 
                            dims.use = 1:20, resolution = .6#, #k.param = 30, print.output = TRUE, save.SNN = TRUE
)

#do UMAP
CVIDagg6 <- RunUMAP(object = CVIDagg6, reduction.use = "pca", 
		dims= 1:20, min_dist = .75)

png(file="scRNACluster_beforeIntergration.png",width=1000, height=1000)
ggarrange(
 DimPlot(CVIDagg6, reduction = "umap", label=T)
,DimPlot(CVIDagg6, reduction = "umap", label=T, group.by="cvid")
,DimPlot(CVIDagg6, reduction = "umap", label=T, group.by="sub")
,DimPlot(CVIDagg6, reduction = "umap", label=T, group.by="cells")
,nrow=2, ncol=2

)
dev.off()

#save the seurat object since it takes too long to run scale
README<-"the sc seurate object before doing integration/batch correct.\n"
README<-paste0(README,"Saved in CVIDagg6_SeuratePipeline_FF_v1.0.R\n")
README<-paste0(README," the object is still named CVIDagg6\n")
save(file="scObjectSeurate.before.RData", CVIDagg6, README)
load(file="scObjectSeurate.before.RData")
####now start doing the subsetting and then integrating 
#   by subsetting. (Note: intergrating among different subject
#		and each subject combining different cell populations )


#start doing the subsetting and then integrating
#   CVIDagg6.orig  this one is normalized!!!
CVIDagg6.ZA03<-subset(CVIDagg6.orig, subset= sub==1)
CVIDagg6.CC33<-subset(CVIDagg6.orig, subset= sub==2)
CVIDagg6.CC22<-subset(CVIDagg6.orig, subset= sub==3)
CVIDagg6.ZA61<-subset(CVIDagg6.orig, subset= sub==4)
CVIDagg6.CC27<-subset(CVIDagg6.orig, subset= sub==5)
CVIDagg6.ZA63<-subset(CVIDagg6.orig, subset= sub==6)

features.names=c("TNFRSF13B","TNFRSF13C", "NFKB1","NFKB2","RELA","RELB", "BCL2", "BCL2L1"
                 ,"BCL2L2", "BCL2A1","MCL1","ICOSLG")
#now we have done the normalization and findVariableFeatures, do directly integration
ancs<-SelectIntegrationFeatures(list(CVIDagg6.ZA03,CVIDagg6.CC33,CVIDagg6.CC22,CVIDagg6.ZA61,CVIDagg6.CC27,CVIDagg6.ZA63), nfeatures=2500) #find the variable features for doing integration
ancs<-union(ancs, features.names) #add ones that we are interested.
CVIDagg6.anchors <- FindIntegrationAnchors(object.list = list(CVIDagg6.ZA03,CVIDagg6.CC33,CVIDagg6.CC22,CVIDagg6.ZA61,CVIDagg6.CC27,CVIDagg6.ZA63), dims = 1:20, anchor.features=ancs)
CVIDagg6.int<- IntegrateData(anchorset = CVIDagg6.anchors, dims = 1:20)

CVIDagg6.int@meta.data %>% View("message")

#do scale and PCA
CVIDagg6.int<-ScaleData(CVIDagg6.int, verbose =T, vars.to.regress = c("nCount_RNA", "percent.mito"), features=rownames(CVIDagg6.int))
CVIDagg6.int<-RunPCA(CVIDagg6.int, npcs = 100, verbose = T)
print(CVIDagg6.int[["pca"]], dims = 1:10, nfeatures = 5)
#do clustering
ElbowPlot(object = CVIDagg6.int, ndims = 50)
CVIDagg6.int<- FindNeighbors(CVIDagg6.int, reduction = "pca", dims = 1:20)
CVIDagg6.int<- FindClusters(object = CVIDagg6.int, #reduction.type = "pca", 
                            dims= 1:20, resolution = .6#, #k.param = 30, print.output = TRUE, save.SNN = TRUE
)

#do UMAP
CVIDagg6.int <- RunUMAP(object = CVIDagg6.int, reduction = "pca", dims= 1:20, min.dist = .75)
DimPlot(CVIDagg6.int, reduction = "umap", label=T)
DimPlot(CVIDagg6.int, reduction = "umap", label=T, group.by="sub")
DimPlot(CVIDagg6.int, reduction = "umap", label=T, group.by="cvid")
DimPlot(CVIDagg6.int, reduction = "umap", label=T, group.by="cells")

png(file="scRNACluster_IntergrateBySun.png",width=1000, height=1000)
DimPlot(CVIDagg6.int, reduction = "umap", label=T, group.by=c("cvid","cells","sub", "seurat_clusters"))
dev.off()
save(CVIDagg6.int, file="scObjectSeurate.integrateSub.RData")
load(file="scObjectSeurate.integrateSub.RData")
##############################
#	now integrate by cells ###
##############################
CVIDagg6<-CVIDagg6.orig
#subsetting cells by subject ID
ZA03.high<-subset(CVIDagg6, subset= sub==1&cells=="CD38_high" )
ZA03.mid<-subset(CVIDagg6, subset= sub==1&cells=="CD38_mid" )
ZA03.low<-subset(CVIDagg6, subset= sub==1&cells=="CD38_low" )
CC33.high<-subset(CVIDagg6, subset= sub==2&cells=="CD38_high")
CC33.mid<-subset(CVIDagg6, subset= sub==2&cells=="CD38_mid")
CC33.low<-subset(CVIDagg6, subset= sub==2&cells=="CD38_low")

CC22.high<-subset(CVIDagg6, subset= sub==3&cells=="CD38_high")
CC22.mid<-subset(CVIDagg6, subset= sub==3&cells=="CD38_mid")
CC22.low<-subset(CVIDagg6, subset= sub==3&cells=="CD38_low")

ZA61.high<-subset(CVIDagg6, subset= sub==4&cells=="CD38_high")
ZA61.mid<-subset(CVIDagg6, subset= sub==4&cells=="CD38_mid")
ZA61.low<-subset(CVIDagg6, subset= sub==4&cells=="CD38_low")

CC27.high<-subset(CVIDagg6, subset= sub==5&cells=="CD38_high")
CC27.mid<-subset(CVIDagg6, subset= sub==5&cells=="CD38_mid")
CC27.low<-subset(CVIDagg6, subset= sub==5&cells=="CD38_low")

ZA63.high<-subset(CVIDagg6, subset= sub==6&cells=="CD38_high")
ZA63.mid<-subset(CVIDagg6, subset= sub==6&cells=="CD38_mid")
ZA63.low<-subset(CVIDagg6, subset= sub==6&cells=="CD38_low")

#now do integration
#now we have done the normalization and findVariableFeatures, do directly integration
ancs<-SelectIntegrationFeatures(
	list(ZA03.high, ZA03.mid, ZA03.low, CC33.high,CC33.mid,CC33.low,
		CC22.high,CC22.mid,CC22.low, ZA61.high, ZA61.mid, ZA61.low,
		CC27.high, CC27.mid, CC27.low, ZA63.high, ZA63.mid, ZA63.low	
			)
	, nfeatures=2500) #find the variable features for doing integration
ancs<-union(ancs, features.names) #add ones that we are interested.
CVIDagg6.anchors <- FindIntegrationAnchors(object.list = 
		list(ZA03.high, ZA03.mid, ZA03.low, CC33.high,CC33.mid,CC33.low,
		CC22.high,CC22.mid,CC22.low, ZA61.high, ZA61.mid, ZA61.low,
		CC27.high, CC27.mid, CC27.low, ZA63.high, ZA63.mid, ZA63.low	
			), 
		dims = 1:20, anchor.features=ancs)
CVIDagg6.int<- IntegrateData(anchorset = CVIDagg6.anchors, dims = 1:20)

#do scale and PCA
CVIDagg6.int<-ScaleData(CVIDagg6.int, verbose =T, vars.to.regress = c("nCount_RNA", "percent.mito"), features=rownames(CVIDagg6.int))

CVIDagg6.int<-RunPCA(CVIDagg6.int, npcs = 100, verbose = T)
print(CVIDagg6.int[["pca"]], dims = 1:10, nfeatures = 5)
#do clustering
ElbowPlot(object = CVIDagg6.int, ndims = 50)
CVIDagg6.int<- FindNeighbors(CVIDagg6.int, reduction = "pca", dims = 1:20)
CVIDagg6.int<- FindClusters(object = CVIDagg6.int, #reduction.type = "pca", 
                            dims= 1:20, resolution = .6#, #k.param = 30, print.output = TRUE, save.SNN = TRUE
)

#do UMAP
CVIDagg6.int <- RunUMAP(object = CVIDagg6.int, reduction = "pca", dims= 1:20, 
		min.dist = .75)

png(file="scRNACluster_IntergrateByCells.png",width=1000, height=1000)
DimPlot(CVIDagg6.int, reduction = "umap", label=T, group.by=c("cvid","cells","sub", "seurat_clusters"))
dev.off()
README="Batch correction/integration by cells with seurat"
README =paste0(README,"See code and details in CVIDagg6_SeuratPipeline_FF_v1.0.R")
save(CVIDagg6.int, README, file="scObjectSeurate.integrateCells.RData")
load(file="scObjectSeurate.integrateCells.RData")
###############################
###             MNN batch correction
#what to do next? need to do using scRNASeq and bachelor
#It can be used on seurat objects
#

CVIDagg6.int<-RunFastMNN(object.list=
	list(ZA03.high, ZA03.mid, ZA03.low, CC33.high,CC33.mid,CC33.low,
		CC22.high,CC22.mid,CC22.low, ZA61.high, ZA61.mid, ZA61.low,
		CC27.high, CC27.mid, CC27.low, ZA63.high, ZA63.mid, ZA63.low	
			),features = 5000  #originally this was by default 2000
	)
CVIDagg6.int <-RunUMAP(CVIDagg6.int, reduction = "mnn", 
	dims = 1:20, min.dist=0.75)
#saveRDS(file=here("Output","CVIDagg6_MNN_20K.rds"), CVIDagg6.int)
#show the diagnostics for the batch correction

CVIDagg6.int@tools$RunFastMNN$merge.info$lost.var %>% View("ff")
#http://bioconductor.org/books/3.15/OSCA.multisample/correction-diagnostics.html#mnn-specific-diagnostics
#This is returned via the lost.var field in the metadata of mnn.out, which contains a matrix of the variance lost in each batch (column) at each merge step (row).
#Large proportions of lost variance (>10%) suggest that correction is removing genuine biological heterogeneity. This would occur due to violations of the assumption of orthogonality between the batch effect and the biological subspace (Haghverdi et al. 2018). 

CVIDagg6.int <- FindNeighbors(CVIDagg6.int, reduction = "mnn", dims = 1:20)
CVIDagg6.int <- FindClusters(CVIDagg6.int,dims= 1:20, resolution = .6)
png(file="scRNACluster_IntergrateByCellsMNN.png",width=1000, height=1000)
DimPlot(CVIDagg6.int, group.by = c("cvid", "cells", "sub","ident") )
dev.off()
README="Batch correction/integration by cells with seurat using MNN bachelor"
README =paste0(README,"See code and details in CVIDagg6_SeuratPipeline_FF_v1.0.R")
save(CVIDagg6.int, README, 
	file=here("Output","scObjectSeurate.integrateCellsMNN.RData"))
load(file="scObjectSeurate.integrateCellsMNN.RData")
#now let's do cell cyle.
#cc.genes is in seurat library now.

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(CVIDagg6.int)<-"mnn.reconstructed"
CVIDagg6.int<-CellCycleScoring(CVIDagg6.int, 
	s.features = s.genes, 
	g2m.features = g2m.genes, 
	set.ident = TRUE)
head(CVIDagg6.int)
png(file="scRNACluster_MNN_cellcycle.png")
DimPlot(CVIDagg6.int,label=T, 
		group.by = c("cvid", "cells", "sub","Phase") )
dev.off()

png(file="scRNACluster_MNN_clusters.png")
DimPlot(CVIDagg6.int,label=T, 
		group.by = c("seurat_clusters"
		#,"RNA_snn_res.0.6"
		) )
		#RNA_snn_res.0.6 is the same as seurat_clusters
dev.off()
#note: it has been saved into above object. CVIDagg6.int MNN saved

CVIDagg6.int@meta.data$cells<-factor(CVIDagg6.int@meta.data$cells,
	levels=c("CD38_low","CD38_mid","CD38_high"))

png(file="scRNACluster_MNN_clustersByCells.png",
	width=1200, height=500)
DimPlot(CVIDagg6.int,label=T, 
		group.by = c("seurat_clusters"
		#,"RNA_snn_res.0.6"
		) ,split.by=c("cells"))
dev.off()


####doing alluvial gram

mdat<-CVIDagg6.int@meta.data
mdat.sub.cluster<-aggregate(mdat$sub, by=list(mdat$sub, mdat$cells, mdat$seurat_clusters), FUN=length, drop=F)
mdat.sub.cluster[is.na(mdat.sub.cluster$x),"x"]<-0
names(mdat.sub.cluster)<-c("sub", "cells","clusters","Freq")
mdat.sub.cluster$CVID<-"CVID"

mdat.sub.cluster[mdat.sub.cluster$sub==1,"CVID"]<-"HC"
mdat.sub.cluster[mdat.sub.cluster$sub==4,"CVID"]<-"HC"
mdat.sub.cluster[mdat.sub.cluster$sub==6,"CVID"]<-"HC"
library(ggalluvial)
###more example from https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
png(file="figure_alluvial_1.png", width=800, height=650)
ggplot(mdat.sub.cluster,
       aes(y = Freq, axis1 = sub, axis2 = clusters)) +  #need to change this line for different categories
  geom_alluvium(aes(fill = cells), width = 1/12) +  #change this line for different flows
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sub", "clusters"), expand = c(.05, .05)) +  #name the columns or categories
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Figure_Alluvial_1")
dev.off()
png(file="figure_alluvial_2.png", width=800, height=650)
ggplot(mdat.sub.cluster,
       aes(y = Freq, axis1 = sub, axis2=cells, axis3 = clusters)) +
  geom_alluvium(aes(fill = CVID), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sub", "cells","clusters"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Figure_Alluvial_2")
dev.off()


#now we need to draw a figure different showing stacked for each cluster
# show the composition of the cluster.

ggplot(mdat.sub.cluster, 
		aes(fill=CVID, y=Freq, x=clusters)) + 
    geom_bar(position="fill", stat="identity")

ggplot(mdat.sub.cluster, 
		aes(fill=cells, y=Freq, x=clusters)) + 
    geom_bar(position="fill", stat="identity")





#####################left over#####################
CVIDagg6.int.CD38low<-subset(CVIDagg6.int, subset= cells=="CD38_low")
#CVIDagg6.int.CD38low<-RunUMAP(object = CVIDagg6.int.CD38low, reduction.use = "pca", dims= 1:20, min.dist = .75)
low<-DimPlot(CVIDagg6.int.CD38low, reduction = "umap", label=T, 
	group.by=
		"seurat_clusters"
	)#+theme(title=element_blank())
	
CVIDagg6.int.CD38high<-subset(CVIDagg6.int, subset= cells=="CD38_high")
#CVIDagg6.int.CD38high<-RunUMAP(object = CVIDagg6.int.CD38high, reduction.use = "pca", dims= 1:20, min_dist = .75)
high<-DimPlot(CVIDagg6.int.CD38high, reduction = "umap", 
		label=T
		#, group.by=c("seurat_clusters")
		)

CVIDagg6.int.CD38mid<-subset(CVIDagg6.int, subset= cells=="CD38_mid")
#CVIDagg6.int.CD38mid<-RunUMAP(object = CVIDagg6.int.CD38mid, reduction.use = "pca", dims= 1:20, min_dist = .75)
mid<-DimPlot(CVIDagg6.int.CD38mid, reduction = "umap", 
		label=T#, group.by=c("seurat_clusters"))+
	) 

ggarrange(low, mid, high, nrow=2, ncol=2,
	labels=c("CD38_low","CD38_mid","CD38_high"),legend="right",
	common.legend=T)


#subset CVIDagg6.int by CD38 high/mid/low populations, then r

#Find Markers different between CVID and HC B cells within each CD38 high/mid/low population
CVIDagg6.CD38high.cvid.markers <- FindMarkers(CVIDagg6.int.CD38high, ident.1 = "CVID", group.by = "cvid", min.pct = 0.1)
WriteXLS(CVIDagg6.CD38high.cvid.markers, "Documents/PhD Thesis Work/Data/Single Cell/CVIDagg6.high.cvid.xlsx", row.names = TRUE, col.names = TRUE)

CVIDagg6.CD38mid.cvid.markers <- FindMarkers(CVIDagg6.int.CD38mid, ident.1 = "CVID", group.by = "cvid", min.pct = 0.1)
WriteXLS(CVIDagg6.CD38mid.cvid.markers, "Documents/PhD Thesis Work/Data/Single Cell/CVIDagg6.mid.cvid.xlsx", row.names = TRUE, col.names = TRUE)

CVIDagg6.CD38low.cvid.markers <- FindMarkers(CVIDagg6.int.CD38low, ident.1 = "CVID", group.by = "cvid", min.pct = 0.1)
WriteXLS(CVIDagg6.CD38low.cvid.markers, "Documents/PhD Thesis Work/Data/Single Cell/CVIDagg6.low.cvid.xlsx", row.names = TRUE, col.names = TRUE)


