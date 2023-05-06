# R code to do cell type annotation
#-----
# Start doing manual annotation
#	- plot expression of following markers
#		- CD19, CD20, CD21, CD24, CD27, CD38
#			CD79a, CD10, CD138, IgD, IgM
#	based on Sanz, et al,Front. Immunol., 18 October 2019
#		https://doi.org/10.3389/fimmu.2019.02458
#	updated 2/20/2023
#
#https://pubmed.ncbi.nlm.nih.gov/36380692/
#CD21 https://www.pnas.org/doi/full/10.1073/pnas.0901984106
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4955005/
#https://ashpublications.org/blood/article/115/3/519/27130/Differential-expression-of-CD21-identifies

#=========================
#	see automated cell annotation in 
#	annotation_v1.0.R
#	

library(scCATCH)
library(SCINA)
library(here)
library(dplyr)
library(Seurat)

library(SingleR)
library(celldex)
library(SingleCellExperiment)


#load data, assuming the 3 clusters are in the 
#Data folder.

data.dir<-"Data"

cvid.combined<-readRDS(file=here("Output","CVIDagg6_MNN_2K.rds"))
DefaultAssay(cvid.combined)<-"RNA"
VlnPlot(cvid.combined, c("CD27", "CD38", "CD24", "CR2", "CD79A","CD19"),
	 ncol = 3)

#########################
#	let's show markers
#########################
#options(repr.plot.height = 5, repr.plot.width = 18)
# CR2 = CD21
# CD79A = IgA
# SDC1 = CD138
# MS4A1 = CD20
# MME = CD10

png(file="general_B cellmark.png", width=1000, height=650)
VlnPlot(cvid.combined, c("CD19", "MS4A1", "CD27", 
		 "CD79A"), ncol = 4)#+xlab("General B cell marker")
dev.off()
png(file="Naive_TransitionMarker.png", width=1200, height=650)
VlnPlot(cvid.combined, c("CD38", "CR2", "CD24", 
		 "MME", "IGHD","IGHM","SDC1" 
		), ncol = 3)
dev.off()
png(file="IGs_marker.png", width=1200, height=650)
VlnPlot(cvid.combined, c("IGHG1","IGHG2","IGHG3","IGHG4",
		"IGHA1","IGHA2","IGHE", "IGHD", "IGHM" 
		), ncol = 3)
dev.off()
DefaultAssay(cvid.combined) <- "RNA"

cvid.combined@meta.data <- cvid.combined@meta.data %>%
	mutate(orig.clusters=seurat_clusters)

#now calling scTransform incorporated normalization, scale and variable genes
cvid.combined<-SCTransform(cvid.combined, vars.to.regress = c("sub", "nCount_RNA"))
#cvid.combined <- FindVariableFeatures(cvid.combined, nfeatures=5000)

cvid.combined<-RunPCA(cvid.combined, verbose = T)
ElbowPlot(cvid.combined, ndims = 50)

cvid.combined<-RunUMAP(cvid.combined, dims = 1:35, verbose = T)

cvid.combined <- FindNeighbors(cvid.combined, dims = 1:35, verbose = T)
cvid.combined <- FindClusters(cvid.combined, verbose = T, resolution = 0.8)


colours <- c(RColorBrewer::brewer.pal(n = 8, "Set2"), RColorBrewer::brewer.pal(n = 9, "Set1"))     

png(file="newClusters.png", width=650, height=650)
DimPlot(cvid.combined, label = F, group.by="orig.clusters")
dev.off()

png(file="features38_24_21_138.png", width=1000, height=1000)
FeaturePlot(cvid.combined,label=F, feature=c("CD38","CD24","CR2"
		, "MS4A1" 
		))
dev.off()

png(file="features_IgDMAE.png",width=1000, height=1000)
FeaturePlot(cvid.combined,label=F, feature=c("IGHD","IGHM","IGHA1","IGHE"), 
	pt.size=1.2, blend=F)
dev.off()

png(file="features_IgG1_4.png",width=1000,height=1000)
FeaturePlot(cvid.combined,label=F, feature=c("IGHG1","IGHG2","IGHG3","IGHG4"), 
	pt.size=1.2, blend=F)
dev.off()

png(file="featureCoExp_IgHMD.png", width=1000, height=600)
FeaturePlot(cvid.combined,label=F, feature=c("IGHM","IGHD"), 
	pt.size=3, blend=T)
dev.off()

png(file="featureCoExp_IgHMG3.png", width=1000, height=600)
FeaturePlot(cvid.combined,label=F, feature=c("IGHM","IGHG3"), 
	pt.size=3, blend=T)
dev.off()

png(file="featureCoExp_IgHDG3.png", width=1000, height=600)
FeaturePlot(cvid.combined,label=F, feature=c("IGHD","IGHG3"), 
	pt.size=3, blend=T)
dev.off()

###let's do do plots
#IGHM vs IGHD
dd<-GetAssayData(object = cvid.combined[["integrated"]], slot = "data")
dd<-as.matrix(dd)
dd<-t(dd)
dd<-as.data.frame(dd)
temp<-dd[,c("IGHM","IGHD","IGHG1","IGHG2",
		"IGHG3")]

#y<-dd["IGHD",]
x1<-ggplot(temp, aes(x=IGHM, y=IGHD))+
	geom_point()

x2<-ggplot(temp, aes(x=IGHM, y=IGHG3))+
	geom_point()

x3<-ggplot(temp, aes(x=IGHD, y=IGHG3))+
	geom_point()
png(file="2d_play_coexpression.png", width=900, height=900)
ggarrange(x1,x2,x3, ncol=2, nrow=2)
dev.off()