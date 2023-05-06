# R code to do cell type annotation
#  -- NOTE: see the 
#"file:///home/feng/Feng/hg/scravid/Output/annotation_v1.1_allNew.Rmd"
# used to annotate the newly clustered data
# (4/25/2023) for all clusters
#  The data were generated in 
#../../scravid2/CVIDagg6_SeuratPipeline_FF_v1.0.R
#
#-------------------------------
#   --- minor change
#=========================
# using different methods
#  - singleR
#  - scina  , marker based
#  - scCATCH , cell match
#  - scDeepSort, deep learning and python.
#
#		- started 2/10/2023  by Feng
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

#cluster 2, scRNASeq data.
#load(here(data.dir,"scObject.intergrateMNN.RData"))
		#loaded 
		# the data were generated in ../../scravid2/CVIDagg6_SeuratPipeline_FF_v1.0.R
#now we load a newly saved data which has more features
CVIDagg6.int<-readRDS(file=here("Output","CVIDagg6_MNN_2K.rds"))

#########################
#	do marker based annotation
#########################
#scCATCH first
#now we need to make a maker list (cellmatch)

obj <- createscCATCH(data = CVIDagg6.int[["RNA"]]@data, 
		cluster = as.character(CVIDagg6.int$seurat_clusters))
obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch, tissue = "Peripheral blood")

obj <- findcelltype(object = obj)
obj@celltype %>% View("m")
saveRDS(file=here("Output","scCATCH_annotation.Rds"),obj)
#scina
#need to make a signature (markers)
signature.db<- cellmatch %>% select(-pmid) %>%
	filter(species=="Human", tissue=="Peripheral blood",
		cancer=="Normal",condition=="Normal cell"#,
		#resource=="Single-cell sequencing"
		) %>% #filter(celltype=="B Cell")
		group_by(celltype, subtype1, subtype2, subtype3) %>%
		reframe(gene) 
signature.db <- signature.db %>% 
		mutate(celltype.new=paste(celltype,subtype1, subtype2, subtype3)) %>%
		mutate(celltype.new=gsub(x=celltype.new,
	pattern=" NA",replacement="", fixed=T))		

#make celltypes by combining the 
 celltypes<-unique(signature.db$celltype.new)

signatures.list<-list()
for(i in 1:length(celltypes))
{
	#get the list of all markers
	celltype.name<-celltypes[i]
	
	temp<-
		signature.db %>% filter(celltype.new==celltype.name)

	signatures.list[[ celltype.name ]]<-temp$gene

}

signatures.list<-lapply(signatures.list, FUN=unique)
#s<-list(signatures.list[[1]], signatures.list[[7]], signatures.list[[8]])

s<-signatures.list[unlist(
		lapply(signatures.list, function(x) {length(x)>4}) )  ]
lapply(s, length)
#s<-s[1:10]
#get exp
#exp<-cvid.combined[['RNA']]@data
exp<-CVIDagg6.int[['RNA']]@data
exp<-as.matrix(exp)
results = SCINA(exp, s, max_iter = 100, convergence_n = 10, 
    convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
    rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')

View(results$cell_labels)
#show summary of the labels. still not correct
table(results$cell_labels) %>% View("m3")
#View(results$probabilities)

# get rid of the signature genes not in the exp expression array
s.rm<-lapply(s, FUN=function(x){
		x<-x[is.element(x,rownames(exp))]
	})

##BE CAREFUL!! this will take very long.
plotheat.SCINA(exp, results, s.rm)

#doing singleR

#in this celldex package there are a bunch of reference data set.
#           see http://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html
#           for more information. Please at least try two of the below to see 
#           which one gives the most reasonable assignment.
hpca.se <- HumanPrimaryCellAtlasData()
imm.ref<-celldex::MonacoImmuneData()
DefaultAssay(CVIDagg6.int)<-"RNA"
#turn into a sce to do singleR
cvid6.c2<-as.SingleCellExperiment(CVIDagg6.int) #cvid.combined)

#run assignment
pred.cvid6.c2 <- SingleR(test = cvid6.c2, ref = hpca.se, assay.type.test=1,
    labels = hpca.se$label.main)

pred.cvid6.c2.fine <- SingleR(test = cvid6.c2, ref = hpca.se, assay.type.test=1,
    labels = hpca.se$label.fine)
#plot to see the quality of prediction
#ref:http://bioconductor.org/books/3.15/OSCA.basic/cell-type-annotation.html
plotScoreHeatmap(pred.cvid6.c2)
plotScoreHeatmap(pred.cvid6.c2.fine)

table(pred.cvid6.c2.fine$labels) %>% View("m3")
table(pred.cvid6.c2$labels) %>% View("m3")
#now add back the labels to the original seurat object

CVIDagg6.int@meta.data$label.main<-pred.cvid6.c2[rownames(CVIDagg6.int@meta.data),"labels"]
CVIDagg6.int@meta.data$label.fine<-pred.cvid6.c2.fine[rownames(CVIDagg6.int@meta.data),"labels"]

library(pheatmap)
#ref http://bioconductor.org/books/3.15/OSCA.basic/cell-type-annotation.html
tab <- table(Assigned=pred.cvid6.c2$labels, Cluster=CVIDagg6.int@meta.data$seurat_clusters))
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))


#do this to other 2 clusters please!!!
##############################################
#
#   NOTE: check the .Rmd file for the functions/code 
#file:///home/feng/Feng/hg/scravid/Output/annotation_v1.1_allNew.Rmd
# it combines all methods, singleR and scDeepSort
################################################

