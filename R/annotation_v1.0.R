# R code to do cell type annotation

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
load(here(data.dir,"CVIDagg6_Low_Cluster2.rdata"))
		#loaded the CVIDagg6.Low.Cluster2

load(here(data.dir,"CVIDagg6_Low_Cluster6.rdata"))
		#loaded the CVIDagg6.Low.Cluster6

load(here(data.dir,"CVIDagg6_Low_Cluster7.rdata"))
		#loaded the CVIDagg6.Low.Cluster2
cvid.combined<-merge(CVIDagg6.Low.Cluster2, 
		y=c(CVIDagg6.Low.Cluster6,CVIDagg6.Low.Cluster7),
		project="CVID.merge")

#########################
#	do marker based annotation
#########################
#scCATCH first
#now we need to make a maker list (cellmatch)

obj <- createscCATCH(data = cvid.combined[['RNA']]@data, 
		cluster = as.character(Idents(cvid.combined)))
obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch, tissue = "Peripheral blood")

obj <- findcelltype(object = obj)
obj@celltype %>% View("m")

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
exp<-cvid.combined[['RNA']]@data
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

#turn into a sce to do singleR
cvid6.c2<-as.SingleCellExperiment(CVIDagg6.Low.Cluster2) #cvid.combined)

#run assignment
pred.cvid6.c2 <- SingleR(test = cvid6.c2, ref = hpca.se, assay.type.test=1,
    labels = hpca.se$label.main)

pred.cvid6.c2.fine <- SingleR(test = cvid6.c2, ref = hpca.se, assay.type.test=1,
    labels = hpca.se$label.fine)

table(pred.cvid6.c2.fine$labels) %>% View("m3")
#now add back the labels to the original seurat object

CVIDagg6.Low.Cluster2@meta.data$label.main<-pred.cvid6.c2[rownames(CVIDagg6.Low.Cluster2@meta.data),"labels"]
CVIDagg6.Low.Cluster2@meta.data$label.fine<-pred.cvid6.c2.fine[rownames(CVIDagg6.Low.Cluster2@meta.data),"labels"]

#do this to other 2 clusters please!!!
