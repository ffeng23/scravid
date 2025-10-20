#R code to do analysis after pagoda pathway 
# first round: 
# to show the pathways differentially up- or down-regulated
#

#here we started from pathway pca analysis.
# we do by creating group and test 
# for difference by T or wilcoxn.

#library



#load data

library(scde)
library(here)
library(dplyr)
library(org.Hs.eg.db)
library(Seurat)
library(SeuratObject)
library(presto)
library(GOfuncR)

#start loading data
# for scde/pagoda we need raw count (integer counts)

output.dir<-"Output"

cat("Read data..........\n")

pwpca<- readRDS(file=here(output.dir, "Output_allCells_1st_pwpca.Rds"))

#also need to load in the seurat object of RNASeq 
# we need the meta.data for clustering
# 
#data.dir<-"Data"

cvid.combined<-readRDS(file=here(output.dir,"CVIDagg6.int_2K_singleR.Rds"))
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

#get cluster info
mdat<-cvid.combined@meta.data
varinfo<-readRDS(file=here("Output","varinfo_allCells_1st.Rds"))

#
#x<-pwpca[[1]]$xp$scores

y<-lapply(pwpca, FUN=function(x){x$xp$scores})
#y<-varinfo$mat#lapply(varinfo, FUN=function(x){x$mat})

x<-data.frame(y)
x<-t(x)
remove(y)
remove(cvid.combined)
gc()

y<-wilcoxauc(x ,
	mdat[colnames(x),"seurat_clusters"])
y$feature<-sub(x=y$feature, pattern=".",":",fixed=T)

y_c2<- y %>% dplyr::filter(group==2)

y_c2$GOname<-get_names(y_c2$feature)$go_name

#do heatmap
library(pheatmap)
y_c2.sig<- y_c2 %>% dplyr::filter(padj<0.0001)%>% 
dplyr::filter(logFC>3 | logFC< -3)

#get the data matrix
x.dat<-x
rownames(x.dat)<-sub(x=rownames(x.dat),pattern=".", ":",fixed=T)

x.dat<-x.dat[y_c2.sig$feature,]

mdat$c2<-0
mdat[mdat$seurat_clusters==2,"c2"]<-1
mdat.2<-mdat[colnames(x),]
mdat.2<-rbind(mdat.2[mdat.2$c2==0,][1:2000,], mdat.2[mdat.2$c2==1,])

x.dat<-as.data.frame(x.dat) %>% dplyr::select(rownames(mdat.2))

pheatmap(x.dat,cluster_rows=T, cluster_cols=F,show_colnames=F, 
	scale="row" )
#<-


#doing pca 
# translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids 
# convert GO lists from ids to gene names
gos.interest <- unique(ls(org.Hs.egGO2ALLEGS))#[1:100],"GO:0022008","GO:0048699", "GO:0000280", "GO:0007067")) 
go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment

tam <- pagoda.top.aspects(pwpca,  n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
# determine overall cell clustering
#hc <- pagoda.cluster.cells(tam, varinfo)
#saveRDS(hc, file=here("Output","Output_hc_allcells_1st.Rds"))
hc<-readRDS(file=here("Output","Output_hc_allcells_1st.Rds"))

pagoda.show.pathways(c("GO:0030658","GO:0030666"), varinfo, go.env, cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)

get_names("GO:0030658")
get_names("GO:0030666")