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
library(readr)
library(PerformanceAnalytics)


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

stat.cell<-CVIDagg6@meta.data %>% 
	group_by(sub,cvid) %>%
	summarize(cell_counts=n())
stat.seq<-CVIDagg6@meta.data %>% 
	group_by(sub,cvid) %>%
	summarize(total_counts=sum(nFeature_RNA), 
		unique_counts=sum(nCount_RNA))

stat<-stat.cell %>% left_join(stat.seq,by=c("sub","cvid"))
chart.Correlation(stat[,c("cell_counts", 
	"total_counts","unique_counts")])
#get number of cells and sequences for each subject
muts<-read_csv("/home/feng/Windows/windowsD/feng/LAB/MSI/maglione/SComatic_output/counts_mutations.csv")
stat <-stat %>% mutate(sub=as.integer(sub))
muts.stat<-stat %>% left_join(muts,by=c("sub"))

muts.stat <- muts.stat %>% 
	mutate(TACI_muts_normed= TACI_muts/cell_counts*1000,
		total_C17_normed=Total_muts_passed_C17/cell_counts*1000,
		total_C22_normed=Total_muts_passed_C22/cell_counts*1000
		)

ggplot(data=muts.stat, aes(y=TACI_muts_normed,x=sub,colour=as.factor(cvid)))+
	geom_boxplot() + xlab("Subjects")+ylab("Frequency/1K cells")+
	theme(legend.position="none")

ggplot(data=muts.stat, aes(y=total_C17_normed,x=sub,colour=as.factor(cvid)))+
	geom_boxplot()+ xlab("Subjects")+ylab("Frequency/1K cells")+
	theme(legend.position="none")

ggplot(data=muts.stat, aes(y=total_C22_normed,x=sub,colour=as.factor(cvid)))+
	geom_boxplot()+ xlab("Subjects")+ylab("Frequency/1K cells")+
	theme(legend.position="none")


