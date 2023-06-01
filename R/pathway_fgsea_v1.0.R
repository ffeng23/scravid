#R code to do pathway analysis using differential 
#  gene expression and then fgsea 
#ref: https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html#gene_set_enrichment_(gsea)_analysis
#--- started 5/23/2023
# ref: https://reactome.org/

#library(scde)
library(here)
library(dplyr)
library(org.Hs.eg.db)
library(Seurat)
library(ReactomeGSA)
#install_github('immunogenomics/presto')
library(fgsea)
library(presto)
library(msigdbr)
library(ggplot2)
library(tidyverse)

#start loading data
# for scde/pagoda we need raw count (integer counts)

data.dir<-"Data"

cvid.combined<-readRDS(file=here("Output","CVIDagg6.int_2K_singleR.Rds"))
DefaultAssay(cvid.combined)<-"RNA"

#check to make sure the identification is the clusters
#Idents(cvid.combined)
 
#start calling wilcoxauc (presto), this is fast 
# and can do group vs non-group comparison 
# we first call on cvid vs healthy

de.cvid <- wilcoxauc(cvid.combined, group_by='cvid')
head(de.cvid)

dplyr::count(de.cvid, group)
#

#human gene, with category C7 is immunological genes
# C7 5219 gene sets in total
m_df<- msigdbr(species = "Homo sapiens", category = "C7")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

head(fgsea_sets)

fgsea_sets$GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_UP

#get a rank so to do gsea
# select only the feature and auc columns for fgsea, which statistics to use is an open question
cluster0.genes<- de.cvid %>%
  dplyr::filter(group == "CVID") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(cluster0.genes)

head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks,
		eps=0, scoreType="pos"
	)#, nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()
write_tsv(x=as.data.frame(fgseaRes),
	file=here("Output","fGSEA_results_cvid_hc.tsv"))

#plot
# only plot the top 20 pathways
ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% 
head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plotEnrichment(fgsea_sets[["GSE10325_CD4_TCELL_VS_MYELOID_UP"]],
               ranks) + labs(title="GSE10325 CD4 TCELL VS MYELOID UP")

### doing cell population differences
de.cells <- wilcoxauc(cvid.combined, group_by='cells')
head(de.cells)

dplyr::count(de.cells, group)
#
# select only the feature and auc columns for fgsea, which statistics to use is an open question
cluster0.genes<- de.cells %>%
  dplyr::filter(group == "CD38_high") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(cluster0.genes)

head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks,
		eps=0, scoreType="pos"
	)#, nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()
write_tsv(x=as.data.frame(fgseaRes),
	file=here("Output","fGSEA_results_cd38_hilowmid.tsv"))

#plot
# only plot the top 20 pathways
ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% 
head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plotEnrichment(fgsea_sets[["GSE22886_IGG_IGA_MEMORY_BCELL_VS_BM_PLASMA_CELL_UP"]],
               ranks) + labs(title="GSE22886_IGG_IGA_MEMORY_BCELL_VS_BM_PLASMA_CELL_UP")


########now we need to do clusters.
### doing cell population differences
de.clusters <- wilcoxauc(cvid.combined, group_by='seurat_clusters')
head(de.clusters)

dplyr::count(de.clusters, group)
#
# select only the feature and auc columns for fgsea, which statistics to use is an open question
cluster0.genes<- de.clusters %>%
  dplyr::filter(group == "2") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


ranks<- deframe(cluster0.genes)

head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks,
		eps=0, scoreType="pos"
	)#, nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()

write_tsv(x=as.data.frame(fgseaRes),
	file=here("Output","fGSEA_results_clusters_0To34.tsv"))

#plot
# only plot the top 20 pathways
ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% 
head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plotEnrichment(fgsea_sets[["GSE22886_IGG_IGA_MEMORY_BCELL_VS_BM_PLASMA_CELL_UP"]],
               ranks) + labs(title="GSE22886_IGG_IGA_MEMORY_BCELL_VS_BM_PLASMA_CELL_UP")
