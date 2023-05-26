#R code to do pathway analysis using ReactomeGSA 
#
#--- started 5/23/2023
# ref: https://reactome.org/

#library(scde)
library(here)
library(dplyr)
library(org.Hs.eg.db)
library(Seurat)
library(ReactomeGSA)


#start loading data
# for scde/pagoda we need raw count (integer counts)

data.dir<-"Data"

cvid.combined<-readRDS(file=here("Output","CVIDagg6.int_2K_singleR.Rds"))
DefaultAssay(cvid.combined)<-"RNA"

#check to make sure the identification is the clusters
Idents(cvid.combined)


##
gsa_result <- analyse_sc_clusters(cvid.combined, verbose = TRUE)

saveRDS(gsa_result, here("Output","pathway_Reactome_gsa_result.Rds"))

gsa_result

gsa_pathways <- pathways(gsa_result)

#find the max diff for each pathway
# find the maximum differently expressed pathway
max_diff <- do.call(rbind, apply(gsa_pathways, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_diff$diff<-max_diff$max - max_diff$min
max_diff <- max_diff[order(max_diff$diff, decreasing = T), ]

head(max_diff)

#plotting for single pathway
plot_obj <- plot_gsva_pathway(gsa_result, "R-HSA-1169091")

#heatmap
selected_pathways <- c(
                       "R-HSA-1169091", "R-HSA-1168372", "R-HSA-983705", # BCR signalling
                       "R-HSA-9034013", "R-HSA-9024909", "R-HSA-9025046" # NTRK3 associated pathways 
                       )

# sort the cells aplhabetically in the pathway result
org_cols <- colnames(gsa_result@results[[1]]$pathways)
org_cols <- c(org_cols[1:2], sort(org_cols[3:length(org_cols)]))

gsa_result@results[[1]]$pathways <- gsa_result@results[[1]]$pathways[, org_cols]

# replace the "." from the cell names
colnames(gsa_result@results[[1]]$pathways) <- gsub("\\.", " ", colnames(gsa_result@results[[1]]$pathways))

# display the heatmap
#options(repr.plot.width = 15, repr.plot.height = 7)
plot_gsva_heatmap(gsa_result, truncate_names = F, pathway_ids = selected_pathways, 
                  dendrogram = "none", 
                  rowsep = 1:length(selected_pathways), colsep = 1:34, sepcolor = "black", sepwidth = c(0.01, 0.0001),
                  scale = "row",
                  margins = c(13, 25),
                  cexRow = 1,
                  xlab = "Cell clusters",
                  ylab = "Pathways",
                  Colv = F)
     
