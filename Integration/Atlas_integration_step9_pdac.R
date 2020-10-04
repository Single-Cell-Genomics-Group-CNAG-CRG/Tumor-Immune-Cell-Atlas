# set working directory
rm(list=ls())
setwd("/home/devel/pnieto/TIL_Atlas")

# load libraries
library(tidyverse)
library(Seurat)
#library(future)
#library(sctransform)
#library(ggsci)
#library(pals)
#library(RColorBrewer)
#source("myfunctions.R")
`%notin%` <- Negate(`%in%`)
#values = as.vector(polychrome())

# Enable parallelization
#options(future.globals.maxSize = +Inf)
#plan("multiprocess", workers = 24)

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

atlas <- readRDS(file = paste(out_path, "Atlas_patient_filtered_k6.rds", sep = ''))
atlas$general_cell_type <- as.character(atlas$general_cell_type)
atlas$specific_cell_type <- as.character(atlas$specific_cell_type)

atlas$general_cell_type[atlas$general_cell_type %in% c("CD4", "CD8")] <- "T cell"
atlas$general_cell_type[atlas$general_cell_type %in% c("mDC", "cDC", "pDC")] <- "DC"

# go back to RNA assay and normalize!
DefaultAssay(atlas) <- "RNA"
# Normalize RNA data for visualization purposes
#atlas <- NormalizeData(atlas, verbose = T, assay = "RNA")

atlas$general_cell_type <- as.factor(atlas$general_cell_type)
Idents(atlas) <- "general_cell_type"
table(atlas$general_cell_type)
atlas_general_markers <- FindAllMarkers(atlas, only.pos = T, assay = "RNA", test.use = "MAST")
print(table(atlas_general_markers$cluster))
atlas_general_markers %<>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
saveRDS(atlas_general_markers, paste0(out_path, "top200_general_markers_logFC.rds"))

Idents(atlas) <- "specific_cell_type"
table(atlas$specific_cell_type)
atlas$general_cell_type <- as.factor(atlas$specific_cell_type)
atlas_specific_markers <- FindAllMarkers(atlas, only.pos = T, assay = "RNA", test.use = "MAST")
print(table(atlas_specific_markers$cluster))
atlas_specific_markers %<>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
saveRDS(atlas_specific_markers, paste0(out_path, "top200_specific_markers_logFC.rds"))

atlas$general_cell_type <- as.character(atlas$general_cell_type)
atlas <- subset(atlas, general_cell_type == "T cell")
table(atlas$specific_cell_type)
Idents(atlas) <- "specific_cell_type"
atlas$general_cell_type <- as.factor(atlas$specific_cell_type)
atlas_tcell_markers <- FindAllMarkers(atlas, only.pos = T, assay = "RNA", test.use = "MAST")
print(table(atlas_tcell_markers$cluster))
atlas_tcell_markers %<>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
saveRDS(atlas_tcell_markers, paste0(out_path, "top200_tcell_markers_logFC.rds"))

q()

#saveRDS(integrated, file = paste(out_path, "Atlas_patient_filtered_k6.rds", sep = ''))

#unique(integrated$state)
#integrated
#Idents(integrated) <- "specific_cell_type"
#integrated <- subset(integrated, downsample = 500)
#write.table(integrated@active.ident, file=paste0(out_path, 'Cells_label.tsv'), quote=FALSE, sep='\t', col.names = TRUE)
#write.table(integrated@assays[["RNA"]]@counts, file=paste0(out_path, 'Gene_Count_per_Cell.tsv'), quote=FALSE, sep='\t', col.names = TRUE)

#integrated$cluster_kmeans_k4[integrated$cluster_kmeans_k4 %notin% c(1,2,3,4)] <- "NA"
#integrated <- subset(integrated, cluster_kmeans_k4 == "NA", invert = T)
#table(integrated$cluster_kmeans_k4)
#saveRDS(integrated, file = paste(out_path, "Atlas_complete_annotation_k4_clusters.rds", sep = ''))

#print("general cell type")
#integrated$general_cell_type <- as.factor(integrated$general_cell_type)
#Idents(integrated) <- "general_cell_type"
#print(unique(Idents(integrated)))

#markers <- FindAllMarkers(integrated, verbose = T, only.pos = T, logfc.threshold = 0, min.pct = 0)
#markers <- pull(top_n(group_by(markers, cluster), n = 150, wt = avg_logFC), gene)
#markers <- matchSCore2::cut_markers(levels(markers$cluster), markers, ntop=150) 
#saveRDS(markers, file = paste(out_path, "Atlas_general_markers_full_list.rds", sep = ''))
#print("done with general markers")

###

print("specific_cell_type")
integrated$specific_cell_type <- as.factor(integrated$specific_cell_type)
Idents(integrated) <- "specific_cell_type"
print(unique(Idents(integrated)))

markers <- FindAllMarkers(integrated, verbose = T, only.pos = T, logfc.threshold = 0, min.pct = 0)
#markers <- pull(top_n(group_by(markers, cluster), n = 150, wt = avg_logFC), gene)
#markers <- matchSCore2::cut_markers(levels(markers$cluster), markers, ntop=150) 
saveRDS(markers, file = paste(out_path, "Atlas_specific_markers_full_list.rds", sep = ''))
print("done with specific markers")

###

integrated <- subset(integrated, cluster_kmeans_k4 %in% c(1,2,3,4))

print("kmeans clusters")
integrated$cluster_kmeans_k4 <- as.factor(integrated$cluster_kmeans_k4)
Idents(integrated) <- "cluster_kmeans_k4"
print(unique(Idents(integrated)))

markers <- FindAllMarkers(integrated, verbose = T, only.pos = T, logfc.threshold = 0, min.pct = 0)
saveRDS(markers, file = paste(out_path, "Atlas_kmeans_clusters_markers_full_list.rds", sep = ''))
print("done with cluster markers")

q()

pancreas <- subset(integrated, source == "pancreas")
#pancreas@meta.data %<>% select("nCount_RNA", "nFeature_RNA", "percent.mt", "cell_type")
saveRDS(pancreas, file = paste(out_path, "pancreas_new_celltypes.rds", sep = '')) # for Marc
pancreas <- subset(pancreas, downsample = 100)
#table(pancreas$new_assigned_cell_type)
saveRDS(pancreas, file = paste(out_path, "pancreas_new_celltypes_100.rds", sep = ''))  # for Marc
unique(Idents(pancreas))

