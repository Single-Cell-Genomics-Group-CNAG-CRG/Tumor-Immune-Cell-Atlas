rm(list=ls())
setwd("/home/devel/pnieto/TIL_Atlas")

# load libraries
library(tidyverse)
library(Seurat)
library(sctransform)
library(ggsci)
library(pals)
library(RColorBrewer)
library(matchSCore2)
library(magrittr)
source("myfunctions.R")
library(future)

# Enable parallelization
options(future.globals.maxSize = +Inf)
plan("multiprocess", workers = 24)

c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
         "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown") # distinct color palette

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

#integrated <- readRDS("output/Atlas/Atlas_integrated_processed.rds")

# CLUSTERING 
#integrated <- FindNeighbors(integrated, dims = 1:30)
#integrated <- FindClusters(integrated, verbose = T)

#saveRDS(integrated, file = paste(out_path, "Atlas_integrated_clustered.rds"))
integrated <- readRDS(file = paste(out_path, "Atlas_integrated_clustered.rds"))

# go back to RNA assay and normalize!
DefaultAssay(integrated) <- "RNA"
# Normalize RNA data for visualization purposes
integrated <- NormalizeData(integrated, verbose = T, assay = "RNA")

integrated$seurat_clusters <- as.factor(integrated$seurat_clusters)
Idents(integrated) <- integrated$seurat_clusters

#saveRDS(integrated, file = paste(out_path, "Atlas_integrated_clustered.rds"))

options(future.globals.maxSize = +Inf)
plan("multiprocess", workers = 24)

markers <- FindAllMarkers(integrated, verbose = T, only.pos = T, test.use = "MAST", assay = "RNA")
cluster_markers <- cut_markers(levels(markers$cluster), markers, ntop=25)
print("MAST")
print(cluster_markers)
markers <- filter(markers, p_val_adj < 0.05)
saveRDS(markers, file = paste(out_path, "Atlas_markers_MAST.rds"))

markers <- FindAllMarkers(integrated, verbose = T, only.pos = T, assay = "RNA")
cluster_markers <- cut_markers(levels(markers$cluster), markers, ntop=25)
print("wilcox")
print(cluster_markers)
markers <- filter(markers, p_val_adj < 0.05)
saveRDS(markers, file = paste(out_path, "Atlas_markers_wilcox.rds"))

#saveRDS(integrated, file = paste(out_path, "Atlas_integrated_clustered.rds"))