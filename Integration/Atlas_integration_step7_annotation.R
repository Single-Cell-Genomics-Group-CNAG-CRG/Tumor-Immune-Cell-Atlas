### in this script we will integrate the small (subsampeled) immune datasets
### more datasets will be added

# set working directory
rm(list=ls())
setwd("/home/devel/pnieto/TIL_Atlas")

# load libraries
library(tidyverse)
library(Seurat)
library(sctransform)
library(ggsci)
library(pals)
library(RColorBrewer)
source("myfunctions.R")

c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
         "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown") # distinct color palette

values = as.vector(polychrome())

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

integrated <- readRDS("output/Atlas/Atlas_integrated_clustered.rds")
print(integrated)

FeaturePlot(integrated, features = c("IL7R", "CD4", "CD8A", "CD8B", "CD68", "CD79A", "TAGLN", "CD1C", "CLEC4C"), label = T, combine=TRUE)
ggsave(paste(plot_path, date, "_lineage_markers.png", sep = ''), width = 15)

DotPlot(integrated, features = c("IL7R", "CD4", "CD8A", "CD8B", "CD68", "CD79A", "TAGLN", "CD1C", "CLEC4C")) + RotatedAxis()
ggsave(paste(plot_path, date, "_lineage_markers_dot.png", sep = ''))

###

DimPlot(integrated, cols = values, group.by = "seurat_clusters") + theme_void()
ggsave(paste(plot_path, date, "_Atlas_clusters.png", sep = ''))

DimPlot(integrated, cols = values, group.by = "seurat_clusters", split.by = "source", ncol = 3) + theme_void()
ggsave(paste(plot_path, date, "_Atlas_clusters_split_by_source.png", sep = ''))

DimPlot(integrated, cols = values, group.by = "source") + theme_void()
ggsave(paste(plot_path, date, "_Atlas_source.png", sep = ''))

###

Idents(integrated) <- "seurat_clusters"
new_names <- c("CD8", "CD4", "CD4", "NK", "CD8", "B cell", "CD4", "Myeloid", "CD4", "Myeloid", "Myeloid", #10
		"CD4", "CD4", "Myeloid", "B cell", "CD8", "Myeloid", "DC", "B cell", "CD8", "NK", #20
		"CD4", "Myeloid", "CD4", "NK", "Myeloid", "DC", "Myeloid", "Myeloid", "B cell", "CD4", #30
		"B cell", "Other", "B cell")
names(new_names) <- levels(integrated)
integrated <- RenameIdents(integrated, new_names)

###

markers <- FindAllMarkers(integrated, test.use = "MAST", verbose = T, min.pct = 0, min.diff.pct = 0)
markers <- matchSCore2::cut_markers(levels(markers$cluster), markers, ntop=100) 
saveRDS(markers, file = paste(out_path, date, "_Atlas_markers.rds", sep = ''))
q()

LabelClusters(DimPlot(integrated, pt.size = 1, cols = c25), "ident", fontface = "bold") + 
  guides(color=guide_legend(ncol=1)) + 
  theme_void()
ggsave(paste(plot_path, date, "named_clusters.png", sep = ''), width = 12, height = 10)

p1 <- DimPlot(integrated, pt.size = 1, cols = values, label = T) + theme_void()
p2 <- DimPlot(integrated, pt.size = 1, cols = values, group.by = "seurat_clusters", label = T) + theme_void()
CombinePlots(plots = c(p1, p2), ncol = 2)
ggsave(paste(plot_path, date, "_Atlas_combined_plots.png", sep = ''))

p1 <- DimPlot(integrated, pt.size = 1, cols = values, label = T, spli.by = "source", ncol = 3) + theme_void()
p2 <- DimPlot(integrated, pt.size = 1, cols = values, group.by = "seurat_clusters", split.by = "source", label = T, ncol = 3) + theme_void()
CombinePlots(plots = c(p1, p2), ncol = 2)
ggsave(paste(plot_path, date, "_Atlas_combined_plots_split.png", sep = ''))


# SAVE DATASETS --------------
#integrated <- readRDS(file = paste(out_path, date, "_Atlas_integrated_full.rds", sep = ''))
integrated$new_cell_type <- Idents(integrated)
integrated <- subset(integrated, new_cell_type %in% c("CD4", "CD8", "NK", "B cell", "Myeloid", "DC"))
saveRDS(integrated, file = paste(out_path, date, "_Atlas_integrated_full.rds", sep = ''))

CD4 = subset(integrated, idents = c("CD4"))
saveRDS(CD4, file = paste(out_path, date, "_Atlas_CD4.rds", sep = ''))
rm(CD4)

CD8 = subset(integrated, idents = c("CD8", "NK"))
saveRDS(CD8, file = paste(out_path, date, "_Atlas_CD8.rds", sep = ''))
rm(CD8)

B_cells = subset(integrated, idents = c("B cell"))
saveRDS(B_cells, file = paste(out_path, date, "_Atlas_Bcells.rds", sep = ''))
rm(B_cells)

macro = subset(integrated, idents = c("Myeloid", "DC"))
saveRDS(macro, file = paste(out_path, date, "_Atlas_Myeloid.rds", sep = ''))
rm(macro)

rest = subset(integrated, idents = c("Other"))
saveRDS(rest, file = paste(out_path, date, "_Atlas_Rest.rds", sep = ''))
rm(rest)
