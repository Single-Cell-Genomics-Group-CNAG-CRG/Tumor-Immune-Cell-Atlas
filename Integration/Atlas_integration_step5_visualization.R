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
library(matchSCore2)
library(cowplot)
library(magrittr)
source("myfunctions.R")

c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
         "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown") # distinct color palette

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

integrated <- readRDS("output/Atlas/integrated_object.rds")
integrated
integrated <- subset(integrated, state == "metastasis", invert = T)
integrated

counts <- GetAssayData(integrated, assay = "RNA")
grep(pattern = "^ERCC-", x = rownames(counts), value = T)
counts <- counts[-(which(rownames(counts) %in% grep(pattern = "^ERCC", x = rownames(counts), value = T))),]
integrated <- subset(integrated, features = rownames(counts))
integrated

integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, verbose = T)
integrated <- RunUMAP(integrated, dims = 1:30)
integrated$seurat_clusters <- NULL

# save integrated dataset (w/o cluster info)
saveRDS(integrated, file = paste(out_path, "Atlas_integrated_processed.rds", sep = ''))

q()

### VISUALIZATION --------
DimPlot(integrated, group.by = "source", pt.size =1, reduction = "tsne", cols = c25)
ggsave(paste(plot_path, date, "_integrated_tsne.png", sep = ''), width = 10)

DimPlot(integrated, group.by = "source", label = T, cols = as.vector(polychrome()))
ggsave(paste(plot_path, date, "_integrated_original_cell_types.png", sep = ''), width = 20)

DimPlot(integrated, pt.size = 1, cols = as.vector(polychrome()), split.by = "source", ncol = 3, label = T) +
  ggplot2::theme(legend.position = "bottom")
ggsave(paste(plot_path, date, "_integrated_split.png", sep = ''))