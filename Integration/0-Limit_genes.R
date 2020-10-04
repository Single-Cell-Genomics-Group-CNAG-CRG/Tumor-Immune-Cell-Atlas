# this script is used to compute the dataset-specific markers
# we will substract these markers from the anchoring set
# this will make the integration less noisy and more robust
# the input this script takes is a list of the seurat objects that will be integrated

# load libraries
library(tidyverse)
library(Seurat)

# the list of seurat objects
obj_list <- readRDS("path to list")

obj <- merge(obj_list[[1]], obj_list[2:length(obj_list)], merge.data = F)
rm(obj_list)

# source variable contains the dataset of origin of each cell
Idents(obj) <- "source"

markers <- FindAllMarkers(obj, verbose = T, only.pos = T, max.cells.per.ident = 5000)
saveRDS(markers, file = "selected_markers.rds")
