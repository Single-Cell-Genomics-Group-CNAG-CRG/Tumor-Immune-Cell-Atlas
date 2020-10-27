# this script prepares de datasets (seurat objects) for integration
# by normalizing and finding the highly variable genes (HVG) of each dataset

# load libraries
library(tidyverse)
library(Seurat)

# load data (only immune cells) into a list of objects
obj_list <- character()
file_names <- c("list filenames of the seurat objects")
for (item in file_names) {
  temp <- readRDS(item)
  DefaultAssay(temp) <- "RNA"
  obj_list <- c(obj_list, temp)
}

# prepare integration by normalizing and finding the HVG of each dataset
obj_list <- sapply(obj_list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = 2000, verbose = FALSE, selection.method = "vst")
})

saveRDS(obj_list, "object_list.rds")
