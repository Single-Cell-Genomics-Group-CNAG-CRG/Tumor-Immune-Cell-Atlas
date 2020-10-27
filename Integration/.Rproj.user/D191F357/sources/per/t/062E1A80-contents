# this script finds the anchors for the integration of all the datasets

# load libraries
library(tidyverse)
library(Seurat)

features <- readRDS("features.rds")
obj_list <- readRDS("processed_object_list.rds")

# use reciprocal principal components method
anchors <- FindIntegrationAnchors(obj_list,
  normalization.method = "LogNormalize",
  anchor.features = features,
  reduction = "rpca",
  verbose = T,
  dims = 1:20
)
rm(features)
saveRDS(anchors, file = "anchors.rds")
