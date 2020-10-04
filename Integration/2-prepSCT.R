# in this script we further process the datasets prior to integration
# by calculating the integration features
# then scaling and calculationg the principal components of each dataset


# load libraries
library(tidyverse)
library(Seurat)

obj_list <- readRDS("object_list.rds")
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000, verbose = TRUE)

# substract source-specific features
selected_genes <- readRDS("selected_markers.rds")
features <- setdiff(features, selected_genes)
saveRDS(features, file = "features.rds")

obj_list <- lapply(X = obj_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs = 30)
})

saveRDS(obj_list, file = "processed_object_list.rds")
