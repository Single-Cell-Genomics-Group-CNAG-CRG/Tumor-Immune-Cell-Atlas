# find the anchors to perform the integration

library(Seurat)
library(tidyverse)

source("code/functions/load_objects.R")

obj_list <- load_objects(
  directory = "output/",
  extension = "_scaled.rds"
)

# load selected integration features

features <- readRDS("output/integration_features.rds")

anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                  normalization.method = "LogNormalize",
                                  anchor.features = features, 
                                  reduction = "rpca")

# save anchors
saveRDS(anchors, "output/anchors.rds")
