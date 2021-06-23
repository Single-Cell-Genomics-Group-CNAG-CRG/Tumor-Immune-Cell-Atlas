# scale each dataset, compute its PCs and find the integration features

library(Seurat)
library(tidyverse)

source("code/functions/load_objects.R")

obj_list <- load_objects(directory="output/",
                         extension="_normalized.rds")

# select features that are repeatedly variable across datasets for integration 
# run PCA on each # dataset using these features
features <- SelectIntegrationFeatures(obj_list, 
                                      nfeatures = 5000)
print(head(features))

# but remove dataset specific
batch_markers <- readRDS("output/markers/batch_specific_markers.rds")
batch_markers <-batch_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 150, wt = avg_log2FC)

table(batch_markers$gene %in% features)

features <- setdiff(features, batch_markers$gene)
print(length(features))

# save selected features
saveRDS(features, "output/integration_features.rds")