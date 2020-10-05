# in this script we annotate the clusters in the atlas (integrated object)
# we curated a cluster names list with aid of experts in our team

# load libraries
library(tidyverse)
library(Seurat)

integrated <- readRDS("integrated_object_clustered.rds")

Idents(integrated) <- "seurat_clusters"

# we use the list of cluster names created with experts's help
cluster_names <- readRDS("names_list.rds")
names(cluster_names) <- levels(integrated)
integrated <- RenameIdents(integrated, cluster_names)

# save annotation in independent variable
integrated$cell_type <- Idents(integrated)

# we remove unknown/mixed and fibroblast clusters
integrated <- subset(integrated, 
                     cell_type %in% c("Unknown", "CAFs"), 
                     invert = TRUE)
saveRDS(integrated, file = "integrated_object_annotated.rds")
