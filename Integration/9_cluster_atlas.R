# cluster atlas with optimal resolution obtained from the random forest oob

library(Seurat)
library(tidyverse)
# library(future)
# plan("multiprocess")
# options(future.globals.maxSize = Inf)

atlas <- readRDS("output/integrated_processed.rds")

atlas

atlas <- FindClusters(atlas, 
                      resolution = 1.2)

saveRDS(atlas, "output/integrated_clustered.rds")

p1 <- DimPlot(atlas, 
              reduction = "umap", 
              group.by = "source", 
              label = TRUE, 
              repel = TRUE)

p2 <- DimPlot(atlas, 
              reduction = "umap", 
              label = TRUE, 
              repel = TRUE)

plots <- p1 + p2

ggsave(plot = plots, 
       filename = "output/cluster_plots.png",
       dpi = 300,
       width = 11)

p3 <- DimPlot(atlas, 
              reduction = "umap", 
              split.by = "source",
              ncol = 3)

ggsave(plot = p3, 
       filename = "output/cluster_plots_split.png",
       dpi = 300,
       height = 12,
       width = 12)

table(atlas$seurat_clusters, atlas$source)