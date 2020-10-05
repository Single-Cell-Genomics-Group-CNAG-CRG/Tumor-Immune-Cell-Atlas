# this scripts calculates the clusters and the cluster markers of the atlas
# we use MAST and Wilcoxon Rank Sum Test and then combine the results

# load libraries
library(tidyverse)
library(Seurat)
library(magrittr)

integrated <- readRDS("integrated_object_processed.rds")

# find neighboors and cluster
integrated %<>% FindNeighbors(dims = 1:30) %>%
  FindClusters(resolulution = 1.2)

# go back to RNA assay and normalize (this is the recommended pipeline
# to find cluster markers after integration, according to Seurat developers)
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated,
  verbose = TRUE,
  assay = "RNA"
)

integrated$seurat_clusters <- as.factor(integrated$seurat_clusters)
Idents(integrated) <- integrated$seurat_clusters

saveRDS(integrated, file = "integrated_object_clustered.rds")

# for more robustness we are going to use two statistical
# tests to find the cluster markers
markers_mast <- FindAllMarkers(integrated,
  only.pos = TRUE,
  test.use = "MAST",
  assay = "RNA"
)

# only keep markers with adjusted p-value beow 0.05
markers_mast %<>% filter(p_val_adj < 0.05)
saveRDS(markers_mast, file = "cluster_markers_MAST.rds")

markers_wilcox <- FindAllMarkers(integrated,
  only.pos = TRUE,
  assay = "RNA"
)

# only keep markers with adjusted p-value beow 0.05
markers_wilcox %<>% filter(p_val_adj < 0.05)
saveRDS(markers_wilcox, file = "cluster_markers_wilcox.rds")

# get intersection marker list of both methods
markers_intersection <- data.table::data.table()
for (i in levels(markers_wilcox$cluster)) {
  # calculate intersection cluster by cluster
  gene_list <- intersect(
    markers_wilcox$gene[markers_wilcox$cluster == i],
    markers_mast$gene[markers_mast$cluster == i]
  )
  # add selected genes to the table
  markers_intersection <- rbind(
    markers_intersection,
    arrange(
      filter(markers_wilcox, gene %in% gene_list & cluster == i),
      desc(avg_logFC)
    )
  )
}
rownames(markers_intersection) <- NULL
saveRDS(markers_intersection, file = "cluster_markers_mast_wilcox_intersection.rds")
