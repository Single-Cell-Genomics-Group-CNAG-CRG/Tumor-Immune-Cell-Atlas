# This script clusters the cells in the TICAtlas at varying resolutions using
# the Harmony-corrected PC as features.


# Load packages
library(tidyverse)
library(Seurat)


# Load Seurat object
tica <- readRDS("results/tica_with_harmony.rds")


# Cluster
tica <- FindClusters(tica, resolution = c(0.01, 0.05, 0.075, 0.1, 0.2))


# Save UMAP coordinates and clusters
harmony_df <- data.frame(
  cell_barcode = colnames(tica),
  UMAP1 = tica@reductions$umap@cell.embeddings[, "UMAP_1"],
  UMAP2 = tica@reductions$umap@cell.embeddings[, "UMAP_2"],
  cell_type = tica$cell_type,
  cancer_subtype = tica$subtype,
  cluster_res_0.01 = tica$RNA_snn_res.0.01,
  cluster_res_0.05 = tica$RNA_snn_res.0.05,
  cluster_res_0.075 = tica$RNA_snn_res.0.075,
  cluster_res_0.1 = tica$RNA_snn_res.0.1,
  cluster_res_0.2 = tica$RNA_snn_res.0.2
)
saveRDS(harmony_df, "results/umap_harmony_with_clusters.rds")
