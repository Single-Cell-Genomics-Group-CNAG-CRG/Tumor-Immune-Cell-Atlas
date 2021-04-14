# This script computes harmony for the TICA atlas and saves the resulting
# batch-corrected PCs and Seurat object


# Load packages
library(Seurat)
library(tidyverse)
library(harmony)
library(SeuratWrappers)


# Load data
tica <- readRDS("data/TICAtlas.rds")


# Save Seurat-corrected PCA coordinates
saveRDS(tica@reductions$pca@cell.embeddings, "results/Seurat_v3_corrected_pca_coordinates.rds")


# Run Harmony
DefaultAssay(tica) <- "RNA"
tica <- tica %>%
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE, selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA()
print(tica@reductions)
tica <- RunHarmony(tica, dims = 1:50, reduction = "pca", group.by.vars = "source")
tica <- RunUMAP(tica, dims = 1:50, reduction = "harmony")
tica <- FindNeighbors(tica, dims = 1:50, reduction = "harmony")


# Save
saveRDS(tica@reductions$pca@cell.embeddings, "results/uncorrected_pca_coordinates.rds")
write_csv(as.data.frame(tica@reductions$pca@cell.embeddings), "results/uncorrected_pca_coordinates.csv", col_names = TRUE)
saveRDS(tica@reductions$harmony@cell.embeddings, "results/harmony_corrected_pca_coordinates.rds")
write_csv(rownames_to_column(tica@meta.data, var = "cell_barcode"), "results/tica_metadata.csv", col_names = TRUE)
saveRDS(tica, "results/tica_with_harmony.rds")


