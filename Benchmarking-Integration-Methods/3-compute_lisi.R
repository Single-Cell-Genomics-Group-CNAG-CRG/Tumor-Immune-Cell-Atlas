# This script computes the Local Inverse Simpson Index for the following cases:
# 1. Unintegrated PCA
# 2. Seurat v3-corrected PCA
# 3. Harmony-corrected PCA
# 4. Scanorama-corrected PCA


# Load packages
library(tidyverse)
library(lisi)


# Load data
tica_metadata <- read_csv("results/tica_metadata.csv", col_names = TRUE)
uncorrected_coords <- readRDS("results/uncorrected_pca_coordinates.rds")
seurat_coords <- readRDS("results/Seurat_v3_corrected_pca_coordinates.rds")
harmony_coords <- readRDS("results/harmony_corrected_pca_coordinates.rds")
scanorama_coords <- read_csv("results/Scanorama_corrected_pca_coordinates3.csv")
selected_cols <- str_subset(colnames(scanorama_coords), "^Scanorama_")
scanorama_coords_mat <- scanorama_coords[, selected_cols]
scanorama_coords_mat <- as.matrix(scanorama_coords_mat)
rownames(scanorama_coords_mat) <- scanorama_coords$cell_barcode
scanorama_coords_mat <- scanorama_coords_mat[rownames(harmony_coords), ]


# Define technical and biological variables
meta_data <- data.frame(source = tica_metadata$source)
rownames(meta_data) <- tica_metadata$cell_barcode
if (all(rownames(meta_data) == rownames(uncorrected_coords))) {
  print("row names are equal")
} else {
  warning("row names are not equal!")
}


# Compute LISI
dim_red_mats <- list(
  uncorrected_coords,
  seurat_coords,
  harmony_coords,
  scanorama_coords_mat
)
names(dim_red_mats) <- c("uncorrected", "Seurat v3", "Harmony",
                         "Scanorama")
lisi_scores <- purrr::map(dim_red_mats, function(mat) {
  scores <- compute_lisi(
    X = mat,
    meta_data = meta_data,
    label_colnames = "source"
  )
})
lisi_scores <- bind_rows(lisi_scores, .id = "correction")


# Save
saveRDS(lisi_scores, "results/lisi_scores.rds")
