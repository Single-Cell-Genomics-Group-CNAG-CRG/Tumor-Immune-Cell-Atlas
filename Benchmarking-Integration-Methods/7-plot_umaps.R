# This script plots the UMAP for Harmony and Scanorama integration of the
# TICAtlas


# Load packages
library(Seurat)
library(tidyverse)


# Read data
umap_harmony_df <- readRDS("results/umap_harmony_with_clusters.rds")
umap_scanorama_df <- read_csv("results/umap_scanorama_with_clusters.csv")
umap_scanorama_df <- as.data.frame(umap_scanorama_df)
rownames(umap_scanorama_df) <- umap_scanorama_df$cell_barcode
umap_scanorama_df <- umap_scanorama_df[rownames(umap_harmony_df), ]
all(rownames(umap_harmony_df) == rownames(umap_scanorama_df))
umap_scanorama_df$cell_type <- umap_harmony_df$cell_type


# UMAP colored by cell type
cell_type_palette <- readRDS("data/cell_type_palette.rds")
dfs <- list(harmony = umap_harmony_df, scanorama = umap_scanorama_df)
umaps_cell_type <- purrr::map(dfs, function(df) {
  p <- df %>%
    ggplot(aes(UMAP1, UMAP2, color = cell_type)) +
      geom_point(size = 0.1) +
      scale_color_manual(values = cell_type_palette) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  p
})


# UMAP colored by cancer subtype
cancer_subtype_palette <- readRDS("data/cancer_subtype_palette.rds")
cancer_subtype_palette <- c(cancer_subtype_palette, CM = "#ee6a50")
umaps_cancer_subtype <- purrr::map(dfs, function(df) {
  p <- df %>%
    ggplot(aes(UMAP1, UMAP2, color = cancer_subtype)) +
      geom_point(size = 0.1) +
      scale_color_manual(values = cancer_subtype_palette) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  p
})


# UMAPs colored by cluster
# Harmony
cluster_vars <- str_subset(colnames(dfs$harmony), "cluster_res_")
umap_cluster_harmony <- purrr::map(cluster_vars, function(x) {
  p <- umap_harmony_df %>%
    ggplot(aes_string("UMAP1", "UMAP2", color = x)) +
      geom_point(size = 0.1) +
      ggtitle(str_remove(x, "cluster_res_")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  p
})

fct_cells_cell_type <- umap_harmony_df %>%
  select("cell_type", "cluster_res_0.05") %>%
  group_by(cell_type, cluster_res_0.05) %>%
  summarize(n_cells = n()) %>%
  ungroup() %>% 
  group_by(cluster_res_0.05) %>%
  mutate(n_cells_total = sum(n_cells), fct_cells = n_cells / n_cells_total)

fct_cells_cell_type_gg <- fct_cells_cell_type %>%
  ggplot(aes(cluster_res_0.05, fct_cells, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = cell_type_palette) +
  labs(x = "Cluster (Harmony)", y = "", fill = "") +
  theme_classic()


# Scanorama


# Save
for (x in names(dfs)) {
  # Save umap colored by cell type
  path_save_umap1 <- str_c(
    "results/plots/TICAtlas_umap_colored_cell_type_",
    x,
    ".pdf",
    sep = ""
  )
  ggsave(
    path_save_umap1,
    plot = umaps_cell_type[[x]],
    width = 9,
    height = 9,
    units = "cm"
  )
  
  
  # Save umap colored by cell type
  path_save_umap2 <- str_c(
    "results/plots/TICAtlas_umap_colored_cancer_subtype_",
    x,
    ".pdf",
    sep = ""
  )
  ggsave(
    path_save_umap2,
    plot = umaps_cancer_subtype[[x]],
    width = 9,
    height = 9,
    units = "cm"
  )
}


