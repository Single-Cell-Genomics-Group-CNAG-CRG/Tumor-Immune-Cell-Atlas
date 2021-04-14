# This script plots the Local Inverse Simpson Index (LISI)


# Load packages
library(tidyverse)
library(ggridges)


# Load data frame
lisi_df <- readRDS("results/lisi_scores.rds")


# Plot LISI
sorted_corrections <- c("uncorrected", "Seurat v3", "Harmony", "Scanorama")
palette <- c("#999999", "#92e7df", "#612c63", "#e5624f")
lisi_scores_gg <- lisi_df %>%
  mutate(correction = factor(correction, levels = rev(sorted_corrections))) %>%
  ggplot(aes(source, correction, fill = correction)) +
    geom_density_ridges(
      quantile_lines = TRUE,
      quantile_fun = median
    ) +
    # geom_violin() +
    scale_fill_manual(values = palette) +
    labs(x = "iLISI", y = "") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(color = "black", size = 13),
      axis.text.y = element_text(color = "black", size = 12),
      axis.text.x = element_text(size = 11)
    )


# Save
ggsave(
  filename = "results/plots/iLISI_TICAtlas.pdf",
  plot = lisi_scores_gg,
  width = 14,
  height = 8,
  units = "cm"
)
