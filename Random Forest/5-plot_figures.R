# This script plots the results of the random forest


# Load packages
library(ggpubr)
library(pheatmap)
library(tidyverse)


# Load data
accuracy_df <- readRDS("tmp/accuracy_random_forest_dataframe.rds")
conf_mat_probs_l <- readRDS("tmp/confusion_matrix_probabilities_list.rds")


# Plot accuracy
accuracy_gg <- accuracy_df %>%
  ggplot(aes(type, accuracy, col = type)) +
  geom_jitter() +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(values = c("gray50", "limegreen")) +
  labs(x = "", y = "Accuracy") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 11, color = "black"))

kappa_gg <- accuracy_df %>%
  ggplot(aes(type, kappa, col = type)) +
  geom_jitter() +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(values = c("gray50", "limegreen")) +
  labs(x = "", y = "Kappa") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 11, color = "black"))

arranged <- ggarrange(plotlist = list(accuracy_gg, kappa_gg), ncol = 2)
ggsave(
  plot = arranged,
  filename = "results/accuracy_kappa_random_forest.pdf",
  width = 17.5,
  height = 8,
  units = "cm"
)


# Plot heatmap confusion matrix
colors_function <- colorRampPalette(colors = c("white", "red"))
colors <- colors_function(100)
heatmap <- pheatmap(
  conf_mat_probs_l$fold2,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "none",
  color = colors,
  angle_col = 315,
  legend = FALSE,
  fontsize_row = 9,
  fontsize_col = 9
)
save_pheatmap_pdf <- function(x, filename, width = 17.5 * 0.394, height = 15 * 0.394) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmap, "results/confusion_matrix_probabilities.pdf")

