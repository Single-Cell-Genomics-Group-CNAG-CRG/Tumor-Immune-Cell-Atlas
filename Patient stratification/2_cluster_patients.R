# this script calculates the 6 patient clusters
# using the hierarchical k-means algorithm

## load packages
library(tidyverse)
library(magrittr)
library(ggsci)
library(pals)
library(RColorBrewer)
library(patchwork)
library(clValid)
library(factoextra)
library(NbClust)
library(GGally)
library(dendextend)
library(Rtsne)
library(ComplexHeatmap)
library(circlize)
library(kohonen)
library(mclust)

df <- readRDS("output/atlas_proportion_dataset.rds")
tica_pal <- readRDS("output/tica_palettes.rds")

df$subtype[df$source=="liver2"] <- "HCC"

df <- mutate(df,
             `B cells` = `Activated B cells` + `Memory B cells` + `Naive B cells` + `Unswitched Memory B cells`,
              )

df[, c("Activated B cells", "Memory B cells", "Naive B cells", "Unswitched Memory B cells")] <- NULL

# optimal number of clusters (from 1 to 20)
clust_num <- NbClust(
  data = df[, 4:ncol(df)], distance = "euclidean", min.nc = 2,
  max.nc = 20, method = "kmeans", index = "alllong"
)
fviz_nbclust(clust_num)

# HIERACHICAL K-MEANS
hkmeans_cluster <- hkmeans(x = df[, 4:(ncol(df))], hc.metric = "euclidean", hc.method = "ward.D2", k = 6)
fviz_cluster(object = hkmeans_cluster, pallete = "jco", repel = TRUE)
df$cluster_kmeans_k6 <- hkmeans_cluster$cluster

# visualize as heatmap
df <- df[order(df$cluster_kmeans_k6), ]
mat <- as.matrix(df[, 4:(ncol(df) - 1)])
rownames(mat) <- NULL

row_ha <- rowAnnotation(
  cluster = df$cluster_kmeans_k6,
  col = list(cluster = c("1" = "#1F77B4FF", "2" = "#FF7F0EFF", "3" = "#2CA02CFF", "4" = "#D62728FF", "5" = "#9467BDFF", "6" = "#8C564BFF"))
)
Heatmap(mat,
        name = "composition (%)",
        cluster_rows = F, cluster_columns = FALSE,
        left_annotation = row_ha, column_names_max_height = unit(10, "cm"),
        col = colorRampPalette(brewer.pal(8, "Blues"))(50)
)

# visualize cancer types per cluster
ggplot(df, aes(factor(cluster_kmeans_k6))) +
  geom_bar(aes(fill = factor(subtype)), position = "fill") +
  coord_flip() +
  theme_classic() +
  labs(x = "Cluster", fill = "Cancer type", y = "Proportion of patients") +
  scale_fill_manual(values = tica_pal$cancer) +
  theme(text = element_text(size = 30))#, legend.position = "none")

ggplot(df, aes(factor(cluster_kmeans_k6))) +
  geom_bar(aes(fill = factor(subtype))) +
  coord_flip() +
  theme_classic() +
  labs(x = "Cluster", fill = "Cancer type", y = "Number of patients") +
  scale_fill_manual(values = tica_pal$cancer) +
  theme(text = element_text(size = 30))

# calculate PCA
res_pca <- prcomp(df[, 4:(ncol(df) - 1)], scale = F)
fviz_eig(res_pca)
res_ind <- get_pca_ind(res_pca)

# Contributions of variables to PC1 and PC2
p1 <- fviz_contrib(res_pca, choice = "var", axes = 1, top = 10, xtickslab.rt = 0) +
  theme(text = element_text(size = 35), title = element_text(size = 25), axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90)) +
  labs(title = "PC1")
p2 <- fviz_contrib(res_pca, choice = "var", axes = 2, top = 10, xtickslab.rt = 0) +
  theme(text = element_text(size = 35), title = element_text(size = 25), axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90)) +
  labs(title = "PC2")
p1 + p2

# PC contributions table
var <- get_pca_var(res_pca)
ft <- as.data.frame(var$contrib[, 1:2])
ft$Variable <- rownames(ft)
colnames(ft) <- c("PC1 (%)", "PC2 (%)", "Cell type")
ft <- ft[, c("Cell type", "PC1 (%)", "PC2 (%)")]
ft <- arrange(ft, desc(`PC1 (%)`))
ft %>%
  flextable::flextable() %>%
  flextable::theme_vanilla() %>%
  flextable::autofit()

# visualization: TSNE PLOT
tsne <- Rtsne(df[, 4:(ncol(df) - 1)], dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)

df$tsne_x <- tsne$Y[, 1]
df$tsne_y <- tsne$Y[, 2]

p1 <- ggplot(df, aes(tsne_x, tsne_y, color = as.factor(cluster_kmeans_k6))) +
  geom_point(size = 6) +
  labs(color = "cluster", x = "tSNE_1", y = "tSNE_2") +
  scale_color_d3() +
  theme_classic() +
  theme(
    text = element_text(size = 30),
    legend.position = "none"
  ) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Cluster")

p2 <- ggplot(df, aes(tsne_x, tsne_y, color = as.factor(subtype))) +
  geom_point(size = 6) +
  labs(color = "cluster", x = "tSNE_1", y = "tSNE_2") +
  theme_classic() +
  theme(
    text = element_text(size = 30),
    legend.position = "none"
  ) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Cancer type") +
  scale_color_manual(values = tica_pal$cancer)

p1 + p2

saveRDS(df, "output/atlas_proportion_dataset_clustered.rds")
