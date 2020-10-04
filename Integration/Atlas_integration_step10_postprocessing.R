# set working directory
rm(list=ls())
setwd("/home/devel/pnieto/TIL_Atlas")

# load libraries
library(tidyverse)
library(Seurat)
library(future)
library(sctransform)
library(ggsci)
library(pals)
library(RColorBrewer)
#source("myfunctions.R")
`%notin%` <- Negate(`%in%`)
values = as.vector(polychrome())

# Enable parallelization
options(future.globals.maxSize = +Inf)
plan("multiprocess", workers = 24)

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

atlas <- readRDS(file = paste(out_path, "Atlas_new_annotation.rds", sep = ''))

atlas$cluster_kmeans_k6<- as.character(atlas$cluster_kmeans_k6) 
atlas <- subset(atlas, cluster_kmeans_k6 %in% c("1", "2", "3", "4", "5", "6"))
atlas$cluster_kmeans_k6 <- as.factor(atlas$cluster_kmeans_k6)
table(atlas$cluster_kmeans_k6)

Idents(atlas) <- "cluster_kmeans_k6"

atlas <- subset(atlas, downsample = 10000)
table(atlas$cluster_kmeans_k6)

# kmeans clusters markers
markers_clusters <- readRDS(paste0(out_path, "atlas_wilcox_markers_k6_clusters.rds"))
markers_clusters <-  top_n(group_by(markers_clusters, cluster), n = 10, wt = avg_logFC)
table(markers_clusters$cluster)
markers_clusters <- split(markers_clusters$gene, markers_clusters$cluster)
length(unique(unlist(markers_clusters)))

atlas <- ScaleData(atlas, assay = "RNA", features = unlist(markers_clusters))

# heatmap with top 10 genes per cell type by logFC
DoHeatmap(atlas, features = unique(unlist(markers_clusters)), assay = "RNA", group.colors = pal_d3("category20")(20), size = 2, group.bar.height = 0.007) + 
  scale_fill_gradientn(colors = c("#6D9EC1", "white", "#E46726")) +
  guides(color = FALSE) + 
  ggtitle("Top 10 raw cluster markers")
ggsave(paste0(plot_path, "cluster_markers_heatmap_raw.pdf"), dpi = "retina", width = 15, height = 22)

markers_clusters <- readRDS(paste0(out_path, "atlas_wilcox_markers_k6_clusters.rds"))
markers_clusters <- markers_clusters[grep('^MT.', markers_clusters$gene, invert = T),]
markers_clusters <- markers_clusters[grep('^RP', markers_clusters$gene, invert = T),]
markers_clusters <- markers_clusters[grep('^HSP', markers_clusters$gene, invert = T),]
markers_clusters <- markers_clusters[grep('^ENS', markers_clusters$gene, invert = T),]
markers_clusters <-  top_n(group_by(markers_clusters, cluster), n = 10, wt = avg_logFC)
table(markers_clusters$cluster)
markers_clusters <- split(markers_clusters$gene, markers_clusters$cluster)
length(unique(unlist(markers_clusters)))

atlas <- ScaleData(atlas, assay = "RNA", features = unlist(markers_clusters))

# heatmap with top 10 genes per cell type by logFC FILTERED
DoHeatmap(atlas, features = unique(unlist(markers_clusters)), assay = "RNA", group.colors = pal_d3("category20")(20), size = 2, group.bar.height = 0.007) + 
  scale_fill_gradientn(colors = c("#6D9EC1", "white", "#E46726")) +
  guides(color = FALSE) + 
  ggtitle("Top 10 filtered cluster markers")
ggsave(paste0(plot_path, "cluster_markers_heatmap_filtered.pdf"), dpi = "retina", width = 15, height = 22)


q()

atlas <- SCTransform(atlas, vars.to.regress = "source", assay = "RNA", do.scale = TRUE, do.center = FALSE, conserve.memory = T)

saveRDS(atlas, file = paste(out_path, "Atlas_SCT.rds", sep = ''))

q()

atlas <- readRDS(file = paste(out_path, "Atlas_new_annotation.rds", sep = ''))

atlas$cluster_kmeans_k6 <- as.character(atlas$cluster_kmeans_k6)
atlas <- subset(atlas, cluster_kmeans_k6 %in% c("1", "2", "3", "4", "5", "6"))
DimPlot(atlas, split.by = "cluster_kmeans_k6", ncol = 2, group.by = "new_cell_types", cols = grDevices::adjustcolor(values, alpha.f = 0.5)) + 
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
ggsave(paste0(out_path, "Atlas_cell_types_split_k6_cluster.pdf"), width = 12, dpi = "retina")

rm(atlas)

atlas <- readRDS(file = paste(out_path, "Atlas_new_annotation.rds", sep = ''))
atlas$general_cell_type <- as.character(atlas$general_cell_type)
atlas <- subset(atlas, general_cell_type == "T cell" & specific_cell_type != "T cell")
atlas <- RunUMAP(atlas, dims = 1:20)

saveRDS(atlas, file = paste(out_path, "Atlas_new_annotation_tcells.rds", sep = ''))

DimPlot(atlas, group.by = "new_cell_types", cols = values)  + 
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1)) 
ggsave(paste0("Atlas/plots/", date, "full_atlas_specific_cell_types.pdf"), width = 10, dpi = "retina")

q()

DefaultAssay(atlas) <- "RNA"
atlas$new_cell_types <- as.factor(atlas$new_cell_types)
Idents(atlas) <- "new_cell_types"

markers_wilcox <- FindAllMarkers(atlas, only.pos = T, assay = "RNA", max.cells.per.ident = 5000)
markers_wilcox <- filter(markers_wilcox, p_val_adj < 0.05)
table(markers_wilcox$cluster)

saveRDS(markers_wilcox, paste0(out_path, "atlas_wilcox_markers.rds"))

markers_mast <- FindAllMarkers(atlas, only.pos = T, assay = "RNA", test.use = "MAST", max.cells.per.ident = 5000)
saveRDS(markers_mast, paste0(out_path, "atlas_MAST_markers.rds"))

q()

#saveRDS(atlas, file = paste(out_path, "Atlas_complete_annotation_renamed.rds", sep = ''))

print("marker intersection")
markers <- data.table::data.table()
for (i in levels(markers_wilcox$cluster)){
  print(i)
  gene_list <- intersect(markers_wilcox$gene[markers_wilcox$cluster == i], markers_mast$gene[markers_mast$cluster == i])
  #print(gene_list)
  #xlsx::write.xlsx(arrange(filter(markers_wilcox, gene %in% gene_list & cluster == i), desc(avg_logFC)), file="Atlas/output/tcell_markers_table.xlsx", sheetName = make.names(i), col.names=T, row.names=T, append = T)
  markers <- rbind(markers, arrange(filter(markers_wilcox, gene %in% gene_list & cluster == i), desc(avg_logFC)))
}

rownames(markers) <- NULL
saveRDS(markers, paste0(out_path, "atlas_MAST_wilcox_markers.rds"))


q()

Idents(atlas) <- "general_cell_type"
saveRDS(subset(atlas, downsample = 200), paste0(out_path, "atlas_200_general_cell_type_SCT.rds"))
Idents(atlas) <- "specific_cell_type"
saveRDS(subset(atlas, downsample = 250), paste0(out_path, "atlas_250_specific_cell_type_SCT.rds"))

atlas <- subset(atlas, downsample = 1000)
atlas <- Seurat::SCTransform(atlas, vars.to.regress = "source", do.scale = T, do.center = F)

###

#Idents(atlas) <- "general_cell_type"
#atlas_general_markers <- FindAllMarkers(atlas, only.pos = T, assay = "RNA", test.use = "MAST")
#print(table(atlas_general_markers$cluster))
#atlas_general_markers <- atlas_general_markers[grep('^MT.', atlas_general_markers$gene, invert = T),]
#atlas_general_markers <- atlas_general_markers[grep('^RP', atlas_general_markers$gene, invert = T),]
#atlas_general_markers %<>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
#saveRDS(atlas_general_markers, paste0(out_path, "top200_general_markers_logFC.rds"))

Idents(atlas) <- "specific_cell_type"
atlas_specific_markers <- FindAllMarkers(atlas, only.pos = T, assay = "SCT", logfc.threshold = 0, min.pct = 0, max.cells.per.ident = 5000)
print(table(atlas_specific_markers$cluster))
saveRDS(atlas_specific_markers, paste0(out_path, "raw_specific_markers_logFC.rds"))
atlas_specific_markers <- atlas_specific_markers[grep('^MT.', atlas_specific_markers$gene, invert = T),]
atlas_specific_markers <- atlas_specific_markers[grep('^RP', atlas_specific_markers$gene, invert = T),]
atlas_specific_markers <- atlas_specific_markers[grep('^ERCC', atlas_specific_markers$gene, invert = T),]
print(table(atlas_specific_markers$cluster))
atlas_specific_markers %<>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
saveRDS(atlas_specific_markers, paste0(out_path, "top200_specific_markers_logFC_SCT.rds"))

q()

atlas$general_cell_type <- as.character(atlas$general_cell_type)
atlas <- subset(atlas, general_cell_type == "T cell")
table(atlas$specific_cell_type)
Idents(atlas) <- "specific_cell_type"
atlas$general_cell_type <- as.factor(atlas$specific_cell_type)
atlas_tcell_markers <- FindAllMarkers(atlas, only.pos = T, assay = "RNA", test.use = "MAST", logfc.threshold = 0, min.pct = 0, latent.vars = "source")
print(table(atlas_tcell_markers$cluster))
atlas_tcell_markers %<>% group_by(cluster) %>% top_n(n = 250, wt = avg_logFC)
saveRDS(atlas_tcell_markers, paste0(out_path, "top250_tcell_markers_logFC.rds"))

#atlas <- subset(atlas, general_cell_type == "T cell")
#atlas <- RunUMAP(atlas, dims = 1:20)

#saveRDS(atlas, paste0(out_path, "atlas_Tcells.rds"))
#DimPlot(atlas, group.by = "specific_cell_type", cols = values) + guides(color = guide_legend(override.aes = list(size=3), ncol = 1))
#ggsave(paste0(out_path, "Tcell_UMAP_specific_cell_types.png"))
