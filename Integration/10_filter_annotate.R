# filter out non-immune clusters (clusters 15, 16, 20 and 34) and rename clusters (annotate them)

library(Seurat)
library(tidyverse)
# library(future)
# plan("multicore")
# options(future.globals.maxSize = Inf)

atlas <- readRDS("output/integrated_clustered.rds")

# remove non-immune cells
atlas <- subset(atlas, seurat_clusters %in% c("15", "16", "20", "34"), invert = TRUE)

atlas$seurat_clusters <- as.factor(as.numeric(as.character(atlas$seurat_clusters)))
Idents(atlas) <- "seurat_clusters"
levels(atlas)

new_ids <- c("CD4 memory stem cells", "CD8 exhausted", "CD8 cytotoxic", "CD4 memory stem cells", "CD4 resident effector memory", "T regs", "B cells", "TAMs C1QC", "NK", "CD4 effector memory RA",
             "Tumor infiltrating Monocytes", "CD4 activated", "Proliferative T cells", "CD4 transitional memory", "CD4 effector memory RA", "DC2 CD1C+", "CD8 IFN activated", "CD8 activated T cells",
             "Plasma B cells", "Macrophages SPP1", "Plasma B cells", "DC4 CD1C-", "T helper/Th17", "pDC", "DC3 LAMP3", "Macrophages proliferative", "Macrophages CXCL10",
             "Mast cells", "DC1", "Proliferative B cells", "Plasma blasts"
             )

names(new_ids) <- levels(atlas)

atlas <- RenameIdents(atlas, new_ids)

atlas$new_annot <- Idents(atlas)

ggsave(plot = DimPlot(atlas, cols = as.vector(pals::polychrome())), filename = "atlas_plot.png", width = 10)

saveRDS(atlas, "output/integrated_renamed_filtered.rds")