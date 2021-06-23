# integrated analysis on all cells!
  
library(Seurat)
library(tidyverse)
# library(future)
# plan("multiprocess")
# options(future.globals.maxSize = Inf)

atlas <- readRDS("output/integrated.rds")

DefaultAssay(atlas) <- "integrated"

atlas
  
# run the standard workflow on the integrated assay

atlas <- atlas %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30)

saveRDS(atlas, "output/integrated_processed.rds")