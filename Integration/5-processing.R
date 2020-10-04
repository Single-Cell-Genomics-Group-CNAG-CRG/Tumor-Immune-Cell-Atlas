# this scrippts starts the downstream processing of the integrated dataset (atlas)

# load libraries
library(tidyverse)
library(Seurat)
library(magrittr)

integrated <- readRDS("integrated_object.rds")
integrated %<>% ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

saveRDS(integrated, file = "integrated_object_processed.rds")
