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
source("myfunctions.R")
`%notin%` <- Negate(`%in%`)
values = as.vector(polychrome())

# Enable parallelization
options(future.globals.maxSize = +Inf)
plan("multiprocess", workers = 24)

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

atlas <- readRDS(file = paste(out_path, "Atlas_SCT.rds", sep = ''))
Idents(atlas) <- "new_cell_types"

markers_mast <- FindAllMarkers(atlas, only.pos = T, assay = "SCT", test.use = "MAST", logfc.threshold = 0, min.pct = 0)
markers_mast <- filter(markers_mast, p_val_adj < 0.05)
table(markers_mast$cluster)

saveRDS(markers_mast, paste0(out_path, "atlas_mast_markers_specific_SCT.rds"))