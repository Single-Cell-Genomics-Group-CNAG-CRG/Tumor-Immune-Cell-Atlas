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

atlas <- readRDS(file = paste(out_path, "Atlas_new_annotation.rds", sep = ''))
atlas$new_cell_types <- as.character(atlas$new_cell_types)
atlas <- subset(atlas, new_cell_types != "MAST")
atlas$new_cell_types <- as.factor(atlas$new_cell_types)

Idents(atlas) <- "new_cell_types"

markers_wilcox <- FindAllMarkers(atlas, only.pos = T, assay = "RNA")
markers_wilcox <- filter(markers_wilcox, p_val_adj < 0.05)
table(markers_wilcox$cluster)

saveRDS(markers_wilcox, paste0(out_path, "atlas_wilcox_markers_specific.rds"))