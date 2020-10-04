rm(list=ls())
setwd("/home/devel/pnieto/TIL_Atlas")

# load libraries
library(tidyverse)
library(Seurat)
library(sctransform)
#source("myfunctions.R")
#library(future)

#options(future.globals.maxSize = 16826768384)
#plan("multiprocess", workers = 24)

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

features <- readRDS("output/Atlas/filtered_features.rds")
obj_list <- readRDS("output/Atlas/obj_list_PrepSCT.rds")
print(features)
print(obj_list)

anchors <- FindIntegrationAnchors(obj_list, normalization.method = "LogNormalize", anchor.features = features, reduction = "rpca", verbose = T, dims = 1:20)
rm(features)
saveRDS(anchors, file = "output/Atlas/anchors.rds")