rm(list=ls())
setwd("/home/devel/pnieto/TIL_Atlas")

# load libraries
library(tidyverse)
library(Seurat)
library(sctransform)
library(future)

options(future.globals.maxSize = +Inf)
plan("multiprocess", workers = 24)

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

obj_list <- readRDS("output/Atlas/obj_list_SCT.rds")
#obj_list <- readRDS("output/Atlas/obj_list_SCT_noPCA.rds")


print("features")
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000, verbose = TRUE)
print(length(features))
saveRDS(features, file = "output/Atlas/features.rds")

# substract source-specific features
source_genes <- readRDS("output/Atlas/2020-07-06_Atlas_markers_source.rds")
print(length(source_genes))
features = setdiff(features, source_genes)
print(length(features))
saveRDS(features, file = "output/Atlas/filtered_features.rds")

#print("prepSCT")
#obj_list <- PrepSCTIntegration(obj_list, features, verbose = T, assay = "SCT")

print("PCA + scale")
#obj_list <- lapply(X = obj_list, FUN = RunPCA, verbose = FALSE, features = features)

obj_list <- lapply(X = obj_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE, npcs = 30)
})

saveRDS(obj_list, file = "output/Atlas/obj_list_PrepSCT.rds")
print("Done! :)")