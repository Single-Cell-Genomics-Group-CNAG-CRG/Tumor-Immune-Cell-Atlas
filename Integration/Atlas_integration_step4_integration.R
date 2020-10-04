rm(list=ls())
setwd("/home/devel/pnieto/TIL_Atlas")

# load libraries
library(tidyverse)
library(Seurat)
library(sctransform)
#library(future)

#options(future.globals.maxSize = 26826768384)
#plan("multiprocess", workers = 24)

out_path <- "output/"
plot_path <- "plots/"
date <- Sys.Date()

anchors <- readRDS("output/Atlas/anchors.rds")
integrated <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize", verbose = T, dims = 1:20)
rm(anchors)
saveRDS(integrated, file = "output/Atlas/integrated_object.rds")