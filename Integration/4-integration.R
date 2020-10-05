# this script performs the integration of the datasets
# based on the previously calculated anchors

# load libraries
library(tidyverse)
library(Seurat)
library(sctransform)

anchors <- readRDS("anchors.rds")
integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "LogNormalize",
  verbose = TRUE,
  dims = 1:20
)
rm(anchors)
saveRDS(integrated, file = "integrated_object.rds")
