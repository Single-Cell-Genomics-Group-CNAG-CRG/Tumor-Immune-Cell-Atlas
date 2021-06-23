# integrate the datasets

library(Seurat)

# load anchors
anchors <- readRDS("output/anchors.rds")

integrated <- IntegrateData(anchorset = anchors)

rm(anchors)
gc()

# save integrated dataset
saveRDS(integrated, "output/integrated.rds")