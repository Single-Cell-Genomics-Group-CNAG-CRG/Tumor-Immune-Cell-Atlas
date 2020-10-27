# this script creates the cell type proportion per
# cluster dataset that we will use for the clustering

# load libraries
library(tidyverse)
library(Seurat)

# create proportion table per patient
atlas <- readRDS("atlas.rds")
data <- table(atlas$patient, atlas$cell_type)
data <- as.data.frame.matrix(data)
data$total <- rowSums(data)

for (patient in rownames(data)) {
  # iterate over the cell types
  for (cell in colnames(data)[1:25]) {
    # calculate patient's cell type percentage
    data[patient, cell] <- round(data[patient, cell] / data[patient, "total"] * 100, 3)
  }
}

# we filter out patients with less than 500 total
# cells as this could introduce bias
data <- data[data$total > 500, ]
data <- data[, 1:25]

# get the metadata we are interested in
meta <- atlas@meta.data[, c("source", "gender", "patient", "age", "treatment", "response", "state", "cms", "msi", "subtype")]

# unify metadata to one row per patient
meta <- unique(meta)
rownames(meta) <- meta$patient

data <- merge(meta, data, by = "row.names")
data$Row.names <- NULL

# remove datasets with no patient information
# (they are labeled as having only one patient)
data <- data[data$source %notin% c("lung1", "lung2"), ]

saveRDS(data, "atlas_proportion_dataset.rds")
