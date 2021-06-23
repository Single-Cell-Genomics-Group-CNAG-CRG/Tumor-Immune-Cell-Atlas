# this script creates the cell type proportion per
# cluster dataset that we will use for the clustering

# load libraries
library(tidyverse)
library(Seurat)

# create proportion table per patient
atlas <- readRDS("../output/integrated_renamed_filtered.rds")
data <- table(atlas$patient, atlas$new_annot)
data <- as.data.frame.matrix(data)
data$total <- rowSums(data)

# we filter out patients with less than 500 total
# cells as this could introduce bias
data <- filter(data, total > 500)

for (patient in rownames(data)) {
  # iterate over the cell types
  for (cell in colnames(data)[1:(ncol(data)-1)]) {
    # calculate patient's cell type percentage
    data[patient, cell] <- round(data[patient, cell] / data[patient, "total"] * 100, 3)
  }
}

data <- data[, 1:(ncol(data)-1)]

# get the metadata we are interested in
meta <- atlas@meta.data[, c("source", "patient", "subtype")]

# unify metadata to one row per patient
meta <- unique(meta)
rownames(meta) <- meta$patient

data <- merge(meta, data, by = "row.names")
data$Row.names <- NULL

# remove datasets with no patient information
# (they are labeled as having only one patient)
data <- filter(data, !source %in% c("lung1", "lung2"))

saveRDS(data, ".../output/atlas_proportion_dataset.rds")
