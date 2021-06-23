# process each dataset, preparing them for integration
# this scripts takes as arguments the data path for the dataset
# and the abbreviated name
# it normalizes and finds HVG of the dataset

library(Seurat)
library(tidyverse)

args <- commandArgs(TRUE)
data <- args[1]
name <- args[2]

print(name)

# load dataset
data <- readRDS(data)

print(data)

# normalize and identify variable features

data <- data %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 5000)

print("done!")

saveRDS(data, paste0("output/", name, "_normalized.rds"))
