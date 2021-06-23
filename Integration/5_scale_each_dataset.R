# scale each dataset and compute its PCs
# using the selected integration features

library(Seurat)
library(tidyverse)

args = commandArgs(TRUE)
data = args[1]
name = args[2]

print(name)

# load dataset

data <- readRDS(data)

# load selected integration features

features <- readRDS("output/integration_features.rds")
print(length(features))


# normalize and identify variable features 

data <- data %>% 
  ScaleData(features = features) %>% 
  RunPCA(features = features)

print("done!")

saveRDS(data, paste0("output/", name, "_scaled.rds"))