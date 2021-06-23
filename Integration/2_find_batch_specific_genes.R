# to improve the integration, we limit the genes used to compute the anchors later on
# we will remove the top genes that are specific for each batch
# in our case the batch is in the "source" variable

library(Seurat)
library(tidyverse)

source("code/functions/load_objects.R")

obj_list <- load_objects(directory="output/",
                         extension="_filtered_genes.rds")

merged <- merge(obj_list[[1]],
                obj_list[2:length(obj_list)])

# to keep the consistency we are going to rename the cells

merged <- RenameCells(merged,
                      new.names =paste("cell", 1:ncol(merged),
                                       sep = "_"))

merged <- merged %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData()

# set the batches (source) as idents
Idents(merged) <- "source"

# get markers for each batch, limting cells per ident to reduce computing time
markers <- FindAllMarkers(merged,
                          verbose = TRUE,
                          only.pos = TRUE,
                          max.cells.per.ident = 5000)

# filter by adjusted p_value
markers <- filter(markers, 
                  p_val_adj < 0.01)

saveRDS(markers,
        file = "output/markers/batch_specific_markers.rds")
