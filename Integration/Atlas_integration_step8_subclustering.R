### in this script we will integrate the small (subsampeled) immune datasets
### more datasets will be added

# set working directory
rm(list=ls())
setwd("/home/devel/pnieto/TIL_Atlas")

# load libraries
library(tidyverse)
library(Seurat)
library(sctransform)
library(ggsci)
library(pals)
library(RColorBrewer)
#source("myfunctions.R")

c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
         "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown") # distinct color palette

values = as.vector(polychrome())

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

integrated <- readRDS("output/Atlas/Atlas_integrated_clustered.rds")
print(integrated)

# reannotate integrated atlas ---------------------------------------------
print("reannotate integrated atlas")

CD4 <- readRDS("data/2020-04-02_subclustered_CD4.rds")
CD8 <- readRDS("data/2020-04-03_subclustered_CD8.rds")
Bcell <- readRDS("data/2020-04-03_subclustered_Bcell.rds")
Myeloid <- readRDS("data/2020-04-03_subclustered_Myeloid.rds")

# merge subclustered datasets
print("merge subclustered datasets")

merged <- merge(CD4, c(CD8, Bcell, Myeloid))
rm(CD4, CD8, Bcell, Myeloid)

# load integrated Atlas
print("load integrated Atlas")

integrated <- readRDS("output/Atlas/Atlas_integrated_clustered.rds")
integrated$new_assigned_cell_type <- "Other"
print(unique(integrated$new_assigned_cell_type))
integrated <- AddMetaData(integrated, merged$cell_type, "new_assigned_cell_type")
integrated <- subset(integrated, subset = new_assigned_cell_type != "Other" & new_assigned_cell_type != "Macro/DC")
Idents(integrated)  <- "new_assigned_cell_type"
DimPlot(integrated, label = T)

saveRDS(integrated, file = paste(out_path, date, "_Atlas_subclustered_annotated.rds", sep = ''))

