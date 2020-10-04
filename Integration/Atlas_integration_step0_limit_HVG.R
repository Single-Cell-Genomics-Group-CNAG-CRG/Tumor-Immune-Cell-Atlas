# set working directory
rm(list=ls())
setwd("/home/devel/pnieto/TIL_Atlas")
#setwd("C:/Users/PNIETO/OneDrive - CRG - Centre de Regulacio Genomica/projects/TIL Atlas")

# load libraries
library(tidyverse)
library(Seurat)
library(sctransform)
library(matchSCore2)
library(magrittr)
library(future)

options(future.globals.maxSize = +Inf)
plan("multiprocess", workers = 24)

`%notin%` <- Negate(`%in%`)

c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
         "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown") # distinct color palette

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

#obj_list <- readRDS("/home/devel/pnieto/TIL_Atlas/output/Atlas/2020-07-06_obj_list_RAW.rds")

#obj1 <- merge(obj_list[[1]], obj_list[2:6], merge.data = F)
#obj2 <- merge(obj_list[[7]], obj_list[8:12], merge.data = F)
#obj3 <- merge(obj_list[[13]], obj_list[14:length(obj_list)], merge.data = F)
#rm(obj_list)
#obj <- merge(obj1, c(obj2, obj3), merge.data = F)
#rm(obj1, obj2, obj3)

#saveRDS(obj, "/home/devel/pnieto/TIL_Atlas/output/Atlas/obj_list_merged.rds")

obj <- readRDS("/home/devel/pnieto/TIL_Atlas/output/Atlas/obj_list_merged.rds")
Idents(obj) <- "source"

markers <- FindAllMarkers(obj, verbose = T, only.pos = T, max.cells.per.ident = 5000)
print(head(markers))
saveRDS(markers, file = paste(out_path, date, "_Atlas_markers_source.rds", sep = ''))