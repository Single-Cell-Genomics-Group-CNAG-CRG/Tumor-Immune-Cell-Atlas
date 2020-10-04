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

options(future.globals.maxSize = 26826768384)
plan("multiprocess", workers = 24)

`%notin%` <- Negate(`%in%`)

c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
         "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown") # distinct color palette

out_path <- "output/Atlas/"
plot_path <- "plots/Atlas/"
date <- Sys.Date()

# load data (only immune cells, remove CAFs!) -----
obj_list <- character()
data_list <- c("breast", "colorectal", "liver_10x", "liver_smart", "livercancer", "lung", "lung2", "mcolon", "melanoma", "melanoma2", "mendo", "mlung", "mrenal", "ovary", "pancreas", "skin_SCC", "skin_BCC", "uveal_melanoma")
for (item in data_list){
  print(item)
  temp <- readRDS(paste("data/", item, "_all.rds", sep = ''))
  print(temp)	
  #temp <- subset(temp, subset = cell_type_og %notin% c("Other", "Others", "CAF", "Fibro", "tumor", "CAFs", "Alveolar", "EC", "Epi"))
  #print(unique(temp$cell_type_og))
  print(table(temp$cell_type_og))
  #rownames(temp@assays$RNA) <- make.names(rownames(temp@assays$RNA), unique = T)
  DefaultAssay(temp) <- "RNA"
  #temp[["SCT"]] <- NULL
  obj_list <- c(obj_list, temp)
  print("")
}

# prepare integration
saveRDS(obj_list, paste0("/home/devel/pnieto/TIL_Atlas/output/Atlas/", date, "_obj_list_RAW.rds"))

#obj_list <- sapply(obj_list, function(x) {
	x <- NormalizeData(x, verbose = FALSE)
   	x <- FindVariableFeatures(x, nfeatures = 2000, verbose = FALSE, selection.method = "vst")
})

#obj_list <- readRDS(paste0("/home/devel/pnieto/TIL_Atlas/output/Atlas/obj_list_SCT.rds"))
#saveRDS(obj_list, "/home/devel/pnieto/TIL_Atlas/output/Atlas/obj_list_SCT.rds")
saveRDS(merge(obj_list[[1]], obj_list[2:length(obj_list)]), paste0("/home/devel/pnieto/TIL_Atlas/output/Atlas/", date, "_obj_list_merged.rds"))
#saveRDS(obj_list, "/home/devel/pnieto/TIL_Atlas/output/Atlas/obj_list_SCT_noPCA.rds")