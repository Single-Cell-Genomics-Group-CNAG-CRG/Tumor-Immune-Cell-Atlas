#!/usr/bin/env Rscript

# Train Prostate model
library(SPOTlight)
library(NMF)
library(Seurat)
library(dplyr)

# Set parameters
# trn <- "melanoma"
cl_n <- 100
hvg <- 3000
ntop <- NULL
transf <- "uv"
method <- "nsNMF"
min_cont <- 0
clust_vr <- "new_cell_types"

args <- commandArgs(trailingOnly=TRUE)
# trn <- "melanoma"
# spatial <- "161432"
trn <- args[[1]]
spatial <- args[[2]]

if (is.null(ntop)) {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-NULL_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, transf, method, min_cont)
} else {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, ntop, transf, method, min_cont)
}

source("misc/paths.R")
source("utils/bin.r")
source("utils/spatial_plot_spaniel.R")

# Load data
### Spatial data
st_ls <- readRDS(file = sprintf("%s/%s/processed_st_ls_oropharyngeal.RDS",
                                an_oro, robj_dir))

st_se <- st_ls[[spatial]]

### Immune reference atlas
if (trn == "full") {
  ica_se <- readRDS("data/immune_reference/atlas_250_specific_cell_type.rds")
} else if (trn == "breast") {
  ica_se <- readRDS("data/immune_reference/atlas_complete_annotation_breast.rds")
} else if (trn == "melanoma") {
  ica_se <- readRDS("data/immune_reference/atlas_complete_annotation_melanoma.rds")
}  else if (trn == "pancreas") {
  ica_se <- readRDS("data/immune_reference/atlas_complete_annotation_pancreas.rds")
}

# ica_se <- ica_se[, ica_se$new_cell_types != "MAST"]
ica_se[["specific_cell_type_mod"]] <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                           x = as.character(ica_se@meta.data[, clust_vr]),
                                           perl = TRUE)

### Markers TICA
if (trn == "full") {
  ica_markers <- readRDS(file = "data/immune_reference/ica_markers_all.RDS")
} else if (trn == "breast") {
  ica_markers <- readRDS(file = "data/immune_reference/ica_markers_breast.RDS")
} else if (trn == "melanoma") {
  ica_markers <- readRDS(file = "data/immune_reference/ica_markers_melanoma.RDS")
} else if (trn == "pancreas") {
  ica_markers <- readRDS(file = "data/immune_reference/ica_markers_melanoma.RDS")
}

# ica_markers <- ica_markers %>% filter(cluster != "MAST")

deconv_ls <- SPOTlight::spotlight_deconvolution(se_sc = ica_se,
                                                counts_spatial = st_se@assays$Spatial@counts,
                                                clust_vr = "specific_cell_type_mod",
                                                cluster_markers = ica_markers,
                                                cl_n = cl_n,
                                                hvg = hvg,
                                                ntop = ntop,
                                                transf = transf,
                                                method = method,
                                                min_cont = min_cont,
                                                assay = "RNA",
                                                slot = "counts")

saveRDS(object = deconv_ls,
        file = sprintf("%s/%s/spotlight_deconv_ls_%s_%s.RDS",
                       an_oro, robj_dir, spatial, spotlight_id))
