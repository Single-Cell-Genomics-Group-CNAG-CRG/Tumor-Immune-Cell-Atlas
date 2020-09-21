#!/usr/bin/env Rscript

## Load libraries
library(SPOTlight)
library(dplyr)

## Set parameters
source("misc/paths.R")
dir.create(path = sprintf("%s/%s", an_breast_10x, robj_dir),
           showWarnings = FALSE,
           recursive = TRUE)

dir.create(path = sprintf("%s/%s", an_breast_10x, plt_dir),
           showWarnings = FALSE,
           recursive = TRUE)

cl_n <- 100
hvg <- 3000
ntop <- NULL
transf <- "uv"
method <- "nsNMF"
min_cont <- 0
clust_vr <- "new_cell_types"

### Variable passed from command line to determine the SC and spatial dataset  
args <- commandArgs(trailingOnly=TRUE)
# trn <- "melanoma"
# spatial <- "breast"
trn <- args[[1]]
spatial <- args[[2]]



if (is.null(ntop)) {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-NULL_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, transf, method, min_cont)
} else {
  spotlight_id <- sprintf("trn-%s_cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s_mincont-%s",
                          trn, cl_n, hvg, ntop, transf, method, min_cont)
}


## Load data
### Reference scRNAseq dataset
if (trn == "full") {
  print("Loading full 250 TICA")
  ica_se <- readRDS(file = "data/immune_reference/ica_se_full_processed.RDS")
} else {
  print(sprintf("Loading: data/immune_reference/ica_se_%s_processed.RDS",trn))
  ica_se <- readRDS(file = sprintf("data/immune_reference/ica_se_%s_processed.RDS",
                                   trn))
}

### Cell type markers obtained from Seruat::FindAllMarkers
if (trn == "full") {
  print("Loading full TICA markers")
  ica_markers <- readRDS(file = "data/immune_reference/ica_markers_all.RDS")
} else {
  print(sprintf("Loading: data/immune_reference/ica_markers_%s.RDS",
                trn))
  ica_markers <- readRDS(file = sprintf("data/immune_reference/ica_markers_%s.RDS",
                                        trn))
}

### Spatial data
if (spatial == "breast") {
  print("Loading breast spatial")
  spatial_se <- readRDS(file = sprintf("%s/%s/breast_merged_processed.RDS",
                                    an_breast_10x, robj_dir))
} else if (spatial == "hn") {
  print("Loading hn spatial")
  spatial_se <- readRDS(file = sprintf("%s/%s/hn1_processed.RDS",
                                 an_aussie, robj_dir))
} else if (spatial == "pdac") {
  print("Loading pdac spatial")
  
} else if (spatial == "prostate") {
  print("Loading prostate spatial")
  
}

# Run deconvolution
## Here we set transf = RAW because when we carried out SCTransform on the SC and spatial dataset we already performed the normalization.
## That is why we point directly to the SCT assay and data slot 
decon_mtrx_ls <- SPOTlight::spotlight_deconvolution(
  se_sc = ica_se,
  counts_spatial = spatial_se@assays$SCT@data,
  clust_vr = "specific_cell_type_mod",
  cluster_markers = ica_markers,
  cl_n = cl_n,
  hvg = hvg,
  ntop = ntop,
  transf = "raw", # Set to raw 
  method = method,
  min_cont = min_cont,
  assay = "SCT",
  slot = "data")

saveRDS(object = decon_mtrx_ls,
        file = sprintf("%s/%s/decon_mtrx_TICA_%s_%s.RDS",
                       an_breast_10x, robj_dir, spatial, spotlight_id))


