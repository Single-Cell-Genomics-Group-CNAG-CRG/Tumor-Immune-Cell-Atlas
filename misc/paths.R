# Variables
ver <- "2020-06-04"

## PATHS and common variables
version_dir <- sprintf("%s",ver)
dir.create(path = version_dir, 
           showWarnings = F, 
           recursive = T)

## Create directory for plots
plt_dir <- sprintf("%s/plots_%s", ver, ver) 
# dir.create(path = plt_dir, 
#            showWarnings = FALSE, 
#            recursive = TRUE)

## Create directory for RDS objects
robj_dir <- sprintf("%s/R_objects_%s", ver, ver)
# dir.create(path = robj_dir, 
#            showWarnings = FALSE, 
#            recursive = TRUE)

## Paths to all the folders
an <- "analysis"
an_aussie <- sprintf("%s/aussie_ln", an)
an_oro <- sprintf("%s/aussie_oro", an)
an_melanoma <- sprintf("%s/melanoma", an)
an_prostate <- sprintf("%s/prostate", an)
an_breast_10x <- sprintf("%s/breast_10x", an)
an_pdac <- sprintf("%s/pdac_st", an)
an_epid <- sprintf("%s/epid_27", an)
