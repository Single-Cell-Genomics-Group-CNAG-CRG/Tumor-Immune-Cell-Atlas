# This script clusters the cells in the TICAtlas at varying resolutions using
# the Scanorama-corrected PC as features.


# Load packages
import numpy as np
import pandas as pd
import scanpy as sc
import anndata


# Read data
path_to_data = "results/TICAtlas_scanorama3.h5ad"
tica = sc.read_h5ad(path_to_data)



# Cluster at varying resolutions
for res in [0.25, 0.5, 0.75, 1]:
	sc.tl.leiden(tica, resolution = res)
	tica.obs["leiden_res{}".format(res)] = tica.obs["leiden"]



# Save UMAP coords and clusters to visualize in a downstream script
scanorama_df = {
	"cell_barcode": tica.obs.index,
	"UMAP1": tica.obsm["X_umap"][:, 0],
	"UMAP2": tica.obsm["X_umap"][:, 1],
	"cell_type": tica.obs["cell_type"],
	"cancer_subtype": tica.obs["subtype"],
	"leiden_res0.25": tica.obs["leiden_res0.25"],
	"leiden_res0.5": tica.obs["leiden_res0.5"],
	"leiden_res0.75": tica.obs["leiden_res0.75"],
	"leiden_res1": tica.obs["leiden_res1"]
}
scanorama_df = pd.DataFrame(scanorama_df)
scanorama_df.to_csv("results/umap_scanorama_with_clusters.csv", index = False)
