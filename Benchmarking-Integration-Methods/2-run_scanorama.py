# This script runs scanorama to obtain a batch-corrected dimensionality reduction matrix.
# We will follow the pipeline described in https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scanpy/scanpy_03_integration.html

# Import modules
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scanorama


# Read data
path_to_data = "data/TICAtlas.h5ad"
tica = sc.read_h5ad(path_to_data)


# Find highly variable genes
sc.pp.highly_variable_genes(tica, min_mean = 0.0125, max_mean = 3, min_disp = 0.5, batch_key = 'source')
var_genes = tica.var.highly_variable
var_genes = list(var_genes.index[var_genes])


# Run Scanorama
batches = tica.obs['source'].cat.categories.tolist()
alldata = {}
for batch in batches:
    alldata[batch] = tica[tica.obs['source'] == batch, var_genes]
tica_list = list(alldata.values())
scanorama.integrate_scanpy(tica_list, dimred = 50)
tica_list[0].obsm['X_scanorama'].shape
scanorama_int = [ad.obsm['X_scanorama'] for ad in tica_list]
cell_barcodes = []
for ad in tica_list:
    cell_barcodes.extend(ad.obs_names)
tica = tica[cell_barcodes, :]
all_s = np.concatenate(scanorama_int)
print(all_s.shape)
tica.obsm["Scanorama"] = all_s
#cell_barcodes = []
#for ad in tica_list:
#    cell_barcodes.extend(ad.obs_names)
#tica = tica[cell_barcodes, :]


# Run UMAP
sc.pp.neighbors(tica, n_pcs = 50, use_rep = "Scanorama")
sc.tl.umap(tica)


# Save
tica.__dict__['_raw'].__dict__['_var'] = tica.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
tica.write_h5ad("results/TICAtlas_scanorama3.h5ad")
dimred = tica.obsm["Scanorama"]
dimred_df = pd.DataFrame(dimred)
dimred_df.columns = ["Scanorama_{}".format(str(x + 1)) for x in dimred_df.columns]
dimred_df.insert(0, "cell_barcode", tica.obs_names, True)
dimred_df.to_csv("results/Scanorama_corrected_pca_coordinates3.csv", index = False)
umap_df = pd.DataFrame(tica.obsm["X_umap"])
umap_df.columns = ["UMAP1", "UMAP2"]
umap_df.insert(0, "cell_barcode", tica.obs_names, True)
umap_df.to_csv("results/Scanorama_corrected_umap_coordinates3.csv", index = False)
