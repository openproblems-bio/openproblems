print("Importing libraries")
import openproblems as op
import scanpy as sc
import scipy
import yaml

## VIASH START
par = {
    "id": "pancreas",
    "obs_celltype": "celltype",
    "obs_batch": "tech",
    "obs_tissue": "tissue",
    "layer_counts": "counts",
    "output": "test_data.h5ad",
}
## VIASH END

dataset_funs = {
    'allen_brain_atlas': op.data.allen_brain_atlas.load_mouse_brain_atlas,
    'cengen': op.data.cengen.load_cengen,
    'immune_cells': op.data.immune_cells.load_immune,
    'mouse_blood_olssen_labelled': op.data.mouse_blood_olssen_labelled.load_olsson_2016_mouse_blood,
    'mouse_hspc_nestorowa2016': op.data.mouse_hspc_nestorowa2016.load_mouse_hspc_nestorowa2016,
    'pancreas': op.data.pancreas.load_pancreas,
    # 'tabula_muris_senis': op.data.tabula_muris_senis.load_tabula_muris_senis,
    'tabula_muris_senis_droplet_lung': lambda : op.data.tabula_muris_senis.load_tabula_muris_senis(
        organ_list=["lung"], 
        method_list=["droplet"]
    ),
    'tenx_1k_pbmc': op.data.tenx.load_tenx_1k_pbmc,
    'tenx_5k_pbmc': op.data.tenx.load_tenx_5k_pbmc,
    'tnbc_wu2021': op.data.tnbc_wu2021.load_tnbc_data,
    # 'Wagner_2018_zebrafish_embryo_CRISPR': op.data.Wagner_2018_zebrafish_embryo_CRISPR.load_zebrafish_chd_tyr,
    'zebrafish': op.data.zebrafish.load_zebrafish
}

adata = dataset_funs[par['id']]()

print("Setting .uns metadata")
with open( "metadata.yaml") as md:
    metadata = yaml.safe_load(md)

for key, value in metadata[par["id"]].items():
    adata.uns[key] = value

print("Setting .obs['celltype']")
if par["obs_celltype"]:
    if par["obs_celltype"] in adata.obs:
        adata.obs["celltype"] = adata.obs[par["obs_celltype"]]
    else:
        print(f"Warning: key '{par['obs_celltype']}' could not be found in adata.obs.")

print("Setting .obs['batch']")
if par["obs_batch"]:
    if par["obs_batch"] in adata.obs:
        adata.obs["batch"] = adata.obs[par["obs_batch"]]
    else:
        print(f"Warning: key '{par['obs_batch']}' could not be found in adata.obs.")

print("Setting .obs['tissue']")
if par["obs_tissue"]:
    if par["obs_tissue"] in adata.obs:
        adata.obs["tissue"] = adata.obs[par["obs_tissue"]]
    else:
        print(f"Warning: key '{par['obs_tissue']}' could not be found in adata.obs.")

if par["layer_counts"] and par["layer_counts"] in adata.layers:
    print(f"  Temporarily moving .layers['{par['layer_counts']}'] to .X")
    adata.X = adata.layers[par["layer_counts"]]
    del adata.layers[par["layer_counts"]]

if par["sparse"] and not scipy.sparse.issparse(adata.X):
    print("  Make counts sparse")
    adata.X = scipy.sparse.csr_matrix(adata.X)

print("  Removing empty genes")
sc.pp.filter_genes(adata, min_cells=1)

print("  Removing empty cells")
sc.pp.filter_cells(adata, min_counts=2)

print(f"  Moving .X to .layers['counts']")
adata.layers["counts"] = adata.X
del adata.X

print("Writing adata to file")
adata.write_h5ad(par["output"], compression="gzip")
