import subprocess
import pandas as pd
import tempfile
import scanpy as sc

# VIASH START
par = {
    "input_data": "https://zenodo.org/records/12785822/files/slidetag_human_cortex.tar.gz?download=1",
    "dataset_id": "zenodo_slidetags/human_cortex_slidetags",
    "dataset_name": "slidetag_human_cortex",
    "dataset_url": "https://www.nature.com/articles/s41586-023-06837-4",
    "dataset_summary": "Slide-tags enables single-nucleus barcoding for multimodal spatial genomics",
    "dataset_organism": "Homo sapiens",
    "dataset": "dataset.h5ad",
    "spot_filter_min_genes": 200,
    "gene_filter_min_spots": 50,
    "remove_mitochondrial": True
}
meta = {
    "name": "zenodo_slidetags"
}
# VIASH END

print(f"Downloading data", flush=True)
with tempfile.TemporaryDirectory() as tempdir:
    input_data = "input_data.tar.gz"
    dataset_name = par['dataset_name']
    epx_data = subprocess.run(
        ["wget", "-O", f"{tempdir}/{input_data}", par['input_data']], stderr=subprocess.STDOUT)
    extract_spatial = subprocess.run(
        ["tar", "-xzf", f"{tempdir}/{input_data}", "-C", tempdir, "--strip-components=1"], stderr=subprocess.STDOUT)

    # Read gene expression and create anndata object
    adata = sc.read_10x_mtx(path=tempdir)

    # Read spatial locations
    df = pd.read_csv(f"{tempdir}/spatial.csv", skiprows=1)
    df = df.set_index('TYPE')
    df.columns = ['spatial1', 'spatial2', 'cell_type']

    # add spatial locations to anndata object
    sel_cells = list(set(df.index) & set(adata.obs_names))

    df = df.loc[sel_cells, ]
    adata = adata[sel_cells, ]

    adata.obs = df
    adata.obsm['spatial'] = df[['spatial2', 'spatial1']].values

# Make variable names unique
adata.var_names_make_unique()

sc.pp.calculate_qc_metrics(adata, inplace=True)

print("Filtering spots or genes")
t0 = adata.shape
# remove cells with few counts
if par["spot_filter_min_counts"]:
    sc.pp.filter_cells(
        adata, min_counts=par["spot_filter_min_counts"], inplace=True)
# remove cells with few genes
if par["spot_filter_min_genes"]:
    sc.pp.filter_cells(
        adata, min_genes=par["spot_filter_min_genes"], inplace=True)
# remove genes that have few counts
if par["gene_filter_min_counts"]:
    sc.pp.filter_genes(
        adata, min_counts=par["gene_filter_min_counts"], inplace=True)
# remove genes that are found in few cells
if par["gene_filter_min_spots"]:
    sc.pp.filter_genes(
        adata, min_cells=par["gene_filter_min_spots"], inplace=True)
t1 = adata.shape
print(f"Removed {t0[0] - t1[0]} cells and {(t0[1] - t1[1])} genes.")

if par["remove_mitochondrial"]:
    print("Removing mitochondrial genes")
    non_mito_genes_list = [name for name in adata.var_names if not (
        name.startswith('MT-') or name.startswith('mt-'))]
    adata = adata[:, non_mito_genes_list]


# Rename .var columns
adata.var['feature_name'] = adata.var_names
adata.var.set_index(adata.var['gene_ids'], inplace=True)
adata.var.rename(columns={"gene_ids": "feature_id"}, inplace=True)

# Move counts to .layers
print("Add metadata to uns", flush=True)
adata.layers["counts"] = adata.X
adata.X = None

# Add metadata
print("Add metadata to uns", flush=True)
metadata_fields = ["dataset_id", "dataset_name", "dataset_url",
                   "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]
for key in metadata_fields:
    if key in par:
        print(f"Setting .uns['{key}']", flush=True)
        adata.uns[key] = par[key]

print("Writing adata to file", flush=True)
adata.write_h5ad(par["dataset"], compression="gzip")
