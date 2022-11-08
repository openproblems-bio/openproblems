import scanpy as sc
### VIASH START
par = {
    "input": "resources_test/common/pancreas/dataset.h5ad",
    # "keep_celltype_categories": ["acinar", "beta"],
    # "keep_batch_categories": ["celseq", "inDrop4", "smarter"],
    "even": True,
    "ouput": "toy_data.h5ad"
}
### VIASH END

def filter_genes_cells(adata):
    """Remove empty cells and genes."""
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=2)
    return adata


print(">> Load data")
adata = sc.read(par['input'])

# copy counts to .X because otherwise filter_genes and filter_cells won't work
adata.X = adata.layers["counts"]

if par.get('even'):
    keep_batch_categories = adata.obs["batch"].unique()
    adata_out = None
    n_batch_obs_per_value = 500 // len(keep_batch_categories)
    for t in keep_batch_categories:
        batch_idx = adata.obs["batch"] == t
        adata_subset = adata[batch_idx].copy()
        sc.pp.subsample(adata_subset, n_obs=min(n_batch_obs_per_value, adata_subset.shape[0]))
        if adata_out is None:
            adata_out = adata_subset
        else:
            adata_out = adata_out.concatenate(adata_subset, batch_key="_obs_batch")
    adata_out.uns = adata.uns
    adata_out.varm = adata.varm
    adata_out.varp = adata.varp
    adata = adata_out[:, :500].copy()
else:
    adata = adata[:, :500].copy()

filter_genes_cells(adata)

if par.get('keep_celltype_categories') and par.get('keep_batch_categories'):
    print(">> Selecting celltype_categories {categories}".format(categories=par.get('keep_celltype_categories')))
    print(">> Selecting batch_categories {categories}".format(categories=par.get('keep_batch_categories')))
    keep_batch_idx = adata.obs["batch"].isin(par['keep_batch_categories'])
    keep_celltype_idx = adata.obs["celltype"].isin(par['keep_celltype_categories'])
    adata = adata[keep_celltype_idx & keep_batch_idx].copy()

# Note: could also use 200-500 HVGs rather than 200 random genes
# Ensure there are no cells or genes with 0 counts
sc.pp.subsample(adata, n_obs=min(500, adata.shape[0]))
filter_genes_cells(adata)
adata.uns["dataset_id"] = adata.uns["dataset_id"] + "_subsample"

# remove previously copied .X
del adata.X

print(">> Writing data")
adata.write(par['output'])
