## VIASH START
par = {
    "input": "../test_data.h5ad",
    "celltype_categories": [0, 3],
    "tech_categories": [0, -3, -2],
    "ouput": "./toy_data.h5ad"
}
## VIASH END


import scanpy as sc

def filter_genes_cells(adata):
    """Remove empty cells and genes."""
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=2)

    return adata


print(">> Load data")
adata = sc.read(par['input'])
adata = adata[:, :500].copy()
filter_genes_cells(adata)

print(">> Select indexes")
print(">> Selecting celltype_categories indexes {idx}".format(idx=par.get('celltype_categories')))
keep_celltypes = par.get('celltype_categories') and adata.obs["celltype"].dtype.categories[par['celltype_categories']]
print(">> Selected celltype_categories {}".format(keep_celltypes))
keep_celltype_idx = adata.obs["celltype"].isin(keep_celltypes)

print(">> Selecting tech_categories indexes {idx}".format(idx=par.get('tech_categories')))
keep_techs = par.get('tech_categories') and adata.obs["tech"].dtype.categories[par['tech_categories']]
print(">> Selected tech_categories {}".format(keep_techs))
keep_tech_idx = adata.obs["tech"].isin(keep_techs)

adata = adata[keep_tech_idx & keep_celltype_idx].copy()
sc.pp.subsample(adata, n_obs=500)
# Note: could also use 200-500 HVGs rather than 200 random genes
# Ensure there are no cells or genes with 0 counts
filter_genes_cells(adata)

print(">> Writing data")
adata.write(par['output'])
