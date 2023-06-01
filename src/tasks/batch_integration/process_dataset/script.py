import sys
import scib
import anndata as ad

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'hvgs': 2000,
    'output': 'output.h5ad'
}
meta = {}
## VIASH END

# Remove this after upgrading to Viash 0.7.5
sys.dont_write_bytecode = True

# import helper functions
sys.path.append(meta['resources_dir'])
from subset_anndata import read_config_slots_info, subset_anndata

print('Read input', flush=True)
input = ad.read_h5ad(par['input'])

def compute_batched_hvg(adata, n_hvgs):
    adata = adata.copy()
    adata.X = adata.layers['normalized'].copy()
    if n_hvgs > adata.n_vars or n_hvgs <= 0:
        hvg_list = adata.var_names.tolist()
    else:
        hvg_list = scib.pp.hvg_batch(
            adata,
            batch_key='batch',
            target_genes=n_hvgs,
            adataOut=False
        )
    adata.var['hvg'] = adata.var_names.isin(hvg_list)
    del adata.X
    return adata

print(f'Select {par["hvgs"]} highly variable genes', flush=True)
adata_with_hvg = compute_batched_hvg(input, n_hvgs=par['hvgs'])

print(">> Figuring out which data needs to be copied to which output file", flush=True)
# use par arguments to look for label and batch value in different slots
slot_mapping = {
    "obs": {
        "label": par["obs_label"],
        "batch": par["obs_batch"],
    }
}
slot_info = read_config_slots_info(meta["config"], slot_mapping)

print(">> Create output object", flush=True)
output = subset_anndata(adata_with_hvg, slot_info["output"])

print('Writing adatas to file', flush=True)
output.write(par['output'], compression='gzip')
