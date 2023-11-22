import sys
import anndata as ad

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'hvgs': 2000,
    'obs_label': 'celltype',
    'obs_batch': 'batch',
    'subset_hvg': False,
    'output': 'output.h5ad'
}
meta = {
    "config": "target/nextflow/batch_integration/process_dataset/.config.vsh.yaml",
    "resources_dir": "src/common/helper_functions"
}
## VIASH END

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
        import scib
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

if par['subset_hvg']:
    print('Subsetting to HVG dimensions', flush=True)
    adata_with_hvg = adata_with_hvg[:, adata_with_hvg.var['hvg']].copy()

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
output_dataset = subset_anndata(adata_with_hvg, slot_info["output_dataset"])
output_solution = subset_anndata(adata_with_hvg, slot_info["output_solution"])

print('Writing adatas to file', flush=True)
output_dataset.write(par['output_dataset'], compression='gzip')
output_solution.write(par['output_solution'], compression='gzip')
