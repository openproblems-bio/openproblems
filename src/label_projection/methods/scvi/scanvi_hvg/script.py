## VIASH START
par = {
    'input': '../../../resources/toy_preprocessed_data.h5ad',
    'n_top_genes': 2000,
    'max_epochs': 1,
    'limit_train_batches': 10,
    'span': 0.8,
    'limit_val_batches': 10,
    'output': 'output.scviallgenes.h5ad'
}
## VIASH END
resources_dir="../"

import sys
sys.path.append(resources_dir)
sys.path.append(meta['resources_dir'])
import scvi
import scanpy as sc
from tools import scanvi, hvg

print("Load input data")
adata = sc.read(par['input'])

hvg_kwargs = {
    "flavor": "seurat_v3",
    "inplace": False,
    "n_top_genes": par['n_top_genes'],
    "batch_key": "batch",

}

# check parameters for test exists
par.get("span") and hvg_kwargs.update({"span": par['span']})

train_kwargs = {
    "train_size": 0.9,
    "early_stopping": True,
}

# check parameters for test exists
par.get("max_epochs") and train_kwargs.update({"max_epochs": par['max_epochs']})
par.get("limit_train_batches") and train_kwargs.update({"limit_train_batches": par['limit_train_batches']})
par.get("limit_val_batches") and train_kwargs.update({"limit_val_batches": par['limit_val_batches']})

hvg_df = hvg(adata, **hvg_kwargs)
bdata = adata[:, hvg_df.highly_variable].copy()
adata.obs["celltype_pred"] = scanvi(bdata, par['n_hidden'], par['n_latent'], par['n_layers'],  **train_kwargs)
adata.uns["method_id"] = "scanvi_hvg"

print("Write data")
adata.write(par['output'], compression="gzip")
