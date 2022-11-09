## VIASH START
par = {
    'input': '../../../resources/toy_preprocessed_data.h5ad',
    'n_latent': 10,
    'n_layers': 1,
    'n_hidden': 32,
    'max_epochs': 1,
    'limit_train_batches': 10,
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
from tools import scanvi

print("Load input data")
adata = sc.read(par['input'])

train_kwargs = {
    "train_size": 0.9,
    "early_stopping": True,
}

# check if parameters for test exists
par.get("max_epochs") and train_kwargs.update({"max_epochs": par['max_epochs']})
par.get("limit_train_batches") and train_kwargs.update({"limit_train_batches": par['limit_train_batches']})
par.get("limit_val_batches") and train_kwargs.update({"limit_val_batches": par['limit_val_batches']})

adata.obs["celltype_pred"] = scanvi(adata, par['n_hidden'], par['n_latent'], par['n_layers'],  **train_kwargs)
adata.uns["method_id"] = meta["functionality_name"]

print("Write data")
adata.write_h5ad(par['output'], compression="gzip")
