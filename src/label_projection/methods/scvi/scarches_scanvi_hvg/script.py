## VIASH START
par = {
    'input': '../../../resources/toy_preprocessed_data.h5ad',
    'span': 0.8,
    'n_top_genes': 2000,
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
from tools import scanvi_scarches, hvg

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

model_train_kwargs = {
    "train_size": 0.9,
    "early_stopping": True,
}

query_model_train_kwargs = {
    "max_epochs": 200,
    "early_stopping": True,
}

# check parameters for test exists
par.get("max_epochs") and model_train_kwargs.update({"max_epochs": par['max_epochs']}) and query_model_train_kwargs.update({"max_epochs": par['max_epochs']})
par.get("limit_train_batches") and model_train_kwargs.update({"limit_train_batches": par['limit_train_batches']}) and query_model_train_kwargs.update({"limit_train_batches": par['limit_train_batches']})
par.get("limit_val_batches") and model_train_kwargs.update({"limit_val_batches": par['limit_val_batches']}) and query_model_train_kwargs.update({"limit_val_batches": par['limit_val_batches']})

hvg_df = hvg(adata, **hvg_kwargs)
bdata = adata[:, hvg_df.highly_variable].copy()
adata.obs["celltype_pred"] = scanvi_scarches(bdata, par['n_hidden'], par['n_latent'], par['n_layers'], {'model_train_kwargs': model_train_kwargs,
                                                                                                      'query_model_train_kwargs': query_model_train_kwargs})
adata.uns["method_id"] = meta["functionality_name"]

print("Write data")
adata.write(par['output'], compression="gzip")
