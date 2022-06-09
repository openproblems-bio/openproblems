## VIASH START
par = {
    'input': '../../../data/test_data_preprocessed.h5ad',
    'output': 'output.mlplogcpm.h5ad'
}
## VIASH END
resources_dir="../"

import sys
sys.path.append(resources_dir)
sys.path.append(meta['resources_dir'])
import scvi
import scanpy as sc
from tools import scanvi

###TODO add these as input
### n_hidden, n_latent, n_layers
n_latent, n_layers, n_hidden = (10, 1, 32)

print("Load input data")
adata = sc.read(par['input'])

train_kwargs = {
    "train_size": 0.9,
    "early_stopping": True,
}

# check parameters for test exists
par.get("max_epochs") and train_kwargs.update({"max_epochs": par['max_epochs']})
par.get("limit_train_batches") and train_kwargs.update({"limit_train_batches": par['limit_train_batches']})
par.get("limit_val_batches") and train_kwargs.update({"limit_val_batches": par['limit_val_batches']})

adata.obs["labels_pred"] = scanvi(adata, n_hidden, n_latent, n_layers,  **train_kwargs)
