## VIASH START
par = {
    "input": "./test/data.h5ad",
    "test": False,
    "method": 'batch',
    "output": "./test/preprocess.h5ad"
}
resources_dir = '../../utils/'
## VIASH END

import sys
sys.path.append(resources_dir)
import noise
import preprocess
import numpy as np
import scanpy as sc


def batch(adata):
    adata.obs["labels"] = adata.obs["celltype"]
    adata.obs["batch"] = adata.obs["tech"]

    # Assign training/test
    test_batches = adata.obs["batch"].dtype.categories[[-3, -1]]
    adata.obs["is_train"] = [
        False if adata.obs["batch"][idx] in test_batches else True
        for idx in adata.obs_names
    ]
    return adata


def random(adata):
    adata.obs["labels"] = adata.obs["celltype"]
    adata.obs["batch"] = adata.obs["tech"]

    # Assign training/test
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )

    return adata


def random_with_noise(adata):
    adata.obs["labels"] = adata.obs["celltype"]
    adata.obs["batch"] = adata.obs["tech"]

    # Assign trainin/test
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )

    # Inject label noise
    adata = noise.add_label_noise(adata, noise_prob=0.2)

    return adata


func_map = {'batch': batch,
            'random': random,
            'random_with_noise': random_with_noise}

method_func = func_map[par['method']]
adata = sc.read(par['input'])

if par['test']:
    adata = adata[:, :500].copy()
    preprocess.filter_genes_cells(adata)
    keep_celltypes = adata.obs["celltype"].dtype.categories[[0, 3]]
    keep_techs = adata.obs["tech"].dtype.categories[[0, -3, -2]]
    keep_tech_idx = adata.obs["tech"].isin(keep_techs)
    keep_celltype_idx = adata.obs["celltype"].isin(keep_celltypes)
    adata = adata[keep_tech_idx & keep_celltype_idx].copy()
    sc.pp.subsample(adata, n_obs=500)
    # Note: could also use 200-500 HVGs rather than 200 random genes
    # Ensure there are no cells or genes with 0 counts
    preprocess.filter_genes_cells(adata)
else:
    preprocess.filter_genes_cells(adata)


preprocessed_adata = method_func(adata)
preprocessed_adata.write(par['output'])
