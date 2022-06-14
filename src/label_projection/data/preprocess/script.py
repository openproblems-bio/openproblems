## VIASH START
par = {
    "input": "./test/data.h5ad",
    "method": 'batch',
    "output": "./test/preprocess.h5ad"
}
## VIASH END
import sys
sys.path.append(meta["resources_dir"])
import noise
import numpy as np
import scanpy as sc


def filter_genes_cells(adata):
    """Remove empty cells and genes."""
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=2)

    return adata


# TODO split the functions in different viash components
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

print(">> Load data")
adata = sc.read(par['input'])

print(">> Process data using {} method".format(par['method']))
filter_genes_cells(adata)
method_func = func_map[par['method']]
preprocessed_adata = method_func(adata)

print(">> Writing data")
preprocessed_adata.write(par['output'])
