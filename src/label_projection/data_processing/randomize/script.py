## VIASH START
par = {
    "input": "../../../raw_data.h5ad",
    "method": 'batch',
    "output": "./test/preprocess.h5ad"
}
## VIASH END
import sys
sys.path.append(meta["resources_dir"])
import noise
import numpy as np
import scanpy as sc

print(">> Load data")
adata = sc.read(par['input'])

print(">> Process data using {} method".format(par['method']))
# Remove empty cells and genes.
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.filter_cells(adata, min_counts=2)

if par["method"] == "batch":
    test_batches = adata.obs["batch"].dtype.categories[[-3, -1]]
    adata.obs["is_train"] = [
        False if adata.obs["batch"][idx] in test_batches else True
        for idx in adata.obs_names
    ]
elif par["method"] == "random":
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )
elif par["method"] == "random_with_noise":
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )
    adata = noise.add_label_noise(adata, noise_prob=0.2)

print(">> Writing data")
adata.write(par['output'])
