import numpy as np
import anndata as ad
import random

## VIASH START
par = {
    'input': 'resources_test/label_projection/pancreas/dataset_cpm.h5ad',
    'method': 'batch',
    'seed': None,
    'obs_batch': 'batch',
    'obs_label': 'celltype',
    'output_train': 'train.h5ad',
    'output_test': 'test.h5ad',
    'output_solution': 'solution.h5ad'
}
meta = {
    'resources_dir': 'src/label_projection/data_processing/censoring'
}
## VIASH END

if par["seed"]:
    print(f">> Setting seed to {par['seed']}")
    random.seed(par["seed"])

print(">> Load data")
adata = ad.read_h5ad(par["input"])

print("adata:", adata)

print(f">> Process data using {par['method']} method")

if par["method"] == "batch":
    test_batches = adata.obs[par["obs_batch"]].dtype.categories[[-3, -1]]
    is_test = [
        True if adata.obs[par["obs_batch"]][idx] in test_batches else False
        for idx in adata.obs_names
    ]
elif par["method"] == "random":
    train_ix = np.random.choice(adata.n_obs, round(adata.n_obs * 0.8), replace=False)
    is_test = [ not x in train_ix for x in range(0, adata.n_obs) ]

# create new anndata objects according to api spec
def subset_anndata(adata_sub, layers, obs, uns):
    return ad.AnnData(
        layers={key: adata_sub.layers[key] for key in layers},
        obs=adata_sub.obs[obs.values()].rename({v:n for n,v in obs.items()}, axis=1),
        var=adata.var.drop(adata.var.columns, axis=1),
        uns={key: adata_sub.uns[key] for key in uns}
    )
output_train = subset_anndata(
    adata_sub = adata[[not x for x in is_test]], 
    layers=["counts", "lognorm"], 
    obs={"label": par["obs_label"], "batch": par["obs_batch"]}, 
    uns=["raw_dataset_id", "dataset_id"]
)
output_test = subset_anndata(
    adata[is_test], 
    layers=["counts", "lognorm"], 
    obs={"batch": par["obs_batch"]}, # do NOT copy label to test obs!
    uns=["raw_dataset_id", "dataset_id"]
)
output_solution = subset_anndata(
    adata[is_test], 
    layers=["counts", "lognorm"],
    obs={"label": par["obs_label"], "batch": par["obs_batch"]},
    uns=["raw_dataset_id", "dataset_id"]
)

print(">> Writing data")
output_train.write_h5ad(par["output_train"])
output_test.write_h5ad(par["output_test"])
output_solution.write_h5ad(par["output_solution"])
