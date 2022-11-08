import numpy as np
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/label_projection/pancreas/dataset_subsampled_cpm.h5ad',
    'method': 'batch',
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

print(">> Load data")
adata = sc.read(par["input"])

print("adata:", adata)

print(f">> Process data using {par['method']} method")

if par["method"] == "batch":
    test_batches = adata.obs[par["obs_batch"]].dtype.categories[[-3, -1]]
    is_test = [
        True if adata.obs[par["obs_batch"]][idx] in test_batches else False
        for idx in adata.obs_names
    ]
elif par["method"] == "random":
    is_test = np.random.choice(
        [False, True], adata.shape[0], replace=True, p=[0.8, 0.2]
    )

# create new anndata objects according to api spec
def subset_anndata(adata_sub, layers, obs, uns):
    return sc.AnnData(
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
