## VIASH START
par = {
    'input': 'ouput.h5ad',
    'output': 'output.mv.h5ad'
}
## VIASH END
import numpy as np
import scanpy as sc


print("Load data")
adata = sc.read(par['input'])

print("Add celltype prediction")
celltype_distribution = adata.obs.celltype[adata.obs.is_train].value_counts()
celltype_distribution = celltype_distribution / celltype_distribution.sum()
adata.obs["celltype_pred"] = np.nan
adata.obs.loc[~adata.obs.is_train, "celltype_pred"] = np.random.choice(
    celltype_distribution.index,
    size=(~adata.obs.is_train).sum(),
    replace=True,
    p=celltype_distribution
)


print("Write output to file")
adata.uns["method_id"] = meta["functionality_name"]
adata.write(par["output"], compression="gzip")
