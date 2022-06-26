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
majority = adata.obs.celltype[adata.obs.is_train].value_counts().index[0]
adata.obs["celltype_pred"] = np.nan
adata.obs.loc[~adata.obs.is_train, "celltype_pred"] = majority


print("Write output to file")
adata.uns["method_id"] = "majority_vote"
adata.write(par["output"], compression="gzip")
