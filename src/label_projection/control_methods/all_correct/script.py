## VIASH START
par = {
    'input': '../../../../resources_test/label_projection/pancreas/toy_normalized_log_cpm_data.h5ad',
    'output': 'output.mv.h5ad'
}
## VIASH END
import numpy as np
import scanpy as sc


print("Load data")
adata = sc.read(par['input'])

print("Add celltype prediction")
adata.obs["celltype_pred"] = adata.obs.celltype


print("Write output to file")
adata.uns["method_id"] = meta["functionality_name"]
adata.write(par["output"], compression="gzip")
