## VIASH START
par = {
    'input': 'ouput.h5ad',
    'output': 'output.mv.h5ad'
}
## VIASH END
import numpy as np
import sklearn.preprocessing
import scanpy as sc


print("Load data")
adata = sc.read(par['input'])

print("Get prediction accuracy")
encoder = sklearn.preprocessing.LabelEncoder().fit(adata.obs["celltype"])
test_data = adata[~adata.obs["is_train"]]

test_data.obs["celltype"] = encoder.transform(test_data.obs["celltype"])
test_data.obs["celltype_pred"] = encoder.transform(test_data.obs["celltype_pred"])

accuracy = np.mean(test_data.obs["celltype"] == test_data.obs["celltype_pred"])

print("Store metric value")
adata.uns["metric_id"] = meta["functionality_name"]
adata.uns["metric_value"] = accuracy

print("Writing adata to file")
adata.write_h5ad(par['output'], compression="gzip")
