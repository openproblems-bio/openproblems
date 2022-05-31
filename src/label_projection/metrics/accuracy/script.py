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
encoder = sklearn.preprocessing.LabelEncoder().fit(adata.obs["labels"])
test_data = adata[~adata.obs["is_train"]]

test_data.obs["labels"] = encoder.transform(test_data.obs["labels"])
test_data.obs["labels_pred"] = encoder.transform(test_data.obs["labels_pred"])

accuracy = np.mean(test_data.obs["labels"] == test_data.obs["labels_pred"])

print("Store metric value")
adata.uns["metric_id"] = "accuracy"
adata.uns["metric_value"] = accuracy

print("Writing adata to file")
adata.write(par['output'], compression="gzip")
