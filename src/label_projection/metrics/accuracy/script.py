import numpy as np
import sklearn.preprocessing
import anndata as ad

## VIASH START
par = {
    'input_prediction': 'resources_test/label_projection/pancreas/dataset_cpm_knn.h5ad',
    'input_solution': 'resources_test/label_projection/pancreas/dataset_cpm_solution.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'accuracy'
}
## VIASH END

print("Load data")
input_prediction = ad.read_h5ad(par['input_prediction'])
input_solution = ad.read_h5ad(par['input_solution'])

assert (input_prediction.obs_names == input_solution.obs_names).all()

print("Encode labels")
encoder = sklearn.preprocessing.LabelEncoder().fit(input_solution.obs["label"])
input_solution.obs["label"] = encoder.transform(input_solution.obs["label"])
input_prediction.obs["label_pred"] = encoder.transform(input_prediction.obs["label_pred"])

print("Compute prediction accuracy")
accuracy = np.mean(input_solution.obs["label"] == input_prediction.obs["label_pred"])

print("Store metric value")
input_prediction.uns["metric_id"] = meta["functionality_name"]
input_prediction.uns["metric_value"] = accuracy

print("Writing adata to file")
input_prediction.write_h5ad(par['output'], compression="gzip")
