import numpy as np
import sklearn.preprocessing
import anndata as ad

## VIASH START
par = {
    'input_prediction': 'resources_test/label_projection/pancreas/knn.h5ad',
    'input_solution': 'resources_test/label_projection/pancreas/solution.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'accuracy'
}
## VIASH END

print("Load data", flush=True)
input_prediction = ad.read_h5ad(par['input_prediction'])
input_solution = ad.read_h5ad(par['input_solution'])

assert (input_prediction.obs_names == input_solution.obs_names).all(), "obs_names not the same in prediction and solution inputs"

print("Encode labels", flush=True)
cats = list(input_solution.obs["label"].dtype.categories) + list(input_prediction.obs["label_pred"].dtype.categories)
encoder = sklearn.preprocessing.LabelEncoder().fit(cats)
input_solution.obs["label"] = encoder.transform(input_solution.obs["label"])
input_prediction.obs["label_pred"] = encoder.transform(input_prediction.obs["label_pred"])

print("Compute prediction accuracy", flush=True)
accuracy = np.mean(input_solution.obs["label"] == input_prediction.obs["label_pred"])

print("Store metric value", flush=True)
input_prediction.uns["metric_ids"] = "accuracy"
input_prediction.uns["metric_values"] = accuracy

print("Writing adata to file", flush=True)
input_prediction.write_h5ad(par['output'], compression="gzip")
