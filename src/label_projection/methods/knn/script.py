import anndata as ad
import sklearn.neighbors

## VIASH START
par = {
    'input_train': 'resources_test/label_projection/pancreas/train.h5ad',
    'input_test': 'resources_test/label_projection/pancreas/test.h5ad',
    'output': 'output.h5ad',
    'layer_input': 'counts'
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])

print("Fit to train data")
classifier = sklearn.neighbors.KNeighborsClassifier()
classifier.fit(input_train.obsm["X_pca"], input_train.obs["label"].astype(str))

print("Predict on test data")
input_test.obs["label_pred"] = classifier.predict(input_test.obsm["X_pca"])

print("Write output to file")
input_test.uns["method_id"] = meta["functionality_name"]
input_test.write_h5ad(par['output'], compression="gzip")
