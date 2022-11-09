import anndata as ad
import sklearn.linear_model
import sklearn.pipeline
import sklearn.preprocessing
import sklearn.decomposition
import scipy.sparse

## VIASH START
par = {
    'input_train': 'resources_test/label_projection/pancreas/dataset_cpm_train.h5ad',
    'input_test': 'resources_test/label_projection/pancreas/dataset_cpm_test.h5ad',
    'output': 'output.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])

print("Set up classifier pipeline")
def pca_op(adata_train, adata_test, n_components=100):
    is_sparse = scipy.sparse.issparse(adata_train.X)

    min_components = min(
        [adata_train.shape[0], adata_test.shape[0], adata_train.shape[1]]
    )
    if is_sparse:
        min_components -= 1
    n_components = min([n_components, min_components])
    if is_sparse:
        pca_op = sklearn.decomposition.TruncatedSVD
    else:
        pca_op = sklearn.decomposition.PCA
    return pca_op(n_components=n_components)

classifier = sklearn.linear_model.LogisticRegression(
    max_iter=par["max_iter"]
)
pipeline = sklearn.pipeline.Pipeline(
    [
        ("pca", pca_op(input_train, input_test, n_components=100)),
        ("scaler", sklearn.preprocessing.StandardScaler(with_mean=True)),
        ("regression", classifier),
    ]
)

print("Fit to train data")
classifipipelineer.fit(input_train.layers["lognorm"], input_train.obs["label"].astype(str))

print("Predict on test data")
input_test.obs["label_pred"] = pipeline.predict(input_test.layers["lognorm"])

print("Write output to file")
input_test.uns["method_id"] = meta["functionality_name"]
input_test.write_h5ad(par['output'], compression="gzip")