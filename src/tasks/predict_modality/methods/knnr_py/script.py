import anndata as ad
from scipy.sparse import csc_matrix
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import KNeighborsRegressor

## VIASH START
par = {
    'input_train_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/train_mod1.h5ad',
    'input_train_mod2': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/train_mod2.h5ad',
    'input_test_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/test_mod1.h5ad',
    'distance_method': 'minkowski',
    'output': 'output.h5ad',
    'n_pcs': 4,
    'n_neighbors': 5,
}
meta = { 'functionality_name': 'foo' }
## VIASH END

print('Reading `h5ad` files...', flush=True)
input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])
input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])

input_train = ad.concat(
    {"train": input_train_mod1, "test": input_test_mod1},
    axis=0,
    join="outer",
    label="group",
    fill_value=0,
    index_unique="-"
)

print('Performing dimensionality reduction on modality 1 values...', flush=True)
embedder = TruncatedSVD(n_components=par['n_pcs'])
X = embedder.fit_transform(input_train.layers["normalized"])

# split dimred back up
X_train = X[input_train.obs['group'] == 'train']
X_test = X[input_train.obs['group'] == 'test']
y_train = input_train_mod2.layers["normalized"].toarray()

assert len(X_train) + len(X_test) == len(X)

print('Running KNN regression...', flush=True)

reg = KNeighborsRegressor(
    n_neighbors=par['n_neighbors'],
    metric=par['distance_method']
)

reg.fit(X_train, y_train)
y_pred = reg.predict(X_test)

y_pred = csc_matrix(y_pred)

adata = ad.AnnData(
    layers={"normalized": y_pred},
    obs=input_test_mod1.obs,
    var=input_train_mod2.var,
    uns={
        'dataset_id': input_train_mod1.uns['dataset_id'],
        'method_id': meta["functionality_name"],
    },
)

print('Storing annotated data...', flush=True)
adata.write_h5ad(par['output'], compression = "gzip")
