import anndata as ad
import numpy as np
from scipy.sparse import csc_matrix
from sklearn.decomposition import TruncatedSVD
from sklearn.gaussian_process.kernels import RBF
from sklearn.kernel_ridge import KernelRidge

## VIASH START
par = {
    'input_train_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/normal/train_mod1.h5ad',
    'input_train_mod2': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/normal/train_mod2.h5ad',
    'input_test_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/normal/test_mod1.h5ad',
    'output': 'output.h5ad', 
    'distance_method': 'minkowski', 
    'n_pcs': 50
}
meta = {
    'functionality_name': 'guanlab_dengkw_pm'
}
## VIASH END


## Removed PCA and normalization steps, as they arr already performed with the input data
print('Reading input files', flush=True)
input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])
input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])

batches = input_train_mod1.obs.batch.unique().tolist()
batch_len = len(batches)

# combine the train and test data
input_train = ad.concat(
    {"train": input_train_mod1, "test": input_test_mod1},
    axis=0,
    join="outer",
    label="group",
    fill_value=0,
    index_unique="-"
)

print('Determine parameters by the modalities', flush=True)
mod1_type = input_train_mod1.uns["modality"].upper()
mod2_type = input_train_mod2.uns["modality"].upper()
n_comp_dict = {
    ("GEX", "ADT"): (300, 70, 10, 0.2),
    ("ADT", "GEX"): (None, 50, 10, 0.2),
    ("GEX", "ATAC"): (1000, 50, 10, 0.1),
    ("ATAC", "GEX"): (100, 70, 10, 0.1)
}
print(f"{mod1_type}, {mod2_type}", flush=True)
n_mod1, n_mod2, scale, alpha = n_comp_dict[(mod1_type, mod2_type)]
print(f"{n_mod1}, {n_mod2}, {scale}, {alpha}", flush=True)

# Perform PCA on the input data
print('Models using the Truncated SVD to reduce the dimension', flush=True)

if n_mod1 is not None and n_mod1 < input_train.n_vars:
    embedder_mod1 = TruncatedSVD(n_components=n_mod1)
    mod1_pca = embedder_mod1.fit_transform(input_train.layers["normalized"]).astype(np.float32)
    train_matrix = mod1_pca[input_train.obs['group'] == 'train']
    test_matrix = mod1_pca[input_train.obs['group'] == 'test']
else:
    train_matrix = input_train_mod1.to_df(layer="normalized").values.astype(np.float32)
    test_matrix = input_test_mod1.to_df(layer="normalized").values.astype(np.float32)
  
if n_mod2 is not None and n_mod2 < input_train_mod2.n_vars:
    embedder_mod2 = TruncatedSVD(n_components=n_mod2)
    train_gs = embedder_mod2.fit_transform(input_train_mod2.layers["normalized"]).astype(np.float32)
else:
    train_gs = input_train_mod2.to_df(layer="normalized").values.astype(np.float32)

del input_train

print('Running normalization ...', flush=True)
train_sd = np.std(train_matrix, axis=1).reshape(-1, 1)
train_sd[train_sd == 0] = 1
train_norm = (train_matrix - np.mean(train_matrix, axis=1).reshape(-1, 1)) / train_sd
train_norm = train_norm.astype(np.float32)
del train_matrix

test_sd = np.std(test_matrix, axis=1).reshape(-1, 1)
test_sd[test_sd == 0] = 1
test_norm = (test_matrix - np.mean(test_matrix, axis=1).reshape(-1, 1)) / test_sd
test_norm = test_norm.astype(np.float32)
del test_matrix

print('Running KRR model ...', flush=True)
if batch_len == 1:
    # just in case there is only one batch
    batch_subsets = [batches]
elif mod1_type == "ADT" or mod2_type == "ADT":
    # two fold consensus predictions
    batch_subsets = [
        batches[:batch_len//2],
        batches[batch_len//2:]
    ]
else:
    # leave-one-batch-out consensus predictions
    batch_subsets = [
        batches[:i] + batches[i+1:]
        for i in range(batch_len)
    ]

y_pred = np.zeros((input_test_mod1.n_obs, input_train_mod2.n_vars), dtype=np.float32)
for batch in batch_subsets:
    print(batch, flush=True)
    kernel = RBF(length_scale = scale)
    krr = KernelRidge(alpha=alpha, kernel=kernel)
    print('Fitting KRR ... ', flush=True)
    krr.fit(
        train_norm[input_train_mod1.obs.batch.isin(batch)], 
        train_gs[input_train_mod2.obs.batch.isin(batch)]
    )
    y_pred += (krr.predict(test_norm) @ embedder_mod2.components_)

np.clip(y_pred, a_min=0, a_max=None, out=y_pred)
y_pred /= len(batch_subsets)

# Store as sparse matrix to be efficient. 
# Note that this might require different classifiers/embedders before-hand. 
# Not every class is able to support such data structures.
## Changed from csr to csc matrix as this is more supported.
y_pred = csc_matrix(y_pred)

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  layers = { 'normalized': y_pred },
  obs = input_test_mod1.obs[[]],
  var = input_train_mod2.var[[]],
  uns = {
    'dataset_id': input_train_mod1.uns['dataset_id'],
    'method_id': meta['functionality_name']
  }
)
output.write_h5ad(par['output'], compression='gzip')
