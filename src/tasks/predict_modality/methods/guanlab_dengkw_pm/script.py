import anndata as ad
import numpy as np
import gc
from scipy.sparse import csc_matrix
from sklearn.gaussian_process.kernels import RBF
from sklearn.kernel_ridge import KernelRidge

## VIASH START
par = {
  'input_train_mod1': 'resources_test/predict_modality/neurips2021_bmmc_cite/train_mod1.h5ad',
  'input_train_mod2': 'resources_test/predict_modality/neurips2021_bmmc_cite/train_mod2.h5ad',
  'input_test_mod1': 'resources_test/predict_modality/neurips2021_bmmc_cite/test_mod1.h5ad',
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

dataset_id = input_train_mod1.uns['dataset_id']

pred_dimx = input_test_mod1.shape[0]
pred_dimy = input_train_mod2.shape[1]

feature_obs = input_train_mod1.obs
gs_obs = input_train_mod2.obs

batches = input_train_mod1.obs.batch.unique().tolist()
batch_len = len(batches)

obs = input_test_mod1.obs
var = input_train_mod2.var
dataset_id = input_train_mod1.uns['dataset_id']

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

scale = 10
alpha = 0.1 if (mod1_type == "ATAC" or mod2_type == "ATAC") else 0.2

train_norm = input_train_mod1.to_df(layer="normalized").values.astype(np.float32)
test_norm = input_test_mod1.to_df(layer="normalized").values.astype(np.float32)

train_gs = input_train_mod2.to_df(layer="normalized").values.astype(np.float32)

del input_train_mod1
del input_test_mod1
del input_train_mod2
gc.collect()

print('Running KRR model ...', flush=True)
y_pred = np.zeros((pred_dimx, pred_dimy), dtype=np.float32)
np.random.seed(1000)

for _ in range(5):
  np.random.shuffle(batches)
  for batch in [batches[:batch_len//2], batches[batch_len//2:]]:
    # for passing the test
    if not batch:
      batch = [batches[0]]

    print(batch, flush=True)
    kernel = RBF(length_scale = scale)
    krr = KernelRidge(alpha=alpha, kernel=kernel)
    print('Fitting KRR ... ', flush=True)
    krr.fit(train_norm[feature_obs.batch.isin(batch)], train_gs[gs_obs.batch.isin(batch)])
    y_pred += krr.predict(test_norm)

np.clip(y_pred, a_min=0, a_max=None, out=y_pred)
if mod2_type == "ATAC":
    np.clip(y_pred, a_min=0, a_max=1, out=y_pred)

y_pred /= 10

# Store as sparse matrix to be efficient. 
# Note that this might require different classifiers/embedders before-hand. 
# Not every class is able to support such data structures.
## Changed from csr to csc matrix as this is more supported.
y_pred = csc_matrix(y_pred)

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  layers = {
    'normalized': y_pred
  },
  obs = obs,
  var = var,
  uns = {
    'dataset_id': dataset_id,
    'method_id': meta['functionality_name']
  }
)
output.write_h5ad(par['output'], compression='gzip')
