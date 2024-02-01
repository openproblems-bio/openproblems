import anndata as ad
import logging
import numpy as np
from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD
from sklearn.gaussian_process.kernels import RBF
from sklearn.kernel_ridge import KernelRidge
logging.basicConfig(level=logging.INFO)

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

print('Reading input files', flush=True)
input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])
input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])

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

logging.info('Determine parameters by the modalities')
mod1_type = input_train_mod1.uns["modality"]
mod1_type = mod1_type.upper()
mod2_type = input_train_mod2.uns["modality"]
mod2_type = mod2_type.upper()
n_comp_dict = {
                ("GEX", "ADT"): (300, 70, 10, 0.2),
                ("ADT", "GEX"): (None, 50, 10, 0.2),
                ("GEX", "ATAC"): (1000, 50, 10, 0.1),
                ("ATAC", "GEX"): (100, 70, 10, 0.1)
              }
logging.info(f"{mod1_type}, {mod2_type}")
n_mod1, n_mod2, scale, alpha = n_comp_dict[(mod1_type, mod2_type)]
logging.info(f"{n_mod1}, {n_mod2}, {scale}, {alpha}")

# Perform PCA on the input data
logging.info('Models using the Truncated SVD to reduce the dimension')

if n_mod1 is not None and n_mod1 < input_train.shape[1]:
    embedder_mod1 = TruncatedSVD(n_components=n_mod1)
    mod1_pca = embedder_mod1.fit_transform(input_train.layers["counts"]).astype(np.float32)
    train_matrix = mod1_pca[input_train.obs['group'] == 'train']
    test_matrix = mod1_pca[input_train.obs['group'] == 'test']
else:
    train_matrix = input_train_mod1.to_df(layer="counts").values.astype(np.float32)
    test_matrix = input_test_mod1.to_df(layer="counts").values.astype(np.float32)
  
if n_mod2 is not None and n_mod2 < input_train_mod2.shape[1]:
    embedder_mod2 = TruncatedSVD(n_components=n_mod2)
    train_gs = embedder_mod2.fit_transform(input_train_mod2.layers["counts"]).astype(np.float32)
else:
    train_gs = input_train_mod2.to_df(layer="counts").values.astype(np.float32)

del input_train

logging.info('Running normalization ...')
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

logging.info('Running KRR model ...')
y_pred = np.zeros((pred_dimx, pred_dimy), dtype=np.float32)
np.random.seed(1000)

for _ in range(5):
  np.random.shuffle(batches)
  for batch in [batches[:batch_len//2], batches[batch_len//2:]]:
    # for passing the test
    if not batch:
      batch = [batches[0]]

    logging.info(batch)
    kernel = RBF(length_scale = scale)
    krr = KernelRidge(alpha=alpha, kernel=kernel)
    logging.info('Fitting KRR ... ')
    krr.fit(train_norm[feature_obs.batch.isin(batch)], train_gs[gs_obs.batch.isin(batch)])
    y_pred += (krr.predict(test_norm) @ embedder_mod2.components_)

np.clip(y_pred, a_min=0, a_max=None, out=y_pred)
if mod2_type == "ATAC":
    np.clip(y_pred, a_min=0, a_max=1, out=y_pred)

y_pred /= 10

# Store as sparse matrix to be efficient. 
# Note that this might require different classifiers/embedders before-hand. 
# Not every class is able to support such data structures.
y_pred = csr_matrix(y_pred)

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  layers = {
    'normalized': y_pred
  },
  obs = obs,
  var = var,
  uns = {
    'dataset_id': input_train_mod1.uns['dataset_id'],
    'method_id': meta['functionality_name']
  }
)
output.write_h5ad(par['output'], compression='gzip')
