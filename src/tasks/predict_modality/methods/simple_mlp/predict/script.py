from glob import glob
import sys
import numpy as np
from scipy.sparse import csc_matrix
import anndata as ad
import torch
from torch.utils.data import TensorDataset,DataLoader

## VIASH START
par = {
    'input_train_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod1.h5ad',
    'input_train_mod2': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod2.h5ad',
    'input_test_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/test_mod1.h5ad',
    'input_model': 'output/model',
    'output': 'output/prediction'
}
meta = {
    'resources_dir': 'src/tasks/predict_modality/methods/simple_mlp',
    'cpus': 10
}
## VIASH END

resources_dir = f"{meta['resources_dir']}/resources"
sys.path.append(resources_dir)
from models import MLP
import utils

def _predict(model,dl):
  model = model.cuda()
  model.eval()
  yps = []
  for x in dl:
    with torch.no_grad():
      yp = model(x[0].cuda())
      yps.append(yp.detach().cpu().numpy())
  yp = np.vstack(yps)
  return yp


print('Load data', flush=True)
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])
input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])

# determine variables
mod_1 = input_test_mod1.uns['modality']
mod_2 = input_train_mod2.uns['modality']

task = f'{mod_1}2{mod_2}'

print('Load ymean', flush=True)
ymean_path = f"{par['input_model']}/{task}_ymean.npy"
ymean = np.load(ymean_path)

print('Start predict', flush=True)
if task == 'GEX2ATAC':
    y_pred = ymean*np.ones([input_test_mod1.n_obs, input_test_mod1.n_vars])
else:
    folds = [0, 1, 2]

    ymean = torch.from_numpy(ymean).float()
    yaml_path=f"{resources_dir}/yaml/mlp_{task}.yaml"
    config = utils.load_yaml(yaml_path)
    X = input_test_mod1.layers["normalized"].toarray()
    X = torch.from_numpy(X).float()
    
    te_ds = TensorDataset(X)
    
    yp = 0
    for fold in folds:
        # load_path = f"{par['input_model']}/{task}_fold_{fold}/version_0/checkpoints/*"
        load_path = f"{par['input_model']}/{task}_fold_{fold}/**.ckpt"
        print(load_path)
        ckpt = glob(load_path)[0]
        model_inf = MLP.load_from_checkpoint(
            ckpt,
            in_dim=X.shape[1],
            out_dim=input_test_mod1.n_vars,
            ymean=ymean,
            config=config
        )
        te_loader = DataLoader(
            te_ds,
            batch_size=config.batch_size,
            num_workers=0,
            shuffle=False,
            drop_last=False
        )
        yp = yp + _predict(model_inf, te_loader)

    y_pred = yp/len(folds)

y_pred = csc_matrix(y_pred)

adata = ad.AnnData(
    layers={"normalized": y_pred},
    shape=y_pred.shape,
    uns={
        'dataset_id': input_test_mod1.uns['dataset_id'],
        'method_id': meta['functionality_name'],
    },
)

print('Write data', flush=True)
adata.write_h5ad(par['output'], compression = "gzip")