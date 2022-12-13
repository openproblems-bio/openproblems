import anndata as ad
import scanpy as sc
import yaml
import torch
from neuralee.embedding import NeuralEE
from neuralee.dataset import GeneExpressionDataset

## VIASH START
par = {
    'input': 'resources_test/dimensionality_reduction/pancreas/train.h5ad',
    'output': 'reduced.h5ad',
    'no_pca': False,
}
meta = {
    'functionality_name': 'neuralee',
    'config': 'src/dimensionality_reduction/methods/neuralee/config.vsh.yaml'
}
## VIASH END

print("Load input data")
input = ad.read_h5ad(par['input'])

print('Add method and normalization ID')
with open(meta['config'], 'r') as config_file:
    config = yaml.safe_load(config_file)

input.uns['normalization_id'] = config['functionality']['info']['preferred_normalization']
input.uns['method_id'] = meta['functionality_name']

if input.uns['normalization_id'] == 'counts':
    print('Select top 500 high variable genes')
    idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:500]
    dataset = GeneExpressionDataset(input.layers['counts'][:, idx])
    dataset.log_shift()
    dataset.standardscale()
elif input.uns['normalization_id'] == 'log_cpm':
    print('Select top 1000 high variable genes')
    idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:1000]
    dataset = GeneExpressionDataset(input.layers['normalized'][:, idx])

# 1000 cells as a batch to estimate the affinity matrix
dataset.affinity_split(N_small=min(1000, input.n_obs))
NEE = NeuralEE(dataset, d=2, device=torch.device("cpu"))
fine_tune_kwargs = dict(verbose=False)
fine_tune_kwargs["maxit"] = par['maxit']
fine_tune_kwargs["maxit"] = 10
res = NEE.fine_tune(**fine_tune_kwargs)

input.obsm["X_emb"] = res["X"].detach().cpu().numpy()

print("Delete layers and var")
del input.layers
del input.var

print("Write output to file")
input.write_h5ad(par['output'], compression="gzip")