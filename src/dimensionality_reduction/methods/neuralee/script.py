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

print("Load input data", flush=True)
input = ad.read_h5ad(par['input'])

print('Add method ID', flush=True)
input.uns['method_id'] = meta['functionality_name']

if input.uns['normalization_id'] == 'counts':
    print('Select top 500 high variable genes', flush=True)
    # idx = input.var['hvg_score'].to_numpy().argsort()[-500:]
    # dataset = GeneExpressionDataset(input.layers['counts'][:, idx])
    dataset = GeneExpressionDataset(input.layers['counts'])
    dataset.log_shift()
    dataset.subsample_genes(500)
    dataset.standardscale()
elif input.uns['normalization_id'] == 'log_cpm':
    n_genes = 1000
    print('Select top 1000 high variable genes', flush=True)
    idx = input.var['hvg_score'].to_numpy().argsort()[-n_genes:]
    input = input[:, idx].copy()
    dataset = GeneExpressionDataset(input.layers['normalized'])


# 1000 cells as a batch to estimate the affinity matrix
dataset.affinity_split(N_small=min(1000, input.n_obs))
NEE = NeuralEE(dataset, d=2, device=torch.device("cpu"))
fine_tune_kwargs = dict(verbose=False)
fine_tune_kwargs["maxit"] = par['maxit']
fine_tune_kwargs["maxit"] = 10
res = NEE.fine_tune(**fine_tune_kwargs)

input.obsm["X_emb"] = res["X"].detach().cpu().numpy()

print('Copy data to new AnnData object', flush=True)
output = ad.AnnData(
    obs=input.obs[[]],
    obsm={"X_emb": input.obsm["X_emb"]},
    uns={key: input.uns[key] for key in ["dataset_id", "normalization_id", "method_id"]}
)

print("Write output to file", flush=True)
output.write_h5ad(par['output'], compression="gzip")