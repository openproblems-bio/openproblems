import anndata as ad
from scvi.model import SCVI

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
    'n_hvg': 2000,
    'n_latent': 30,
    'n_hidden': 128,
    'n_layers': 2,
    'max_epochs': 400
}
meta = {
    'functionality_name' : 'scvi',
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

if par["n_hvg"]:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = adata.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    adata = adata[:, idx].copy()

print("Processing data", flush=True)
SCVI.setup_anndata(adata, layer="counts", batch_key="batch")

print("Run scVI", flush=True)
model_kwargs = {
    key: par[key]
    for key in ["n_latent", "n_hidden", "n_layers"]
    if par[key] is not None
}

vae = SCVI(adata, **model_kwargs)

vae.train(max_epochs=par["max_epochs"], train_size=1.0)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        "X_emb": vae.get_latent_representation(),
    },
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["functionality_name"],
    },
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
