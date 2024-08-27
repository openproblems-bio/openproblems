import sys
import anndata as ad
from scvi.model import SCVI, SCANVI

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
    'n_hvg': 2000,
    'n_latent': 30,
    'n_hidden': 128,
    'n_layers': 2,
    'max_epochs_scvi': 20,
    'max_epochs_scanvi': 20
}
meta = {
    'functionality_name' : 'scanvi',
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    X='layers/counts',
    obs='obs',
    var='var',
    uns='uns'
)

if par["n_hvg"]:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = adata.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    adata = adata[:, idx].copy()

print("Processing data", flush=True)
SCVI.setup_anndata(adata, batch_key="batch")

print("Run scVI", flush=True)
model_kwargs = {
    key: par[key]
    for key in ["n_latent", "n_hidden", "n_layers"]
    if par[key] is not None
}

vae = SCVI(adata, **model_kwargs)

vae.train(max_epochs=par["max_epochs_scvi"], train_size=1.0)

print('Run SCANVI', flush=True)
scanvae = SCANVI.from_scvi_model(
    scvi_model=vae,
    labels_key="label",
    unlabeled_category="UnknownUnknown", # pick anything definitely not in a dataset
)
scanvae.train(max_epochs=par["max_epochs_scanvi"], train_size=1.0)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        "X_emb": scanvae.get_latent_representation(),
    },
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["functionality_name"],
    },
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
