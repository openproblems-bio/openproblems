import anndata as ad
import scarches as sca
import pandas as pd

# followed procedure from here:
# https://scarches.readthedocs.io/en/latest/scanvi_surgery_pipeline.html

## VIASH START
par = {
    'input_train': 'resources_test/label_projection/pancreas/train.h5ad',
    'input_test': 'resources_test/label_projection/pancreas/test.h5ad',
    'output': 'output.h5ad',
    'num_hvg': 2000
}
meta = {
    'functionality_name': 'scanvi'
}
## VIASH END

print("Load input data", flush=True)
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])

if par["num_hvg"]:
    print("Subsetting to HVG", flush=True)
    hvg_idx = input_train.var['hvg_score'].to_numpy().argsort()[:par["num_hvg"]]
    input_train = input_train[:,hvg_idx]
    input_test = input_test[:,hvg_idx]

print("Concatenating train and test data", flush=True)
input_train.obs['is_test'] = False
input_test.obs['is_test'] = True
input_test.obs['label'] = "Unknown"
adata = ad.concat([input_train, input_test], merge = "same")
del input_train

print("Create SCANVI model and train it on fully labelled reference dataset", flush=True)
sca.models.SCVI.setup_anndata(
    adata, 
    batch_key="batch", 
    labels_key="label",
    layer="counts"
)

vae = sca.models.SCVI(
    adata,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)

print("Create the SCANVI model instance with ZINB loss", flush=True)
scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")

print("Train SCANVI model", flush=True)
scanvae.train()

print("Make predictions", flush=True)
preds = scanvae.predict(adata)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=pd.DataFrame(
        {"label_pred": preds[adata.obs['is_test'].values]},
        index=input_test.obs.index,
    ),
    var=input_test.var[[]],
    uns={
        "dataset_id": input_test.uns["dataset_id"],
        "normalization_id": input_test.uns["normalization_id"],
        "method_id": meta["functionality_name"],
    },
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
