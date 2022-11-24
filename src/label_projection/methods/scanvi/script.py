import anndata as ad
import scarches as sca

# followed procedure from here:
# https://scarches.readthedocs.io/en/latest/scanvi_surgery_pipeline.html

## VIASH START
par = {
    'input_train': 'resources_test/label_projection/pancreas/train.h5ad',
    'input_test': 'resources_test/label_projection/pancreas/test.h5ad',
    'output': 'output.h5ad',
    'hvg': True
}
meta = {
    'functionality_name': 'scanvi'
}
## VIASH END

print("Load input data", flush=True)
input_train_orig = ad.read_h5ad(par['input_train'])
input_test_orig = ad.read_h5ad(par['input_test'])

if par["hvg"]:
    print("Subsetting to HVG", flush=True)
    input_train = input_train_orig[:,input_train_orig.var['hvg']]
    input_test = input_test_orig[:,input_test_orig.var['hvg']]
else:
    input_train = input_train_orig
    input_test = input_test_orig

print("Concatenating train and test data", flush=True)
input_train.obs['is_test'] = False
input_test.obs['is_test'] = True
input_test.obs['label'] = "Unknown"
adata = ad.concat([input_train, input_test], merge = "same")

print("Create SCANVI model and train it on fully labelled reference dataset", flush=True)
sca.models.SCVI.setup_anndata(
    adata, 
    batch_key="batch", 
    labels_key="label",
    layer="normalized"
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
input_test_orig.obs["label_pred"] = preds[adata.obs['is_test'].values]

print("Write output to file", flush=True)
input_test_orig.uns["method_id"] = meta["functionality_name"]
input_test_orig.write_h5ad(par['output'], compression="gzip")

