
import anndata as ad
import scvi

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

print("Load input data")
input_train_orig = ad.read_h5ad(par['input_train'])
input_test_orig = ad.read_h5ad(par['input_test'])

print("Subsetting to HVG")
input_train = input_train_orig[:,input_train_orig.var['hvg']]
input_test = input_test_orig[:,input_test_orig.var['hvg']]

print("Concatenating train and test data")
input_train.obs['is_test'] = False
input_test.obs['is_test'] = True
input_test.obs['label'] = "Unknown"
adata = ad.concat([input_train, input_test], merge = "same")

print("Setting up adata object")
adata.obs["scanvi_label"] = adata.obs["label"].to_numpy()
scvi.model.SCVI.setup_anndata(
    adata, 
    batch_key="batch",
    labels_key="scanvi_label",
    layer="normalized"
)

print("Train SCVI model")
train_kwargs = dict(
    train_size=0.9,
    early_stopping=True,
)
scvi_model = scvi.model.SCVI(adata)
scvi_model.train(**train_kwargs)

print("Train SCANVI model")
model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
model.train(**train_kwargs)

print("Make predictions")
preds = model.predict(adata)
input_test_orig.obs["label_pred"] = preds[adata.obs['is_test'].values]

print("Write output to file")
input_test_orig.uns["method_id"] = meta["functionality_name"]
input_test_orig.write_h5ad(par['output'], compression="gzip")

