import anndata as ad
import numpy as np
import scvi

## VIASH START
par = {
    "input_train": "resources_test/label_projection/pancreas/train.h5ad",
    "input_test": "resources_test/label_projection/pancreas/test.h5ad",
    "output": "output.h5ad",
    "n_latent": 30,
    "n_layers": 2,
    "n_hidden": 128,
    "dropout_rate": 0.2,
    "max_epochs": 200,
}
meta = {"functionality_name": "scanvi_xgboost"}
## VIASH END

print("Reading input files", flush=True)
input_train = ad.read_h5ad(par["input_train"])
input_test = ad.read_h5ad(par["input_test"])
input_train.X = input_train.layers["counts"]
input_test.X = input_test.layers["counts"]

print("Train model", flush=True)
unlabeled_category = "Unknown"

scvi.model.SCVI.setup_anndata(input_train, batch_key="batch", labels_key="label")

# specific scArches parameters
arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=par["dropout_rate"],
    n_hidden=par["n_hidden"],
    n_layers=par["n_layers"],
    n_latent=par["n_latent"],
)
scvi_model = scvi.model.SCVI(input_train, **arches_params)
train_kwargs = dict(
    train_size=0.9,
    early_stopping=True,
)
scvi_model.train(**train_kwargs)
model = scvi.model.SCANVI.from_scvi_model(
    scvi_model, unlabeled_category=unlabeled_category
)
model.train(**train_kwargs)

query_model = scvi.model.SCANVI.load_query_data(input_test, model)
train_kwargs = dict(max_epochs=par["max_epochs"], early_stopping=True)
query_model.train(plan_kwargs=dict(weight_decay=0.0), **train_kwargs)

print("Generate predictions", flush=True)
input_test.obs["label"] = "Unknown"
input_test.obs["label_pred"] = query_model.predict(input_test)

print("Write output AnnData to file", flush=True)
input_test.uns["method_id"] = meta["functionality_name"]
input_test.write_h5ad(par["output"], compression="gzip")
