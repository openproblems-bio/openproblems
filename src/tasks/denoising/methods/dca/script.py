import anndata as ad
from dca.api import dca

## VIASH START
par = {
    'input_train': 'resources_test/denoising/pancreas/train.h5ad',
    'output': 'output_dca.h5ad',
    'epochs': 300,
}
meta = {
    'functionality_name': 'dca',
}
## VIASH END

print("load input data", flush=True)
input_train = ad.read_h5ad(par['input_train'], backed="r")

print("Remove unneeded data", flush=True)
output = ad.AnnData(
    X=input_train.layers["counts"],
    obs=input_train.obs[[]],
    var=input_train.var[[]],
    uns={
        "dataset_id": input_train.uns["dataset_id"],
        "method_id": meta["functionality_name"]
    }
)

del input_train

print("Run DCA", flush=True)
dca(output, epochs=par["epochs"])

print("Move output to correct location", flush=True)
output.layers["denoised"] = output.X
del output.X

print("Writing data", flush=True)
output.write_h5ad(par["output"], compression="gzip")
