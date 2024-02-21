import knn_smooth
import anndata as ad

## VIASH START
par = {
    'input_train': 'resources_test/denoising/pancreas/train.h5ad',
    'output': 'output_knn.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data", flush=True)
input_train = ad.read_h5ad(par["input_train"], backed="r")

print("Remove unneeded data", flush=True)
X = input_train.layers["counts"].astype(float).transpose().toarray()

# Create output AnnData for later use
output = ad.AnnData(
    obs=input_train.obs[[]],
    var=input_train.var[[]],
    uns={
        "dataset_id": input_train.uns["dataset_id"],
        "method_id": meta["functionality_name"]
    }
)

del input_train

print("Run KNN smoothing", flush=True)
X = knn_smooth.knn_smoothing(X, k=10).transpose()

print("Process data", flush=True)
output.layers["denoised"] = X

print("Writing data", flush=True)
output.write_h5ad(par["output"], compression="gzip")
