import anndata as ad
import scanpy as sc
import sklearn.metrics
import scprep

## VIASH START
par = {
    'input_test': 'resources_test/denoising/pancreas/test.h5ad',
    'input_denoised': 'resources_test/denoising/pancreas/magic.h5ad',
    'output': 'output_mse.h5ad'
}
meta = {
    'functionality_name': 'mse'
}
## VIASH END

print("Load data", flush=True)
input_denoised = ad.read_h5ad(par['input_denoised'], backed="r")
input_test = ad.read_h5ad(par['input_test'], backed="r")

test_data = ad.AnnData(X=input_test.layers["counts"], dtype="float")
denoised_data = ad.AnnData(X=input_denoised.layers["denoised"], dtype="float")

print("Normalize data", flush=True)

# scaling and transformation
target_sum = 10000

sc.pp.normalize_total(test_data, target_sum)
sc.pp.log1p(test_data)

sc.pp.normalize_total(denoised_data, target_sum)
sc.pp.log1p(denoised_data)

print("Compute mse value", flush=True)
error = sklearn.metrics.mean_squared_error(
    scprep.utils.toarray(test_data.X), scprep.utils.toarray(denoised_data.X)
)

print("Store mse value", flush=True)
output = ad.AnnData(
    uns={ key: val for key, val in input_test.uns.items() },
)

output.uns["method_id"] = input_denoised.uns["method_id"]
output.uns["metric_ids"] = meta['functionality_name']
output.uns["metric_values"] = error

print("Write adata to file", flush=True)
output.write_h5ad(par['output'], compression="gzip")

