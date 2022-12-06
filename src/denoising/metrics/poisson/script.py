import anndata as ad
import scprep
import numpy as np

## VIASH START
par = {
    'input_denoised': 'output_magic.h5ad',
    'input_test': 'output_test.h5ad',
    'output': 'output_poisson.h5ad'
}
meta = {
    'functionality_name': 'poisson'
}
## VIASH END

print("Load Data")
input_denoised = ad.read_h5ad(par['input_denoised'])
input_test = ad.read_h5ad(par['input_test'])

test_data = input_test.layers["counts"].toarray()
denoised_data = input_denoised.layers["denoised"].toarray()

print("Compute metric value")
# scaling
initial_sum = input_denoised.layers["counts"].sum()
target_sum = test_data.sum()
denoised_data = denoised_data * target_sum / initial_sum

# from molecular_cross_validation.mcv_sweep import poisson_nll_loss
# copied from: https://github.com/czbiohub/molecular-cross-validation/blob/master/src/molecular_cross_validation/mcv_sweep.py
def poisson_nll_loss(y_pred: np.ndarray, y_true: np.ndarray) -> float:
    return (y_pred - y_true * np.log(y_pred + 1e-6)).mean()

error = poisson_nll_loss(scprep.utils.toarray(test_data), denoised_data)

print("Store poisson value")
input_denoised.uns["metric_ids"] = meta['functionality_name']
input_denoised.uns["metric_values"] = error

print("Write adata to file")
input_denoised.write_h5ad(par['output'], compression="gzip")
