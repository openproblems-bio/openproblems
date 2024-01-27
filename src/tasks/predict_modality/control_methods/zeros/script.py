import anndata
from scipy.sparse import csc_matrix
import numpy as np

# VIASH START
par = {
    "input_train_mod1": "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/train_mod1.h5ad",
    "input_test_mod1": "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/test_mod1.h5ad",
    "input_train_mod2": "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/train_mod2.h5ad",
    "output": "output.h5ad",
}

meta = {
    "functionality_name": "foo"
}
# VIASH END

print("Reading h5ad files", flush=True)
ad_mod1_test = anndata.read_h5ad(par["input_test_mod1"])
ad_mod2 = anndata.read_h5ad(par["input_train_mod2"])

print("create output objects", flush=True)
prediction = csc_matrix((ad_mod1_test.n_obs, ad_mod2.n_vars), dtype = np.float32)

out = anndata.AnnData(
    layers={"normalized": prediction},
    shape=prediction.shape,
    obs=ad_mod1_test.obs,
    var=ad_mod2.var,
    uns={
        "dataset_id": ad_mod2.uns["dataset_id"],
        "method_id": meta["functionality_name"],
    }
)

print("write predictions to file", flush=True)
out.write_h5ad(par["output"], compression="gzip")
